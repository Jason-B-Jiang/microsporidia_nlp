# -----------------------------------------------------------------------------
#
# Predict microsporidia species names + hosts from paper titles + abstracts
#
# Jason Jiang - Created: 2022/05/19
#               Last edited: 2022/05/31
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

## Imports

import re
import pandas as pd
import taxonerd
from taxonerd import TaxoNERD
import spacy
from spacy.matcher import PhraseMatcher
from pathlib import Path
from pygbif import species

################################################################################

## Define regex patterns for detecting new species names

SPECIES = r'[A-Z](\.|[a-z]+) [A-Z]?[a-z]+( [a-z]+)?'

NEW_SPECIES_INDICATORS =\
    [r'[Nn](ov|OV)?(\.|,)? ?[Gg](en|EN)?(\.|,)?(.*?)?[Nn](ov|OV)?(\.|,)? ?[Ss][Pp](\.|,)?',
     r'[Gg](en|EN)?(\.|,)? ?[Nn](ov|OV)?(\.|,)?(.*?)?[Ss][Pp](\.|,)? ?[Nn](ov|OV)?(\.|,)?',
     r'[Nn](ov|OV)?(\.|,)? ?[Ss][Pp](\.|,)?',
     r'[Ss][Pp](\.|,)? ?[Nn](ov|OV)?(\.|,)?']

NEW_SPECIES_INDICATORS =  r'{}'.format('(' + '|'.join(NEW_SPECIES_INDICATORS) + ')')
 
NEW_SPECIES = r'{}'.format(  # final regex pattern for new species in text
    (f"{SPECIES},? {NEW_SPECIES_INDICATORS}|[A-Z][a-z]+ [Ss][Pp]\.?")
)

################################################################################

## Function for extracting new species names from titles/abstracts

def get_new_microsp_species(txt: str) -> str:
    """Predict new microsporidia species from a title + abstract.

    txt: title + abstract for a microsporidia species paper

    return: string of predicted species names, separated by '; '
    """
    new_sp = []
    matches = re.finditer(NEW_SPECIES, txt)

    for match in matches:
        start, end = match.span()
        # Remove new species indicator (ex: n. sp.) from predicted microsporidia
        # species, getting only the species name itself
        # Also remove excess trailing whitespaces
        new_sp.append(
            re.sub(NEW_SPECIES_INDICATORS, '', txt[start : end + 1]).strip()
            )
    
    # do list to set conversion to only get unique strings from list
    return '; '.join(list(set(new_sp)))

################################################################################

## Predict microsporidia names from titles + abstracts

# Load in formatted dataframe of microsporidia species data, from
# src/1_format_data/3_misc_cleanup.R
microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Exclude species with >1 papers describing them (for now)
microsp_data = microsp_data[microsp_data['num_papers'] < 2]

# Predict microsporidia species names
microsp_data = microsp_data.assign(
    pred_species = lambda df: df['title_abstract'].map(
        lambda txt: get_new_microsp_species(txt)
    )
)

################################################################################

## Predict hosts associated with microsporidia species

# NOTE: I'm taking a really naive approach, where I'm assuming any identified
#       species name that doesn't correspond to a microsporidia species is a host.
#       I'm also not distinguishing between which hosts belong to which
#       microsporidia species, if there's multiple microsporidia species in a
#       paper.

# Initialize Taxonerd species entity predictor
taxonerd = TaxoNERD(model="en_ner_eco_biobert",  # use more accurate model
                    prefer_gpu=False,
                    with_abbrev=False)

# Helper function for normalizing all predicted/recorded species names w/ GBIF,
# so we can directly compare predicted and recorded species names
def get_gbif_normalized_name(sp: str) -> str:
    """Docstring goes here.
    """
    gbif_hit = species.name_backbone(sp)
    if gbif_hit:
        if gbif_hit['synonym']:
            # get canonical name for this synonymized species name
            return gbif_hit['canonicalName']  # TODO - use speciesKey
        
        return gbif_hit['canonicalName']
    
    return sp  # no hits from GBIF search for species name, return as is


# Write helper function for extracting Taxonerd predictions from text
def predict_taxonerd_hosts(microsp: str, txt: str, taxonerd) -> str:
    """Predict and extract taxonomic/species names for hosts from text, using the
    TaxoNerd model.

    NOTE: a lot of non-specific family names (ex: Simuliidae) are also captured

    microsp: Predicted microsporidia species name(s) in paper
    txt: Title + abstract for a microsporidia paper
    return: string of predicted species, separated by '; '
    """
    taxonerd_pred = taxonerd.find_in_text(txt)

    if not taxonerd_pred.empty:
        # get list of unique species entity predictions from text
        taxonerd_pred = list(set(taxonerd_pred['text'].tolist()))

        # exclude predicted microsporidia species from final list of hosts, and
        # also exclude entries that just say 'Microsporidia' or something
        return '; '.join(filter(lambda s: not microsp in s and not \
                                    re.search('[Mm]icrosp', s),
                         taxonerd_pred))

    return ''

# Make predictions for microsporidia hosts
# probably the most atrocious code i've written
microsp_data['pred_hosts'] = \
    [predict_taxonerd_hosts(*a) for a in tuple(zip(microsp_data['pred_species'],
        microsp_data['title_abstract'], [taxonerd] * len(microsp_data)))]

################################################################################

## New approach for extracting microsporidia names: use microsporidia genuses
## from NCBI to find microsporidia species amongst Taxonerd predictions

microsp_genuses = pd.read_csv('../../data/microsp_genuses/microsporidia_genuses.tsv',
                              sep='\t')['name'].tolist()

nlp = spacy.load('en_core_web_md')

genus_matcher = PhraseMatcher(nlp.vocab)
genus_matcher.add('microsp_genus', [nlp(genus) for genus in microsp_genuses])

def get_microsp_species(taxonerd_species: str) -> str:
    """Docstring goes here.
    """
    taxonerd_species = [nlp(sp) for sp in taxonerd_species.split('; ')]
    microsp_sp = []

    for sp in taxonerd_species:
        if genus_matcher(sp):
            microsp_sp.append(sp.text)

    return '; '.join(microsp_sp)

microsp_data = microsp_data.assign(
    pred_species_2 = lambda df: df['pred_hosts'].map(
        lambda species: get_microsp_species(species), na_action='ignore'
    )
)

################################################################################

## Write resulting dataframe of predictions to results folder

microsp_data[['title_abstract', 'species', 'pred_species', 'pred_species_2',
              'hosts_natural', 'hosts_experimental', 'pred_hosts']].to_csv(
    Path('../../results/microsp_and_host_predictions.csv')
    )
