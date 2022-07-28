# -----------------------------------------------------------------------------
#
# Predict microsporidia species names + hosts from paper titles + abstracts
#
# Jason Jiang - Created: 2022/05/19
#               Last edited: 2022/07/28
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
from typing import Tuple, List

################################################################################

## Initialize language models

# Initialize Taxonerd species entity predictor
taxo = TaxoNERD(model="en_ner_eco_biobert",  # use more accurate model
                    prefer_gpu=False,
                    with_abbrev=True)

# Initialize spaCy medium pipeline
nlp = spacy.load('en_core_web_md')

################################################################################

## Define regex patterns for detecting new species names

SPECIES = r'[A-Z][a-z]+ [a-z]+( [a-z]+)?'

# all microsporidia genus names from NCBI
MICROSP_GENUSES = pd.read_csv(
    '../../data/microsp_genuses/microsporidia_genuses.tsv',
        sep='\t')['name'].tolist()

# conventional ways of indicating new species/genuses
NEW_SPECIES_INDICATORS = \
    [r'( [Nn](ov|OV)?(\.|,) [Gg](en|EN)?(\.|,))? (and |et )?[Nn](ov|OV)?(\.|,) [Ss][Pp](\.|,)',
    r'( [Gg](en|EN)?(\.|,) [Nn](ov|OV)?(\.|,))? (and |et )?[Ss][Pp](\.|,) [Nn](ov|OV)?(\.|,)']

NEW_SPECIES_INDICATORS =  r'{}'.format('(' + '|'.join(NEW_SPECIES_INDICATORS) + ')')
 
# final regex pattern for microspridia species in text
MICROSP_SPECIES = r'{}'.format(
    f"(({'|'.join(MICROSP_GENUSES)}) [A-Z]?[a-z]+|{SPECIES}(?=,?{NEW_SPECIES_INDICATORS}))"
    )

################################################################################

## Cache for storing texts already processed by spaCy

SPACY_CACHE = {}

################################################################################

def main() -> None:
    # load in microsporidia phenotype data + associated texts
    microsp_data = pd.read_csv(
        '../../data/manually_format_multi_species_papers.csv'
        )

    # add in column for Microsporidia species, as flagged by entity ruler
    # add in column for host species, which is any taxonomic entity detected by
    # TaxoNERD and is not a Microsporidia species
    # microsp_data = microsp_data.assign(
    #     pred_microsporidia = lambda x: x['title_abstract'].map(
    #         lambda s: get_microsporidia_from_text(s),
    #         na_action='ignore'
    #     )
    # )

    # write modified microsp_data dataframe to results folder
    microsp_data[['title_abstract', 'species', 'pred_microsporidia',
              'hosts_natural', 'hosts_experimental', 'pred_hosts']].to_csv(
    Path('../../results/microsp_and_host_predictions.csv')
    )

################################################################################

## Helper functions
# TODO - vectorize as much of this code as possible

def get_spacy_cached_text(txt: str) -> spacy.tokens.doc.Doc:
    if txt not in SPACY_CACHE:
        SPACY_CACHE[txt] = nlp(txt)
    
    return SPACY_CACHE[txt]


def get_microsporidia_from_doc(doc: spacy.tokens.doc.Doc) -> \
    List[spacy.tokens.span.Span]:
    """Return a list of spans corresponding to putative Microsporidia species
    mentioned in a spaCy document.
    
    Return an empty list if no Microsporidia species are found.
    """
    microsp_matches = []
    for match in re.finditer(MICROSP_SPECIES, doc.text):
        start, end = match.span()
        span = doc.char_span(start, end)

        if span is not None:
            microsp_matches.append(span)

    return microsp_matches


def get_taxonerd_taxons(doc: spacy.tokens.doc.Doc) -> List[spacy.tokens.span.Span]:
    """Return a list of spans corresponding to predicted taxons by TaxoNERD
    from document.
    """
    pred_taxons = taxo.find_in_text(doc.text)

    # filter out any rows in pred_taxons with 'microsporidium' or 'microsporidia'
    pred_taxons = pred_taxons[~pred_taxons['text'].str.contains("microsp")]

    return [doc.char_span(int(start), int(end)) for start, end in \
        [tuple(re.search('\d+ \d+', offset).group(0).split(' ')) \
            for offset in pred_taxons.offsets.tolist()]]



def is_overlapping_span(s1: spacy.tokens.span.Span, s2: spacy.tokens.span.Span) \
    -> bool:
    """Return True if any part of one span, s1, overlaps with the other, s2.
    """
    if len(set(range(s1.start, s1.end + 1)) & \
        set(range(s2.start, s2.end + 1))) > 0:
        return True
    
    return False


def format_predicted_microsp_and_hosts_str(microsp: List[spacy.tokens.span.Span],
    hosts: List[spacy.tokens.span.Span]) -> str:
    """
    """
    pass


def get_microsporidia_and_hosts(txt: str) -> str:
    """Extract possible microsporidia species mentions from a string, txt.
    Return a semi-colon formatted string of microsporidia species and their spans
    within the text (when tokenized by spaCy).
    
    Ex: Microsporidium sp. 1 (2, 5); Microsporidium sp. 2 (6, 9); ...
    """
    doc = get_spacy_cached_text(txt)
    microsporidia = get_microsporidia_from_doc(doc)
    pred_taxons = get_taxonerd_taxons(doc)
    
    # from hosts, add any predicted taxon that isn't a Microsporidia name
    hosts = []
    for taxon in pred_taxons:
        # TODO - add in dependency parsing to evaluate whether a predicted taxon
        # is a true Microsporidia host species
        if not any([is_overlapping_span(taxon, microsp) or taxon.text == microsp.text \
             for microsp in microsporidia]):
            hosts.append(taxon)
    
    # keep non-unique values for now, as we can use that to make weakly labelled
    # host/microsporidia data.
    return (microsporidia, hosts)
    

################################################################################

if __name__ == '__main__':
    main()

################################################################################

## OLD CODE, KEEPING FOR REFERENCE

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