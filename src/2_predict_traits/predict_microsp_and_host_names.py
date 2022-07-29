# -----------------------------------------------------------------------------
#
# Predict microsporidia species names + hosts from paper titles + abstracts
#
# Jason Jiang - Created: 2022/05/19
#               Last edited: 2022/07/29
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

## Imports

import re
from webbrowser import get
import pandas as pd
from taxonerd import TaxoNERD
import spacy
from pathlib import Path
from typing import List, Tuple

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
 
# final regex pattern for microsporidia species in text
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

    microsp_data[['pred_microsporidia', 'pred_microsporidia_spans', 'pred_hosts',
        'pred_hosts_spans']] = \
            microsp_data.apply(
                lambda x: get_microsporidia_and_hosts(x['title_abstract']),
                axis=1,
                result_type='expand')

    # write modified microsp_data dataframe to results folder
    microsp_data[['title_abstract', 'species', 'pred_microsporidia',
                  'pred_microsporidia_spans', 'hosts_natural', 'hosts_experimental',
                  'pred_hosts', 'pred_hosts_spans']].to_csv(
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

    # return list of spans for the predicted taxons, from the original doc
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


def get_span_boundaries(spans: List[spacy.tokens.span.Span]) -> str:
    """Returns a semicolon-separated string of the start and stop boundaries of
    each span in a list of spans (wrt to their original document.)

    Ex: [hi, there] -> '0,1; 1,2'
    """
    return '; '.join([f"{s.start},{s.end}" for s in spans])


def get_microsporidia_and_hosts(txt: str) -> Tuple[str]:
    """Extract possible microsporidia + non-microsporidia species (putative
    microsporidia hosts) mentions from a text.

    Return a tuple of 4 strings:
        1) semi-colon separated list of Microsporidia species
        2) semi-colon separated list of span start + end for each predicted
           microsporidia species
        3) semi-colon separated list of non-Microsporidia species
        4) semi-colon separated list of span start + end for each predicted
           microsporidia species

    Or, return a tuple of 4 empty strings if a nan value is passed in as txt.
    
    Ex: ('Microsporidium sp. 1; Microsporidium sp. 2', '0,3; 3,6',
         'host 1; host 2', '0,2; 2,4')
    """
    if isinstance(txt, float):
        return ('', '', '', '')

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
    return ('; '.join([m.text for m in microsporidia]),
            get_span_boundaries(microsporidia),
            '; '.join([h.text for h in hosts]),
            get_span_boundaries(hosts))

################################################################################

if __name__ == '__main__':
    main()

################################################################################

## OLD CODE, KEEPING FOR REFERENCE

# Helper function for normalizing all predicted/recorded species names w/ GBIF,
# so we can directly compare predicted and recorded species names
#
# def get_gbif_normalized_name(sp: str) -> str:
#     """Docstring goes here.
#     """
#     gbif_hit = species.name_backbone(sp)
#     if gbif_hit:
#         if gbif_hit['synonym']:
#             # get canonical name for this synonymized species name
#             return gbif_hit['canonicalName']  # TODO - use speciesKey
        
#         return gbif_hit['canonicalName']
    
#     return sp  # no hits from GBIF search for species name, return as is