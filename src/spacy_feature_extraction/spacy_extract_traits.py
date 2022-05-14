# -----------------------------------------------------------------------------
#
# Extract Microsporidia traits from paper titles/abstracts with spaCy
#
# Jason Jiang - Created: 2022/05/13
#               Last edited: 2022/05/13
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

import spacy
from spacy.matcher import Matcher, PhraseMatcher
import pandas as pd
import re
from typing import List, Optional, Dict

################################################################################

# Load in dataset of microsporidia traits + paper titles/abstracts, generated by
# src/taxonerd_entity_extraction/match_hosts_to_abstracts.R
microsp_data = pd.read_csv('../../data/abstracts_traits.csv')

# Initialize spaCy English model + matchers, as global variables
NLP = spacy.load('en_core_web_md')
MATCHER = Matcher(NLP.vocab)
PHRASE_MATCHER = PhraseMatcher(NLP.vocab)

################################################################################

## Data clean-up

def clean_text(txt: str) -> str:
    """Replace tabs/newlines in string, txt with spaces, remove left + right
    trailing whitespaces and allow a max of 1 space between words.
    """
    return re.sub(' +', ' ', txt.replace('\t', ' ').replace('\n', ' ')).strip()

microsp_data = microsp_data.assign(
    abstract = lambda df: df['abstract'].map(
        lambda txt: clean_text(txt), na_action='ignore'
    ),
    title = lambda df: df['first_paper_title'].map(
        lambda txt: clean_text(txt), na_action='ignore'
    )
)
################################################################################

## Functions for extracting microsporidia traits

# Extract location data
def extract_locality(txts: Dict[str, Optional[str]]) -> str:
    """Extract locations of discovery of the Microsporidia species.
    Return an empty list if no locations are predicted by spaCy.

    NOTE: this function currently doesn't distinguish between countries/cities/
          etc
    
    txts: dictionary with 'title' and 'abstracts' as keys, and respective texts
          as values

    return: semi-colon (';') separated string of spaCy predicted localities
    """
    docs = {txt_type: NLP(txts[txt_type]) for txt_type in txts if txts[txt_type]}

    localities = set()
    for doc in docs:
        localities.update(
            list(set([ent.text for ent in docs[doc].ents if ent.label_ in ['GPE', 'LOC']]))
        )

    return('; '.join(list(localities)))


# Extract host-parasite names
# Naive classifier for detecting parasite names, assuming the paper describes a
# new parasite species

NEW_SPECIES_INDICATORS = ['n(\.|,)? ?sp(\.|,)?', 'sp(\.|,)? ?nov(\.|,)?',
                          'n(\.|,)? ?gen(\.|,)?.+?n(\.|,)? ?sp(\.|,)?',
                          'gen(\.|,)? ?nov(\.|,)?.+sp(\.|,)? ?nov(\.|,)?']
SPECIES_REGEX = '[A-Z][a-z]+ [a-z]+( [a-z]+)?'
NEW_SPECIES_REGEX = f'{SPECIES_REGEX} ?{NEW_SPECIES_INDICATORS_MERGED}'

# Extract polar tube coils

# Extract polar tube length
# Define regex patterns for extracting polar tube length
PT_REGEX = '([Pp]olar)? ([Ff]ilament|[Tt]ube|[Tt]ubule)'
MICRON_TERMS = ['μm', 'μ m', 'μ', 'u ', 'mkm', 'mk m', 'mk m', 'micron[s]?',
                'pm', 'mum', 'm um', 'mu m', 'mu', 'mu.m', '.mu.m', 'fim', 'um',
                'u m ']
MICRON_REGEX = '(' + '|'.join(MICRON_TERMS) + ')'
PT_LENGTH_REGEX = "[1-9]+\.?[0-9]* ?(-|to)? ?[1-9]*\.?[0-9]*"

def extract_polar_tube_length(txt: str) -> str:
    """Docstring goes here"""
    # This works in R but doesn't work as I'd expect in Python
    # I'm giving up with this for now
    return re.search(
        f'{PT_REGEX}.+{PT_LENGTH_REGEX} ?{MICRON_REGEX}|{PT_LENGTH_REGEX} ?{MICRON_REGEX}.+{PT_REGEX}',
        txt).group(0
        )