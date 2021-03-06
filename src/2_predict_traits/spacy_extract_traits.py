# -----------------------------------------------------------------------------
#
# Extract Microsporidia traits from paper titles/abstracts with spaCy
#
# Jason Jiang - Created: 2022/05/13
#               Last edited: 2022/05/20
#
# Mideo Lab - Microsporidia text mining
#
# PLEASE IGNORE THIS SCRIPT, I'VE MOVED ALL WORK HERE TO OTHER SCRIPTS
#
# -----------------------------------------------------------------------------

import spacy
from spacy.matcher import Matcher, PhraseMatcher
from taxonerd import TaxoNERD
import pandas as pd
import re
from pathlib import Path

################################################################################

# Load in dataset of microsporidia traits + paper titles/abstracts, generated by
# src/taxonerd_entity_extraction/match_hosts_to_abstracts.R
microsp_data = pd.read_csv('../../data/abstracts_traits.csv')

# Remove rows w/out paper titles
microsp_data = microsp_data[microsp_data['first_paper_title'].notnull()]

# Replace missing paper abstracts with empty strings
microsp_data['abstract'] = microsp_data['abstract'].fillna('')

# Initialize spaCy English model + matchers, as global variables
NLP = spacy.load('en_core_web_md')
MATCHER = Matcher(NLP.vocab)
PHRASE_MATCHER = PhraseMatcher(NLP.vocab)

# Initialize Taxonerd species entity predictor
TAXONERD = TaxoNERD(model="en_ner_eco_biobert",  # use more accurate model
                    prefer_gpu=False,
                    with_abbrev=False)

################################################################################

## Data clean-up

# Replace tabs/newlines with single spaces, remove left + right trailing
# whitespaces and excess spaces (i.e: >1 space between words) in titles and
# abstracts

# Then, combine title + abstract into new variable, 'title_abstract_merge',
# where title and abstract are a single string, separated by '. '.
def clean_text(txt: str) -> str:
    """Replace tabs/newlines in txt with spaces, remove left + right
    trailing whitespaces and allow a max of 1 space between words.
    """
    return re.sub(' +', ' ', txt.replace('\t', ' ').replace('\n', ' ')).strip()

microsp_data = microsp_data.assign(
    abstract = lambda df: df['abstract'].map(
        lambda txt: clean_text(txt), na_action='ignore'
    ),
    title = lambda df: df['first_paper_title'].map(
        lambda txt: clean_text(txt), na_action='ignore'
    ),
    title_abstract_merge = 
        microsp_data['first_paper_title'] + '. ' + microsp_data['abstract']
)

################################################################################

## Functions for extracting microsporidia traits

# Extract location data
def extract_locality(txt: str) -> str:
    """Extract locations of discovery of the Microsporidia species.

    NOTE: this function currently doesn't distinguish between countries/cities/
          etc
    
    txt: String representing both paper title + abstract for a microsporidia

    return: semi-colon (';') separated string of spaCy predicted localities.
            Return an empty string if no locations are predicted by spaCy.
    """
    # Use set to remove redundant locations identified in text
    localities = \
        list(set([ent.text for ent in NLP(txt).ents if ent.label_ in ['GPE', 'LOC']]))

    return('; '.join(list(localities)))

microsp_data = microsp_data.assign(
    predicted_locality = lambda df: df['title_abstract_merge'].map(
        lambda txt: extract_locality(txt)
    )
)


# Extract host-parasite names
# Naive classifier for detecting parasite names, assuming the paper describes a
# new parasite species
SPECIES_REGEX = r'[A-Z][a-z]+ [a-z]+'
NEW_SPECIES_INDICATORS = r'{}'.format('(' + '|'.join(['n(\.|,)? ?g(en)?(\.|,)?,?.+?n(\.|,)? ?sp(\.|,)?,?',
                                         'g(en)?(\.|,)? ?n(\.|,)?,?.+sp(\.|,)? ?n(\.|,)?,?',
                                         'g(en)?(\.|,)? ?nov(\.|,)?,?.+sp(\.|,)? ?nov(\.|,)?,?',
                                         'n(\.|,)? ?sp(\.|,)?',
                                         'sp(\.|,)? ?n(\.|,)?',
                                         'sp(\.|,)? ?nov(\.|,)?',
                                         'nov(\.|,)? ?sp(\.|,)?'
                                         ]) + ')')
NEW_SPECIES_REGEX = r'{}'.format(
    f'({SPECIES_REGEX},? ?.?{NEW_SPECIES_INDICATORS}|[A-Z][a-z]+ {NEW_SPECIES_INDICATORS}|[A-Z][a-z]+ sp\\.)')

def extract_microsporidia_name(txt: str) -> str:
    """Extract microsporidia species name from paper title + abstract, using
    a naive regex search.
    
    txt: String representing both paper title + abstract for a microsporidia
    
    return: predicted Microsporidia species name from txt.
            Return empty string if no Microsporidia species name is found by
            regex.
    """
    pred = re.search(NEW_SPECIES_REGEX, txt)
    if pred:  # match for predicted microsporidia species name was found
        # extract ONLY the microsporidia species name from the matching substring
        return re.sub(' ?' + NEW_SPECIES_INDICATORS, '', pred.group())
    
    return ''

def predict_taxonerd_hosts(txt: str) -> str:
    """Predict and extract taxonomic/species names for hosts from text, using the
    TaxoNerd model.

    NOTE: a lot of non-specific family names (ex: Simuliidae) are also captured

    txt: Title + abstract for a microsporidia paper
    return: string of predicted species, separated by '; '
    """
    microsp, txt = txt.split(' ::: ')  # TODO: remove this
    taxonerd_pred = TAXONERD.find_in_text(txt)
    if not taxonerd_pred.empty:
        # get list of unique species entity predictions from text
        taxonerd_pred = list(set(taxonerd_pred['text'].tolist()))

        # exclude predicted parasite species from final list of hosts, and also
        # exclude entries that just say 'Microsporidia' or something similar
        return '; '.join(filter(lambda s: not microsp in s and not \
                                    re.search('[Mm]icrosp', s),
                         taxonerd_pred))

    return ''

microsp_data = microsp_data.assign(
    predicted_microsp = lambda df: df['title_abstract_merge'].map(
        lambda txt: extract_microsporidia_name(txt)
    )
)

microsp_data = microsp_data.assign(
    # I don't know how to pass multiple columns into functions for assign,
    # so this will have to do
    temp_var = microsp_data['predicted_microsp'] + ' ::: ' + \
        microsp_data['title_abstract_merge']
)

microsp_data = microsp_data.assign(
    predicted_hosts = lambda df: df['temp_var'].map(
        lambda txt: predict_taxonerd_hosts(txt)
    )
)


# Extract polar tube length
# Define regex patterns for extracting polar tube length
PT_REGEX = r'([Pp]olar ?)?([Ff]ilament|[Tt]ube|[Tt]ubule)'
MICRON_TERMS = ['??m', '?? m', '??', 'u ', 'mkm', 'mk m', 'mk m', '[Mm]icron[s]?',
                'pm', 'mum', 'm um', 'mu m', 'mu', 'mu.m', '.mu.m', 'fim', 'um',
                'u m ']
MICRON_REGEX = '(' + '|'.join(MICRON_TERMS) + ')'

# [^a-zA-Z\d\s:] = non-alphanumeric character
PT_MEASUREMENT_REGEX = \
    r'{}'.format(
        f'[1-9]+(\.|,)?[0-9]* ?{MICRON_REGEX}? ?(-|to|[^a-zA-Z\d\s:])+ ?([1-9]+(\.|,)?[0-9]* ?{MICRON_REGEX}?)?')

def extract_polar_tube_length(txt: str) -> str:
    """Docstring goes here
    """
    pred = re.search(
        r'{}'.format(  
            f"({PT_REGEX}.+{PT_MEASUREMENT_REGEX})"
            ),
        txt)

    if pred:
        return(re.search(PT_MEASUREMENT_REGEX, pred.group()).group())

    return ''

microsp_data = microsp_data.assign(
    predicted_polar_tube = lambda df: df['title_abstract_merge'].map(
        lambda txt: extract_polar_tube_length(txt), na_action='ignore'
    )
)


# Extract polar tube coils
COIL_REGEX = r'([Cc]oil[s]?|[Ss]pires[s]?|[Tt]urn[s]?|[Tt]wist[s]?)'

def extract_polar_tube_coils(txt: str) -> str:
    """Docstring goes here
    """
    pred = re.search(
        r'{}'.format(
            f"[1-9][0-9]* ?(-|to|[^a-zA-Z\d\s:]) ?[1-9]?[0-9]*.+?{COIL_REGEX}"
        ),
        txt)

    if pred:
        return(re.search(r'[1-9][0-9]* ?(-|to|[^a-zA-Z\d\s:]) ?[1-9]?[0-9]*',
        pred.group()).group())
    
    return ''


# Extract spore length/width data
pass

################################################################################

# Save modified microsp_dataframe with predicted microsporidia attributes
# and evaluate accuracy of predictions in R later
microsp_data.to_csv(Path('../../results/microsp_predictions.csv'))