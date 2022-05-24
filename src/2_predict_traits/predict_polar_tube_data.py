# -----------------------------------------------------------------------------
#
# Predict microsporidia polar tube coils + lengths from papers
#
# Jason Jiang - Created: 2022/05/20
#               Last edited: 2022/05/23
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

## Imports

import spacy
from spacy.matcher import Matcher
import pandas as pd
import re
from pathlib import Path
from typing import List, Optional

################################################################################

## Set up spaCy model

nlp = spacy.load('en_core_web_md')
matcher = Matcher(nlp.vocab)

################################################################################

## Define matcher for extracting polar tube coil data from text

PT_RANGE = r'(to|or|-|and)'  # let this appear one or more times

# 'l' of 'coil' sometimes gets messed up by ocr, so allow any lowercase char there
PT_COIL = r'([Cc]oi[a-z](s|ed)?|[Ss]pire(s|d)?|[Tt]urn[s]?|[Tt]wist(s|ed)?)'

# for excluding measurements in microns, as these aren't measuring number of polar
# tube coils
MICRON_TERMS = r'(μm|μ|mkm|um|mum)'

# NOTE: if a range of polar tube coils is found, the matcher will return both the
# full range match, and the partial range match
PT_COIL_DATA_PATTERN = [{'TEXT': '(', 'OP': '?'},  # allow polar tube coils to be parenthesized
                        {'POS': 'NUM'},
                        {'TEXT': {'REGEX': PT_RANGE}, 'OP': '?'},
                        {'POS': 'NUM', 'OP': '?'},
                        # exclude measurements in microns, as these aren't
                        # counting number of polar tube coils
                        # {'TEXT': {'REGEX': MICRON_TERMS}, 'OP': '!'},
                        # allow for bracket close to come before or after possible
                        # micron term
                        {'TEXT': ')', 'OP': '?'},
                        # allow for a max of 5 words between coil measurement and
                        # mention of coil term
                        # I tried using '*' operator, but didn't work in some cases
                        # {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '*'},
                        # allow for parenthesized text to come after measurement
                        {'TEXT': '(', 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z0-9]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z0-9]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z0-9]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z0-9]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z0-9]+'}, 'OP': '?'},
                        {'TEXT': ')', 'OP': '?'},
                        {'TEXT': {'REGEX': PT_COIL}}]

matcher.add('pt_coil_data', [PT_COIL_DATA_PATTERN])

################################################################################

## Predicting polar tube coils

# Approach: extract sentences containing polar tube data from titles + abstracts,
# then extract polar tube coil data from these sentences

# Helper function for extracting polar tube coil data from sentences containg
# polar tube mentions
def extract_pt_coils(pt_sent: spacy.tokens.span.Span) -> str:
    """For a spaCy span representing a sentence containing polar tubule data,
    return a ' ||| ' separated string of polar tube coil measures extracted
    from the sentence.

    Input:
        pt_sents: spaCy span representing a sentence with polar tube data

    Return:
        ' ||| ' separated string of polar tube coil measures
        Empty string if no polar tube coil data is found
    """
    matches = matcher(pt_sent)

    if not matches:
        # no polar tube coil data was detected in sentence
        return ''

    curr_start = matches[0][1]
    curr_end = matches[0][2]
    for i in range(1, len(matches)):
        if matches[i][2] > curr_end:
            curr_start = matches[i][1]
            curr_end = matches[i][2]

        elif matches[i][2] == curr_end and matches[i][1] < curr_start:
            curr_start = matches[i][1]
            matches[i - 1] = None

        elif matches[i][1] > curr_start:
            matches[i] = None

    matches = list(filter(lambda s: s is not None, matches))
    matches_text = \
        list(filter(
            lambda s: not re.search(MICRON_TERMS, s),
            [pt_sent[start : end].text for id, start, end in matches]
        ))
    
    return ' ||| '.join(matches_text)
        

# Function for extracting polar tube coil measurements from sentences with
# polar tube coil data
PT_REGEX = r'([Pp]olar ?)?([Ff]ilament|[Tt]ub(ul)?e)[s]?( |\.|,)'

def predict_polar_tube_coils(txt: str) -> Optional[str]:
    """Predict polar tube coil measurements from a text, by getting sentences
    containing polar tube data and extracting coil data from these sentences.

    Input:
        txt: Title + abstract of a microsporidia paper.
    """
    doc = nlp(txt)
    pt_sents = [sent for sent in doc.sents if re.search(PT_REGEX, sent.text)]

    if not pt_sents:
        # None if no sentences containing polar tube data are detected
        return None

    pt_coil_preds = list(map(extract_pt_coils, pt_sents))
    pt_coil_preds = list(filter(lambda s: s != '', pt_coil_preds))

    if pt_coil_preds:
        return ' ||| '.join(pt_coil_preds)
    else:
        return None

################################################################################

## Predict microsporidia polar tube coils from paper titles + abstracts

# Load in formatted dataframe of microsporidia species data, from
# src/1_format_data/3_misc_cleanup.R
microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Exclude species with >1 papers describing them (for now)
microsp_data = microsp_data[microsp_data['num_papers'] < 2]

# Add column for predicted polar tube coil measurements
microsp_data = microsp_data.assign(
    pred_pt_coil = lambda df: df['title_abstract'].map(
        lambda txt: predict_polar_tube_coils(txt)
    )
)

################################################################################

# Write predictions to results folder
microsp_data[['species', 'title_abstract', 'pred_pt_coil', 'pt_coils_range',
              'pt_coils_avg']].to_csv(
    Path('../../results/microsp_pt_predictions.csv')
    )