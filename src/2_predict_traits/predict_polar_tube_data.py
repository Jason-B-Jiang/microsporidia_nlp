# -----------------------------------------------------------------------------
#
# Predict microsporidia polar tube coils + lengths from papers
#
# Jason Jiang - Created: 2022/05/20
#               Last edited: 2022/05/20
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

## Imports

import spacy
from spacy.matcher import PhraseMatcher, Matcher
import pandas as pd
import re
from pathlib import Path
from typing import Set, List

################################################################################

## Set up spaCy model

nlp = spacy.load('en_core_web_md')
matcher = Matcher(nlp.vocab)

################################################################################

# FINDING POLAR TUBE DATA
# Split text into sentences
# Look for sentences with mention of polar tube -> tube, tubule, filament
# For sentences that do match, extract coils and length data individually
# coils: number or range, preceded by mention of coils
# polar tube: number or range, followed by term for microns

# Just implement the polar tube search first

## Write function for detecting sentences with polar tube data
PT_REGEX = r'([Pp]olar ?)?([Ff]ilament|[Tt]ub(ul)?e)[s]?( |\.|,)'

def check_polar_tube_sentence(txt: str, nlp: spacy.lang.en.English) ->\
    Set[spacy.tokens.span.Span]:
    """Return set of sentences, each sentence represented as a spaCy span,
    mentioning polar tube in a text.
    """
    doc = nlp(txt)
    pt_sents = {sent for sent in doc.sents if re.search(PT_REGEX, sent.text)}
    return('; '.join(pt_sents))

################################################################################

## Define matcher from extracting polar tube coil data from polar tube sentences

PT_RANGE = r'(to|or|-|and)'  # let this appear one or more times

# 'l' of 'coil' sometimes gets messed up by ocr, so allow any lowercase char there
PT_COIL = r'([Cc]oi[a-z][s]?|[Ss]pire[s]?|[Tt]urn[s]?|[Tt]wist[s]?)'

# for excluding measurements in microns, as these aren't measuring number of polar
# tube coils
MICRON_TERMS = r'(μm|μ|mkm|um|mum)'

# NOTE: if a range of polar tube coils is found, the matcher will return both the
# full range match, and the partial range match
PT_COIL_DATA_PATTERN = [{'TEXT': '(', 'OP': '?'},
                        {'POS': 'NUM'},
                        {'TEXT': {'REGEX': PT_RANGE}, 'OP': '?'},
                        # {'TEXT': {'REGEX': PT_RANGE}, 'OP': '?'},
                        {'POS': 'NUM', 'OP': '?'},
                        # allow for bracket close to come before or after possible
                        # micron term
                        {'TEXT': ')', 'OP': '?'},
                        # {'TEXT': {'REGEX': MICRON_TERMS}, 'OP': '!'},
                        # allow for a max of 5 words between coil measurement and
                        # mention of coil term
                        # I tried using '*' operator, but didn't work in some cases
                        # {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '*'},
                        {'TEXT': ')', 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': '[a-z]+'}, 'OP': '?'},
                        {'TEXT': {'REGEX': PT_COIL}}]

matcher.add('pt_coil_data', [PT_COIL_DATA_PATTERN])

# Notes to self:
# incorporate coil mention coming before number of coils
# use match tuple locations in doc to see if matches overlap

################################################################################
