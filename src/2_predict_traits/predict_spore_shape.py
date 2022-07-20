# -----------------------------------------------------------------------------
#
# Predict microsporidia spore shape
#
# Jason Jiang - Created: 2022/07/13
#               Last edited: 2022/07/19
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

import spacy
from spacy.matcher import PhraseMatcher, Matcher
import pandas as pd
from pathlib import Path
from typing import List, Tuple

################################################################################

## Global variables

# Note: a lot of these variables are reused from predict_spore_nucleus_count.py

nlp = spacy.load('en_core_web_md')

# matcher for microsporidia shape descriptors
# these are terms that are commonly used in Microsporidia literature that I'm
# calling off the top of my head
shape_terms = ('oval', 'ovoid', 'round', 'pyriform', 'ovocylindrical',
               'spherical', 'ellipsoidal', 'ellipsoid', 'rod-shaped', 'rod')

shape_matcher = PhraseMatcher(nlp.vocab, attr='LOWER')
shape_matcher.add('shape', [nlp(term) for term in shape_terms])

# matcher for different types of mature microsporidia spores
# these spores are typically named [prefix]spore, ex: 'meiospores', 'macrospores'
# or, they're just called 'spores'
#
# make sure we don't match on 'exospore' or 'endospore', as these describe parts
# of microsporidia anatomy and not the spores themselves
#
# (same code from predict_spore_nucleus_count.py)
spore_matcher = Matcher(nlp.vocab)
spore_matcher.add('spore',  [[{'TEXT': {'REGEX': '(?<!(xo|do))[Ss]pore'}}]])

################################################################################

def main() -> None:
    microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')
    microsp_data = microsp_data.assign(
        pred_shape = lambda df: predict_spore_shape(df.title_abstract)
    )

    microsp_data.to_csv(Path('../../results/microsp_spore_shape_predictions.csv'))

################################################################################

## Helper functions

def predict_spore_shape(txt: str) -> str:
    """Predict shapes for each type of spore mentioned in some string, txt.
    Return a string in the form of "{shape} ({spore class name}); ..."

    Input:
        txt: abstract or full text for a Microsporidia paper
    """
    doc = nlp(txt)

    # get sentences from text with possible information about spores
    # (and thus spore shapes)
    # make a dictionary, with keys as spore sentences as values as
    # a dictionary of spore shapes + spore types in each sentence
    spore_sents = {sent: {'spore_types': spore_matcher(sent),
                          'spore_shapes': shape_matcher(sent)} \
                              for sent in doc.sents if spore_matcher(sent)}

    spore_shapes = []  # list of tuples, (spore name, spore shape)
    for sent in spore_sents:
        spore_spans = get_match_spans(spore_sents[sent]['spore_types'], sent)
        shape_spans = get_match_spans(spore_sents[sent]['spore_shapes'], sent)

        spore_shapes.extend(get_spore_type_shapes(spore_spans, shape_spans))
    
    return get_spore_shape_string(spore_shapes)


def get_match_spans(matches: List[Tuple[int]], sentence: spacy.tokens.span.Span) -> \
    List[Tuple[str, int, int]]:
    """Docstring goes here.
    """
    match_spans = []
    for match in matches:
        match_spans.append((sentence[match[1] : match[2]].text, match[1], match[2]))

    return match_spans


def get_spore_type_shapes(spore_spans: List[Tuple[str, int, int]],
                          shape_spans: List[Tuple[str, int, int]]) -> \
                              List[Tuple[str, str]]:
    """Docstring goes here.
    """
    shape_spans_copy = shape_spans.copy()  # to keep original list unmodified


def get_spore_shape_string(spore_shapes: List[Tuple[str, str]]) -> str:
    """Docstring goes here.
    """
    pass

################################################################################

if __name__ == '__main__':
    main()
    