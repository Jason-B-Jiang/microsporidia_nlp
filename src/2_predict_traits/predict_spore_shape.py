# -----------------------------------------------------------------------------
#
# Predict microsporidia spore shape
#
# Jason Jiang - Created: 2022/07/13
#               Last edited: 2022/07/13
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

import spacy
from spacy.matcher import PhraseMatcher, Matcher
import pandas as pd
from pathlib import Path

################################################################################

## Global variables

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
    spore_sents = {sent: {'spore_types' = spore_matcher(sent), 'spore_shapes': shape_matcher(sent)} \
        for sent in doc.sents if spore_matcher(sent)}

    return

################################################################################

if __name__ == '__main__':
    main()
    