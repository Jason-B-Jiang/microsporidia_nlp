# -----------------------------------------------------------------------------
#
# Predict microsporidia spore nucleus counts
#
# Jason Jiang - Created: 2022/07/12
#               Last edited: 2022/07/12
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

nucleus_terms = ('nucleus', 'nuclei', 'unikaryotic', 'unikariotic', 'unicaryotic',
                 'unicariotic', 'unikaryon', 'unikarion', 'unicaryon', 'unicarion',
                 'monokaryotic', 'monokariotic', 'monocaryotic', 'monocariotic',
                 'monokaryon', 'monokarion', 'monocaryon', 'monocarion',
                 'diplokaryotic', 'diplokariotic', 'diplocaryotic', 'diplocariotic',
                 'diplokaryon', 'diplokarion', 'diplocaryon', 'diplocarion',
                 'uninucleate', 'mononucleate', 'binucleate', 'uninuclear',
                 'binuclear')

nucleus_count = {'nucleus': 1, 'nuclei': 2, 'unikaryotic': 1, 'unikariotic': 1,
                 'unicaryotic': 1, 'unicariotic': 1, 'unikaryon': 1, 'unikarion': 1,
                 'unicaryon': 1, 'unicarion': 1, 'monokaryotic': 1, 'monokariotic': 1,
                 'monocaryotic': 1, 'monocariotic': 1, 'monokaryon': 1, 'monokarion': 1,
                 'monocaryon': 1, 'monocarion': 1, 'diplokaryotic': 2,
                 'diplokariotic': 2, 'diplocaryotic': 2, 'diplocariotic': 2,
                 'diplokaryon': 2, 'diplokarion': 2, 'diplocaryon': 2, 'diplocarion': 2,
                 'uninucleate': 1, 'mononucleate': 1, 'binucleate': 2, 'uninuclear': 1,
                 'binuclear': 2}

immature_spore_terms = ('sporoblast', 'sporont', 'meront', 'schizont')

# matcher for nucleus data terms
nucleus_matcher = PhraseMatcher(nlp.vocab, attr='LOWER')
nucleus_matcher.add('nucleus', [nlp(term) for term in nucleus_terms])

# matcher for sentences containing info about immature spore structures that
# may also have nucleus data (that we want to exclude)
immature_spore_matcher = PhraseMatcher(nlp.vocab, attr='LEMMA')
immature_spore_matcher.add('immature_spore', [nlp(term) for term in immature_spore_terms])

# matcher for different types of mature microsporidia spores
# these spores are typically named [prefix]spore, ex: 'meiospores', 'macrospores'
# or, they're just called 'spores'
#
# make sure we don't match on 'exospore' or 'endospore', as these describe parts
# of microsporidia anatomy and not the spores themselves
spore_matcher = Matcher(nlp.vocab)
spore_matcher.add('spore',  [[{'TEXT': {'REGEX': '(?<!(xo|do))[Ss]pore'}}]])

################################################################################

def main() -> None:
    """Makes spore nucleus predictions from Table S1 from Murareanu et al. and
    writes the resulting dataframe as a csv to the results folder.
    """
    microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')
    microsp_data['pred_nucleus'] = microsp_data.assign(
        lambda df: predict_spore_nucleus(df.title_abstract),
        axis=1
    )

    microsp_data.to_csv(Path('../../results/microsp_spore_nuclei_predictions.csv'))

################################################################################

## Helper functions

def predict_spore_nucleus(txt: str) -> str:
    """For a string, txt, return a semi-colon separated string of different spore
    types for Microsporidia species and the number of nuclei (1 or 2) for each
    spore type.
    
    Input:
        txt: abstract or full-text for a Microsporidia species paper.
    """
    doc = nlp(txt)
    nucleus_sents = [sent for sent in doc.sents if nucleus_matcher(sent) and \
        not immature_spore_matcher(sent)]

    # assume one sentence = nucleus data for one spore type
    # get first spore name from sentence, or just call 'spore' if no matches
    # get first mention of nucleus data and map to number of nuclei it describes
    # do this for each sentence, and return something like "c1: 1; c2: 2"
    spore_nucleus_info = []
    for sent in nucleus_sents:
        # TODO: maybe further split sentences by commas/semicolons, and
        # treat those as individual "sentences"?
        spore_name = spore_matcher(sent)
        if not spore_name:
            spore_name = 'normal spore'
        else:
            spore_name = sent[spore_name[0][1]].lemma_
            # TODO - clean this part up
            if spore_name == 'spore':
                spore_name = 'normal spore'
        
        nucleus_term = doc[nucleus_matcher(sent)[0][1]].text.lower()
        num_nuclei = str(nucleus_count[nucleus_term])

        spore_nucleus_info.append((spore_name, num_nuclei))

    return '; '.join([f"{num} ({name})" for name, num in spore_nucleus_info])

################################################################################

if __name__ == '__main__':
    main()