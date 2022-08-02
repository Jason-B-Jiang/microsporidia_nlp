# -----------------------------------------------------------------------------
#
# Create labelled microsporidia species + host spans for NER training
#
# Jason Jiang - Created: 2022/08/01
#               Last edited: 2022/08/02
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

import pandas as pd
from typing import List, Tuple, Dict
import re
import spacy
from spacy.tokens import DocBin, Span
import random
from pathlib import Path

################################################################################

## Global variables

random.seed(42)  # for consistent train/test splits
TRAIN_SPLIT = 0.8
nlp = spacy.blank('en')  # load blank spaCy pipeline as we only want tokenizer

################################################################################

def main() -> None:
    train, test = get_labelled_train_test_data()
    convert(train, Path('train.spacy'))
    convert(test, Path('test.spacy'))

################################################################################

## Helper functions

def get_labelled_train_test_data() -> \
    Tuple[List[Tuple[str, Dict[str, List[Tuple[int, int, str]]]]]]:
    """Docstring goes here.
    """
    labelled_ents_df = pd.read_csv('../../../data/formatted_species_host_data.csv')
    labelled_ents_df[['hosts_in_text', 'microsp_in_text']] = \
        labelled_ents_df[['hosts_in_text', 'microsp_in_text']].fillna(value='')

    labelled_ents = []

    for _, row in labelled_ents_df.iterrows():
        labelled_ents.append(
            get_labelled_entities(row['microsp_in_text'],
                                  row['hosts_in_text'],
                                  row['title_abstract'])
            )

    random.shuffle(labelled_ents)
    train_idx = int(round(len(labelled_ents) * TRAIN_SPLIT))

    # return tuple of training, testing data
    return labelled_ents[:train_idx], labelled_ents[train_idx:]


def get_labelled_entities(microsp_in_text: str, hosts_in_text: str, text: str) -> \
    Tuple[str, Dict[str, List[Tuple[int, int, str]]]]:
    """Docstring goes here.
    """
    labelled_ents = remove_overlapping_entities(
        get_entity_spans(microsp_in_text, text, 'MICROSPORIDIA'),
        get_entity_spans(hosts_in_text, text, 'HOST')
    )

    return (text, {'entities': labelled_ents})


def get_entity_spans(ents_in_text: str, text: str, label: str) -> \
    List[Tuple[int, int, str]]:
    """Docstring goes here.
    """
    # convert list of entities to set then back to list to remove
    # duplicate entities (ex: shared hosts between species)
    ents_in_text = list(set([s.strip() for s in re.split(r'(; | \|\|?)', ents_in_text) if \
        '|' not in s and ';' not in s and s != '' and s != ' ']))

    spans = []
    for ent in ents_in_text:
        abbrev = get_abbreviated_species_name(ent)
        if (abbrev == ent):
            spans.extend(get_match_spans(ent, text, label))
        else:
            spans.extend(get_match_spans(ent, text, label) + \
                get_match_spans(abbrev, text, label))
    
    return spans


def get_abbreviated_species_name(species: str) -> str:
    """Docstring goes here.
    """
    species = species.split(' ')
    if len(species) == 1:
        # species name is a single word, no abbreviation possible
        return species[0]
    
    return ' '.join([s[0] + '.' for s in species[:len(species) - 1]]) + \
        ' ' + species[-1]


def get_match_spans(subtext: str, text: str, label: str) -> \
    List[Tuple[int, int, str]]:
    """Docstring goes here.
    """
    # find matches in text w/ re.iter or something
    # convert character match indices to document match indices
    # return list of tuples of each match
    doc = nlp(text)
    match_spans = []

    # find literal matches for subtext in text
    for match in re.finditer(re.escape(subtext), doc.text):
        start, end = match.span()
        span = doc.char_span(start, end)

        if span is not None:
            match_spans.append((span.start, span.end, label))
    
    return match_spans


def remove_overlapping_entities(ents_1: List[Tuple[int, int, str]],
    ents_2: List[Tuple[int, int, str]]) -> List[Tuple[int, int, str]]:
    """Docstring goes here.
    """
    ent_1_overlaps = [ent for ent in ents_1 if \
        any([is_overlapping_ent(ent, other) for other in ents_2])]
    ent_2_overlaps = [ent for ent in ents_2 if \
        any([is_overlapping_ent(ent, other) for other in ents_1])]

    [ents_1.remove(ent) for ent in ent_1_overlaps]
    [ents_2.remove(ent) for ent in ent_2_overlaps]

    return ents_1 + ents_2


def is_overlapping_ent(ent_1: Tuple[int, int, str], ent_2: Tuple[int, int, str]) ->\
    bool:
    """Docstring goes here.
    """
    if len(set(range(ent_1[0], ent_1[1] + 1)) & \
        set(range(ent_2[0], ent_2[1] + 1))) > 0:
        return True
    
    return False


def convert(data: Tuple[List[Tuple[str, Dict[str, List[Tuple[int, int, str]]]]]],
    output_path: Path) -> None:
    """Converts spaCy 2 training lists to spaCy 3 binary training data.
    Inspired by:
    https://github.com/explosion/projects/blob/v3/pipelines/ner_demo/scripts/convert.py"""
    docs = []
    for entry in data:
        doc = nlp(entry[0])
        doc.ents = [Span(doc, i, j, label=label) for i, j, label in entry[1]['entities']]
        docs.append(doc)

    db = DocBin(docs=docs)
    db.to_disk(output_path)

################################################################################

if __name__ == '__main__':
    main()