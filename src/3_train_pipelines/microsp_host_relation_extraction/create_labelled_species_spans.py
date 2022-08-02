# -----------------------------------------------------------------------------
#
# Create labelled microsporidia species + host spans for NER training
#
# Jason Jiang - Created: 2022/08/01
#               Last edited: 2022/08/01
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

import pandas as pd
from typing import List, Tuple, Dict
import json
import re

################################################################################

def get_labelled_entity_spans() -> List[Tuple[str, Dict[int, int, str]]]:
    """Docstring goes here.
    """
    labelled_ents_df = pd.read_csv('../../../data/formatted_species_host_data.csv')
    labelled_ents = []

    for _, row in labelled_ents_df[1:5].iterrows():
        labelled_ents.append(
            (row['title_abstract'],
                {'entities': get_spans(row['host_spans'], 'HOST_SPECIES') + \
                    get_spans(row['microsp_spans'], 'MICROSPORIDIA_SPECIES')})
        )

    # write this list of labelled entity spans for each document as a json file
    # in the same directory as this script
    json.dumps(labelled_ents)


def get_spans(entity_spans: str, label: str) -> List[Tuple[int, int, str]]:
    """Docstring goes here.
    """
    spans = \
        [s.split('-') for s in re.split(r'(; |\|\| )', entity_spans) \
            if s not in ('; ', '|| ', '', ' ')]

    return [(int(i), int(j), label) for i, j in spans]