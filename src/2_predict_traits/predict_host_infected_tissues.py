# -----------------------------------------------------------------------------
#
# Predict microsporidia sites of infection in hosts
#
# Jason Jiang - Created: 2022/06/02
#               Last edited: 2022/06/02
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

import scispacy
import spacy
from typing import List
from scispacy.linking import EntityLinker

################################################################################

## Model initialization

# Use scispaCy en_ner_bionlp13cg_md NER model, as this has been trained to tag
# organs, tissues, etc
nlp = spacy.load("en_ner_bionlp13cg_md")

# Add entity linking from the Unified Medical Language System, for getting
# normalizing tissue/organ/etc names w/ canonical UMLS names
nlp.add_pipe("scispacy_linker",
             config={"resolve_abbreviations": True, "linker_name": "umls"})

linker = nlp.get_pipe("scispacy_linker")

################################################################################

def predict_infected_host_tissues(txt: str) -> List[str]:
    """Extract all 'ORGAN' and 'TISSUE' entities from a text, assuming these
    are microsporidia sites of infection in their hosts.
    
    NOTE: this only extracts these entities from texts, and doesn't associate
    them with any particular hosts.
    """
    doc = nlp(txt)
    return [ent for ent in doc.ents if ent.label_ in ['ORGAN', 'TISSUE']]


def normalize_entity_names(ents: List[spacy.tokens.span.Span]) -> List[str]:
    """Get UMLS normalized names for a list of named entities corresponding
    to tissues, organs, etc.
    """
    pass

# for umls_ent in entity._.kb_ents:
# 	print(linker.kb.cui_to_entity[umls_ent[0]][0])