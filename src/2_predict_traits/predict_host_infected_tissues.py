# -----------------------------------------------------------------------------
#
# Predict microsporidia sites of infection in hosts
#
# Jason Jiang - Created: 2022/06/02
#               Last edited: 2022/06/06
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

import scispacy
import spacy
from spacy.matcher import Matcher
from typing import List
from scispacy.linking import EntityLinker
import re
import pandas as pd
from pathlib import Path


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

## Cache texts as we process them with spaCy, to speed up code
CACHED_TEXT = {}

def get_cached_text(txt: str) -> spacy.tokens.doc.Doc:
    """Docstring goes here.
    """
    if txt not in CACHED_TEXT:
        CACHED_TEXT[txt] = nlp(txt)
    
    return CACHED_TEXT[txt]

################################################################################

## Cleaning up recorded sites of infection

# Manually parsed from inspecting 'Site of Infection' column
IRRELEVANT_PARENTHESIZED_INFO =\
    ['except', 'main', 'site', 'larva', 'adult', 'male', 'female', '\?', 'type',
    'at first', 'after', 'rarely', 'in ', 'spores', 'secondarily', 'between',
    'in advanced cases', 'prominent', 'close', '1', '2', 'chickens', 'of ',
    'does not', 'wall', 'low infection', 'all hosts', 'juvenile', 'embryo',
    'near', 'in one host individual', 'vertical', 'colonies', 'most organs',
    'similar tissues infected for both hosts', 'heavily infected', 'havily infected',
    'most infected' 'free spores', 'primarily', 'of Corethra (Savomyia) plumicornis',
    'anterior end', 'posterior end', 'most heavily infected', 'including',
    'of homo sapiens', 'of athymic mice']

EXCLUSIONS = ['in lamina propria through muscularis mucosae into submucosa of small intestine tract',
              'mainly in duodenum']

INFECTION_SITE_ENTITIES = ['ORGAN', 'TISSUE', 'ORGANISM_SUBDIVISION', \
    'MULTI_TISSUE_STRUCTURE', 'CELL', 'CELLULAR_COMPONENT']


def get_abbrev_species_name(host: str) -> str:
    """Docstring goes here.
    """
    if re.search("[A-Z]\. ", host):
        # host name already abbreviated, return as is
        return host
    
    host = host.split(' ')
    abbrevs = [s[0] + '.' for s in host[:len(host) - 1]]

    return ' '.join(abbrevs + [host[len(host) - 1]])


def get_host_names_and_abbrevs(hosts_natural, hosts_experimental) -> List[str]:
    """Docstring goes here.
    """
    # split hosts and join tgt in list
    # remove parenthesized text from each host
    # get abbreviated name for each host and add to host list
    # return the list
    hosts_natural = [h.split('(')[0].strip() for h in hosts_natural.split('; ')]
    hosts_experimental = [h.split('(')[0].strip() for h in hosts_experimental.split('; ')]
    all_hosts = hosts_natural + hosts_experimental
    all_hosts.extend([get_abbrev_species_name(h) for h in all_hosts])

    return(list(filter(lambda s: s != '', all_hosts)))


def should_remove_subsite(subsite: str, all_host_names: List[str]) -> bool:
    """Docstring goes here.
    """
    return any([term in subsite for term in IRRELEVANT_PARENTHESIZED_INFO + all_host_names]) and \
        not any([term in subsite for term in EXCLUSIONS])


def clean_recorded_infection_sites(sites: str, hosts_natural: str,
                                   hosts_experimental: str) -> str:
    """Format recorded sites of infection by ensuring parenthesized text
    are indicating subcellular/subtissue locations only, and parenthesized
    info are separated by commas.
    """
    sites = [s.strip() for s in sites.replace('(1;2)', '').split(';')]
    sites_formatted = []

    all_host_names = get_host_names_and_abbrevs(hosts_natural, hosts_experimental)

    for site in sites:
        # extract parenthesized info and split into subentries w/ ', '
        # strip away parenthesized text from main entry
        # remove subentries which have anything in irrelevant info as substring
        #   (or any entries in host species)
        # extend entries + subentries to sites_formatted
        if '(' in site:
            subsites = \
                [s.strip().strip(')') for s in re.search("(?<=\().+\)?", site).group(0).split(',')]
        else:
            subsites = []

        subsites = [s for s in subsites if not should_remove_subsite(s, all_host_names)]
        
        sites_formatted.append(site.split('(')[0].strip())
        sites_formatted.extend(subsites)
    
    return '; '.join(sites_formatted)


def get_umls_normalized_name(ent) -> str:
    """Docstring goes here.
    """
    if ent._.kb_ents:
        umls_ent = ent._.kb_ents[0]  # get first umls entry for entity
        return linker.kb.cui_to_entity[umls_ent[0]][1]

    return ent.text  # return entity text as is if no umls entry


def get_site_spans(doc: spacy.tokens.doc.Doc) -> List[spacy.tokens.span.Span]:
    """Docstring goes here
    """
    semicolon_posns = [i for i, tok in enumerate(doc) if tok.text == ';']
    curr_posn = 0
    sites = []

    for i in semicolon_posns:
        sites.append(doc[curr_posn : i])
        curr_posn = i + 1

    if semicolon_posns:
        return sites + [doc[semicolon_posns[-1] + 1:]]

    return [doc[:]]  # return doc as a span and 'as is' if no semicolons


def overlap_with_entity(site: spacy.tokens.span.Span, \
    ents: List[spacy.tokens.span.Span]) -> bool:
    """Docstring goes here.
    """
    for ent in ents:
        if ent.start >= site.start and ent.end <= site.end:
            return True
    
    return False


def normalize_recorded_infection_sites(sites_formatted: str) -> str:
    """Docstring goes here.
    """
    doc = get_cached_text(sites_formatted)
    sites_formatted = get_site_spans(doc)
    normalized_names = []
    ents_of_interest = [ent for ent in doc.ents if ent.label_ in INFECTION_SITE_ENTITIES]

    for ent in ents_of_interest:
        # only want to pick up broader level sites from organs or specific tissues
        normalized_names.append(get_umls_normalized_name(ent))

    # for recorded entries not tagged as entity by spaCy, add their names in
    # 'as is' for their normalized names
    for site in sites_formatted:
        if not overlap_with_entity(site, ents_of_interest):
            normalized_names.append(site.text)
    
    return '; '.join(list(filter(lambda x: x != '', normalized_names)))

################################################################################

## Predicting sites of infection from paper abstracts

INFECTION_LEMMAS = ['find', 'parasitize', 'infect', 'describe', 'localize', 'invade',
                    'appear', 'parasite']
infection_pattern = [{'POS': 'VERB', 'LEMMA': {'IN': INFECTION_LEMMAS}}]
matcher = Matcher(nlp.vocab)
matcher.add('infection_site', [infection_pattern])

MICROSPORIDIA_PARTS = ['polar tube', 'polar filament', 'sporoblast', 'spore'
                       'meront', 'meronts', 'sporoplasm', 'sporophorous vesicles',
                       'sporont', 'sporonts' 'polaroplast', 'anchoring disc',
                       'lamellar body', 'anchoring disk', 'endospore', 'exospore',
                       'posterior vacuole', 'sporoblasts', 'meiospores', 'meiospore',
                       'macrospore', 'macrospores', 'microspore', 'microspores',
                       'basal', 'spores', 'schizogony', 'sporogony']

def predict_infected_sites(txt: str) -> List[str]:
    """Docstring goes here.
    """
    doc = get_cached_text(txt)
    infection_sentences = [sent for sent in doc.sents if matcher(sent)]
    pred_sites = []

    for sent in infection_sentences:
        # extract all tissue/organ/etc entities from sentences possibly describing
        # microsporidia infection, and remove and entities corresponding to
        # microsporidia structures
        pred_sites.extend(
            [ent for ent in sent.ents if ent.label_ in INFECTION_SITE_ENTITIES and \
                not ent.text.lower() in MICROSPORIDIA_PARTS]
            )
    
    # Return UMLS normalized name (if available) for each detected entity
    return '; '.join(list(map(get_umls_normalized_name, pred_sites)))

################################################################################

microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Fill missing values in these columns with empty strings
microsp_data[['infection_site', 'hosts_natural', 'hosts_experimental']] = \
    microsp_data[['infection_site', 'hosts_natural', 'hosts_experimental']].fillna('')

# clean and normalize recorded infection sites
microsp_data['infection_site_formatted'] = microsp_data.apply(
    lambda df: clean_recorded_infection_sites(df.infection_site, df.hosts_natural,
                                              df.hosts_experimental),
    axis=1)

microsp_data['infection_site_normalized'] = microsp_data.apply(
    lambda df: normalize_recorded_infection_sites(df.infection_site_formatted),
    axis=1
)

# make predictions for infection sites from each paper title + abstract
microsp_data = microsp_data.assign(
    pred_infection_site = lambda df: df['title_abstract'].map(
        lambda txt: predict_infected_sites(txt), na_action='ignore'
    )
)

microsp_data[['species', 'title_abstract', 'infection_site', 'infection_site_formatted',
'infection_site_normalized', 'pred_infection_site']].to_csv(
    Path('../../results/microsp_infection_site_predictions.csv')
    )