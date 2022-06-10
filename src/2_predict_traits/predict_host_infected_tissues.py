# -----------------------------------------------------------------------------
#
# Predict microsporidia sites of infection in hosts
#
# Jason Jiang - Created: 2022/06/02
#               Last edited: 2022/06/10
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

import scispacy
import spacy
from spacy.matcher import Matcher
from typing import List, Tuple
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


def get_overlapping_entity(site: spacy.tokens.span.Span, \
    ents: List[spacy.tokens.span.Span]) -> spacy.tokens.span.Span:
    """Docstring goes here.
    """
    for ent in ents:
        if ent.start >= site.start and ent.end <= site.end:
            return ent  # return the overlapping entity
    
    return site  # return site 'as is' if no overlapping entity found


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


def get_intersecting_umls_names(recorded_names: List[Tuple[str]], pred_names: List[Tuple[str]]) \
    -> List[str]:
    """Docstring goes here.
    """
    if not recorded_names:
        return []

    # get intersecting entries between recorded_names and pred_names
    return [s[0] for s in recorded_names if s[0] in [t[0] for t in pred_names]]


def resolve_normalized_names(pred_sites: dict, sites_formatted) -> Tuple[str, str]:
    """Docstring goes here.
    """
    doc = get_cached_text(sites_formatted)
    sites_formatted = get_site_spans(doc)
    ents_of_interest = [ent for ent in doc.ents if ent.label_ in INFECTION_SITE_ENTITIES]
    site_names = {}
    
    for i in range(len(sites_formatted)):
        # replace site with predicted entity if it can be associated with a predicted entity
        sites_formatted[i] = get_overlapping_entity(sites_formatted[i], ents_of_interest)

        site_names[sites_formatted[i]] = {'umls_names': sites_formatted[i]._.kb_ents,
                                          'needs_normalization': True,
                                          'canonical_name': sites_formatted[i].text}

        # look up each recorded site in all predicted sites and look for overlaps
        for pred in pred_sites:
            if sites_formatted[i].text.lower() in pred.text.lower():
                # for any predicted sites that have recorded site as a substring, set
                # its normalized canonical name to the text for the recorded site
                site_names[sites_formatted[i]]['needs_normalization'] = False
                pred_sites[pred]['needs_normalization'] = False
                pred_sites[pred]['canonical_name'] = sites_formatted[i].text
                continue

            intersecting_umls_names = get_intersecting_umls_names(
                site_names[sites_formatted[i]]['umls_names'],
                pred_sites[pred]['umls_names']
            )

            if intersecting_umls_names:
                canonical_name = linker.kb.cui_to_entity[intersecting_umls_names[0]][1]
                site_names[sites_formatted[i]]['canonical_name'] = canonical_name
                pred_sites[pred]['canonical_name'] = canonical_name

    return (
        '; '.join(list(set([site_names[site]['canonical_name'] for site in site_names]))),
        '; '.join(list(set([pred_sites[pred]['canonical_name'] for pred in pred_sites])))
    )


def predict_and_normalize_infection_sites(txt, sites_formatted) -> Tuple[str, str]:
    """Docstring goes here.
    """
    doc = get_cached_text(txt)
    infection_sentences = [sent for sent in doc.sents if matcher(sent)]
    pred_sites = {}  # keys = entities, values = {normalized_names: [], needs_normalization = True}

    for sent in infection_sentences:
        infection_ents = [ent for ent in sent.ents if ent.label_ in \
            INFECTION_SITE_ENTITIES and not ent.text.lower() in MICROSPORIDIA_PARTS]

        for ent in infection_ents:
            pred_sites[ent] = {'umls_names': ent._.kb_ents,
                               # if no umls entries for this entity, then don't need
                               # to normalize
                               'needs_normalization': len(ent._.kb_ents) > 0,
                               # set accepted normalized name as plain entity text for now
                               'canonical_name': ent.text}

    return resolve_normalized_names(pred_sites, sites_formatted)

################################################################################

microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Fill missing values in these columns with empty strings
microsp_data[['infection_site', 'hosts_natural', 'hosts_experimental', 'title_abstract']] = \
    microsp_data[['infection_site', 'hosts_natural', 'hosts_experimental', 'title_abstract']].fillna('')

# clean recorded infection sites
microsp_data['infection_site_formatted'] = microsp_data.apply(
    lambda df: clean_recorded_infection_sites(df.infection_site, df.hosts_natural,
                                              df.hosts_experimental),
    axis=1)

# normalize recorded infection sites with umls names (if possible) and get
# predicted infection sites
microsp_data[['infection_site_normalized', 'pred_infection_site']] = \
    [predict_and_normalize_infection_sites(txt, sites_formatted) for \
        txt, sites_formatted in \
            zip(microsp_data.title_abstract, microsp_data.infection_site_formatted)]

microsp_data[['species', 'num_papers', 'title_abstract', 'infection_site',
              'infection_site_formatted', 'infection_site_normalized',
              'pred_infection_site']].to_csv(
    Path('../../results/microsp_infection_site_predictions.csv')
    )