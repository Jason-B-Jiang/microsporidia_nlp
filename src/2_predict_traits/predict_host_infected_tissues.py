# -----------------------------------------------------------------------------
#
# Predict microsporidia sites of infection in hosts
#
# Jason Jiang - Created: 2022/06/02
#               Last edited: 2022/06/03
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

from asyncio import exceptions
import scispacy
import spacy
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

# Manually parsed from inspecting 'Site of Infection' column
irrelevant_parenthesized_info =\
    ['except', 'main', 'site', 'larva', 'adult', 'male', 'female', '\?', 'type',
    'at first', 'after', 'rarely', 'in ', 'spores', 'secondarily', 'between',
    'in advanced cases', 'prominent', 'close', '1', '2', 'chickens', 'of ',
    'does not', 'wall', 'low infection', 'all hosts', 'juvenile', 'embryo',
    'near', 'in one host individual', 'vertical', 'colonies', 'most organs',
    'similar tissues infected for both hosts', 'heavily infected', 'havily infected',
    'most infected' 'free spores', 'primarily', 'of Corethra (Savomyia) plumicornis',
    'anterior end', 'posterior end', 'most heavily infected', 'including',
    'of homo sapiens', 'of athymic mice']

exclusions = ['in lamina propria through muscularis mucosae into submucosa of small intestine tract',
              'mainly in duodenum']


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
    return any([term in subsite for term in irrelevant_parenthesized_info + all_host_names]) and \
        not any([term in subsite for term in exclusions])


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


def get_umls_normalized_name(ent) -> List[str]:
    for umls_ent in ent._.kb_ents:



def normalize_recorded_infection_sites(sites_formatted: str) -> list:
    """Docstring goes here.
    """
    doc = nlp(sites_formatted)
    normalized_names = []

    for ent in doc.ents:
        # only want to pick up broader level sites from organs or specific tissues
        if ent.label_ in ['ORGAN', 'TISSUE', 'ORGANISM_SUBDIVISION',
                          'MULTI_TISSUE_STRUCTURE']:
            normalized_names.extend(get_umls_normalized_name(ent))
    
    return '; '.join(normalized_names)

# for umls_ent in entity._.kb_ents:
# 	print(linker.kb.cui_to_entity[umls_ent[0]][1])

################################################################################

microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Fill missing values in these columns with empty strings
microsp_data[['infection_site', 'hosts_natural', 'hosts_experimental']] = \
    microsp_data[['infection_site', 'hosts_natural', 'hosts_experimental']].fillna('')

microsp_data['infection_site_formatted'] = microsp_data.apply(
    lambda x: clean_recorded_infection_sites(x.infection_site, x.hosts_natural,
                                             x.hosts_experimental),
    axis=1)

microsp_data[['species', 'hosts_natural', 'hosts_experimental', 
              'infection_site', 'infection_site_formatted']].to_csv(
    Path('./site_cleaning_check.csv')
    )