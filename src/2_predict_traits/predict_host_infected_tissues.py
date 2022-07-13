# -----------------------------------------------------------------------------
#
# Predict microsporidia sites of infection in hosts
#
# Jason Jiang - Created: 2022/06/02
#               Last edited: 2022/07/13
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

import spacy
import scispacy
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
nlp_bio = spacy.load("en_ner_bionlp13cg_md")

# Alternatively: use large model to tag all entities, then link to UMLS
nlp_large = spacy.load("en_core_sci_lg")  # use large model for optimal entity tagging

# Add entity linking from the Unified Medical Language System, for getting
# normalizing tissue/organ/etc names w/ canonical UMLS names
#
# Only add entity linking for the en_core_sci_lg pipeline, due to memory
# constraints
nlp_large.add_pipe("scispacy_linker",
                   config={"resolve_abbreviations": True, "linker_name": "umls"})

linker_large = nlp_large.get_pipe("scispacy_linker")

################################################################################

## Cache texts as we process them with spaCy, to speed up code
## Create a separate cache for each spaCy model
CACHED_TEXT_BIO = {}  # bionlp13cg_md
CACHED_TEXT_LARGE = {}  # en_core_sci_lg

def get_cached_text(txt: str, spacy_large: bool = False) -> spacy.tokens.doc.Doc:
    """Retrieve cached results for some string, txt, already processed by
    either the bionlp13cg_md or en_core_sci_lg spaCy model.

    Inputs:
        txt: string to retrieve cached spaCy results for

        spacy_large: bool indicating if we want results from en_core_sci_lg
        model. If False, then get cached results from bionlp13cg_md model
    """
    if spacy_large:
        if txt not in CACHED_TEXT_LARGE:
            CACHED_TEXT_LARGE[txt] = nlp_large(txt)
        
        return CACHED_TEXT_LARGE[txt]
    
    if not txt in CACHED_TEXT_BIO:
        CACHED_TEXT_BIO[txt] = nlp_bio(txt)

    return CACHED_TEXT_BIO[txt]

################################################################################

## Cleaning up recorded sites of infection

# Manually parsed from inspecting 'Site of Infection' column
# "Unimportant" parenthesized information following recorded infection sites, to
# remove from recorded infection sites
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

# Parenthesized information following recorded infection sites that IS informative,
# as it tells us more information about where the Microsporidia infects, so keep
# this info
EXCLUSIONS = ['in lamina propria through muscularis mucosae into submucosa of small intestine tract',
              'mainly in duodenum']

# Named entities from en_ner_bionlp13cg_md (nlp_large) that likely correspond to
# Microsporidia infection sites in texts
INFECTION_SITE_ENTITIES = ['ORGAN', 'TISSUE', 'ORGANISM_SUBDIVISION', \
    'MULTI_TISSUE_STRUCTURE', 'CELL', 'CELLULAR_COMPONENT']


def get_abbrev_species_name(host: str) -> str:
    """For a taxonomic species name, host, return its abbreviated species name.
    Ex: Nosema bombycis -> N. bombycis.
    Ex: Aedes punctor punctor -> A. p. punctor

    Keep already abbreviated names as is
    Ex: N. bombycis -> N. bombycis (unchanged)
    """
    if re.search("[A-Z]\. ", host):
        # host name already abbreviated, return as is
        return host
    
    host = host.split(' ')
    abbrevs = [s[0] + '.' for s in host[:len(host) - 1]]

    return ' '.join(abbrevs + [host[len(host) - 1]])


def get_host_names_and_abbrevs(hosts_natural, hosts_experimental) -> List[str]:
    """For semi-colon separated strings of recorded natural and experimental
    hosts for a microsporidia species, return a list of the abbreviated names
    for all the species.

    Inputs:
        hosts_natural: semi-colon separated string of natural infection hosts

        hosts_experimental: semi-colon separated string of experimental
        infection hosts
    """
    hosts_natural = [h.split('(')[0].strip() for h in hosts_natural.split('; ')]
    hosts_experimental = [h.split('(')[0].strip() for h in hosts_experimental.split('; ')]
    all_hosts = hosts_natural + hosts_experimental
    all_hosts.extend([get_abbrev_species_name(h) for h in all_hosts])

    return(list(filter(lambda s: s != '', all_hosts)))


def should_remove_subsite(subsite: str, all_host_names: List[str]) -> bool:
    """Determine whether to remove a subsite (i.e: parenthesized text) associated
    with some recorded infection site.

    Return True (should remove the recorded subsite) if it corresponds to
    an irrelevant term that doesn't have to do with infection sites, or if it
    corresponds to some taxonomic species name.

    Inputs:
        subsite: string of the parenthesized text/subsite info associated with
        a recorded Microsporidia infection site

        all_host_names: a list of all host names (full + abbreviated names) for
        a Microsporidia species, to sus out subsite info that corresponds to
        taxonomic names (and should be removed)
    """
    return any([term in subsite for term in IRRELEVANT_PARENTHESIZED_INFO + all_host_names]) and \
        not any([term in subsite for term in EXCLUSIONS])


def clean_recorded_infection_sites(sites: str, hosts_natural: str,
                                   hosts_experimental: str) -> str:
    """Format recorded sites of infection by ensuring parenthesized text
    are indicating subcellular/subtissue locations only (i.e: remove
    parenthesized text that corresponds to taxonomic host names, and
    other irrelevant info that doesn't inform us more about infection
    sites).

    Return a semi-colon separated string of all infection sites and
    their relevant subsites, all separated by semi-colons.

    Inputs:
        sites: semi-colon separated string of recorded infection sites for
        some Microsporidia species

        hosts_natural: semi-colon separated string of naturally infected
        hosts for some Microsporidia species

        hosts_experimental: semi-colon separated string of experimentally
        infected hosts for some Microsporidia species
    """
    # replace('(1;2)', '') fixes a corner case in a single entry where a
    # semi-colon was accidentally used to separate subsite entries
    # (semi-colons should only separate full infection site entries)
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
            # subsites are recorded for site, indicated by parenthesized text
            subsites = \
                [s.strip().strip(')') for s in re.search("(?<=\().+\)?", site).group(0).split(',')]
        else:
            subsites = []
        
        # only keep informative subsites for each infection site
        subsites = [s for s in subsites if not should_remove_subsite(s, all_host_names)]
        
        # add sites + subsites together into one list
        sites_formatted.append(site.split('(')[0].strip())
        sites_formatted.extend(subsites)
    
    # return a single string containing all infection sites + subsites,
    # separated by semi-colon
    return '; '.join(sites_formatted)

################################################################################

## Predicting sites of infection from paper abstracts

# Base words in a sentence indicates it has info about Microsporidia infection
# sites
INFECTION_LEMMAS = ['find', 'parasitize', 'infect', 'describe', 'localize', 'invade',
                    'appear', 'parasite']
infection_pattern = [{'POS': 'VERB', 'LEMMA': {'IN': INFECTION_LEMMAS}}]
matcher = Matcher(nlp_large.vocab)
matcher.add('infection_site', [infection_pattern])

# Words that correspond to Microsporidia anatomy, and not host infection sites,
# and so shouldn't be part of predicted infection sites of a Microsporidia species
MICROSPORIDIA_PARTS = ['polar tube', 'polar filament', 'sporoblast', 'spore'
                       'meront', 'meronts', 'sporoplasm', 'sporophorous vesicles',
                       'sporont', 'sporonts' 'polaroplast', 'anchoring disc',
                       'lamellar body', 'anchoring disk', 'endospore', 'exospore',
                       'posterior vacuole', 'sporoblasts', 'meiospores', 'meiospore',
                       'macrospore', 'macrospores', 'microspore', 'microspores',
                       'basal', 'spores', 'schizogony', 'sporogony']


def is_overlapping_spans(s1: spacy.tokens.span.Span, s2: spacy.tokens.span.Span) -> \
    bool:
    """Return True if two spaCy spans have overlapping boundaries in a document.
    """
    if len(set(range(s1.start, s1.end + 1)) & \
        set(range(s2.start, s2.end + 1))) > 0:
        return True
    
    return False


def get_overlapping_entity(site: spacy.tokens.span.Span, \
    ents: List[spacy.tokens.span.Span]) -> spacy.tokens.span.Span:
    """Check if some span has any overlap (i.e: is a substring of) any
    entity in some list of entity spans, ents.

    Inputs:
        site: span to check for overlaps with entity list
        ents: list of spans corresponding to spaCy entities of interest
    """
    for ent in ents:
        if is_overlapping_spans(ent, site):
            return ent
    
    return site  # return site 'as is' if no overlapping entity found


def get_site_spans(doc: spacy.tokens.doc.Doc) -> List[spacy.tokens.span.Span]:
    """For a semi-colon separated string of recorded infection sites,
    return the corresponding list of spans for each item in this
    semi-colon separated string

    Input:
        doc: semi-colon separated string of recorded infection sites that
        has been processed by a spaCy pipeline
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
    """Docstring goes here
    """
    if not recorded_names:
        return []

    # get intersecting entries between recorded_names and pred_names
    return [s[0] for s in recorded_names if s[0] in [t[0] for t in pred_names]]


def is_unlabelled_by_bionlp13cg(ent: spacy.tokens.span.Span,
                                bio_infection_sents: List[spacy.tokens.span.Span]) -> bool:
    """Return True if a particular entity is not a named entity in a set
    of given sentences from a text processed by en_ner_bionlp13cg_md.

    Inputs:
        ent: spaCy span for some entity predicted by en_core_sci_lg

        bio_infection_sents: sentences from a text processed by
        en_ner_bionlp13cg_md, which likely contain Microsporidia
        infection site data
    """
    # use Matcher to find spans in each infection sentence matching ent
    # check if matching spans overlap with any entities in each sentence
    # if yes, return False (is labelled by bionlp13cg)
    # if no spans overlap with any labelled entities, then return True
    tmp_matcher = Matcher(nlp_large.vocab)
    tmp_matcher.add('ent', [[{'LOWER': ent.text.lower()}]])

    for sent in bio_infection_sents:
        matches = tmp_matcher(sent)
        if matches:
            # check all matches for overlap with entities in this sentence
            # if any overlap at all, return False
            for match in matches:
                if any([is_overlapping_spans(sent[match[1] : match[2]], ent) for ent in sent.ents]):
                    return False

    return True


def resolve_normalized_names(pred_sites: dict, sites_formatted: str) -> Tuple[str, str]:
    """Get UMLS normalized names for predicted sites + recorded sites.

    TODO - MODIFY FOR NEW PREDICTION ALGORITHM
    """
    # Process recorded infection sites with spaCy, and extract tagged entities
    doc = get_cached_text(sites_formatted)
    sites_formatted = get_site_spans(doc)
    ents_of_interest = [ent for ent in doc.ents if ent.label_ in INFECTION_SITE_ENTITIES]

    # Create dictionary of recorded infection site entities + UMLS normalized names
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
                # set first intersecting UMLS normalized name for predicted and
                # recorded infection site as the normalzied name for both
                canonical_name = linker_large.kb.cui_to_entity[intersecting_umls_names[0]][1]

                site_names[sites_formatted[i]]['canonical_name'] = canonical_name
                site_names[sites_formatted[i]]['needs_normalization'] = False

                pred_sites[pred]['canonical_name'] = canonical_name
                pred_sites[pred]['needs_normalization'] = False
    
    # Return Tuple of semi-colon separated strings for normalized recorded and predicted
    # infection sites
    return (
        # normalized recorded sites
        '; '.join(list(set([site_names[site]['canonical_name'] for site in site_names]))),

        # normalized predicted sites
        '; '.join(list(set([pred_sites[pred]['canonical_name'] for pred in pred_sites])))
    )


def predict_and_normalize_infection_sites_2(txt, sites_formatted) -> Tuple[str, str]:
    """Docstring goes here.
    """

    ############################################################################

    txt_docs = {'bio': {'doc': get_cached_text(txt)},
                'large': {'doc': get_cached_text(txt, spacy_large=True)}}

    # Get sentences likely containing Microsporidia infection site data for each
    # document processed by the two spaCy models
    txt_docs['bio']['infection_sentences'], txt_docs['large']['infection_sentences'] = \
        [sent for sent in txt_docs['bio']['doc'].sents if matcher(sent)], \
            [sent for sent in txt_docs['large']['doc'].sents if matcher(sent)]

    # Detect probable infection entities in "infection sentences" from
    # en_ner_bionlp13cg_md model
    txt_docs['bio']['infection_entities'] = []
    txt_docs['large']['infection_entities'] = []

    for sent in txt_docs['bio']['infection_sentences']:
        # Extract entities in infection sentences that have tags corresponding
        # to probable infection sites, and exclude entities that are actually
        # Microsporidia anatomy parts.
        txt_docs['bio']['infection_entities'].extend(
            [ent for ent in sent.ents if ent.label_ in \
                INFECTION_SITE_ENTITIES and not ent.text.lower() in MICROSPORIDIA_PARTS]
        )

    # Detect probable infection entities in "infection sentences" for
    # en_core_sci_lg model
    #
    # Pick entities meeting the following criteria:
    # 1) Has UMLS linking, so is probably a biologically relevant term, and
    # 2) Is also a relevant infection site entity found by bionlp13cg, OR
    #    has not been tagged as an entity at all by bionlp13cg
    #
    # Use these entities as our "final" set of predicted infection sites from
    # the text
    for sent in txt_docs['large']['infection_sentences']:
        txt_docs['large']['infection_entities'].extend(
            [ent for ent in sent.ents if len(ent._.kb_ents) > 0 and \
                (any([is_overlapping_spans(ent, bio_ent) for bio_ent in txt_docs['bio']['infection_entities']]) or \
                    is_unlabelled_by_bionlp13cg(ent, txt_docs['bio']['infection_sentences']))]
        )

    # sloppy temporary code
    pred_sites = {}
    for ent in txt_docs['large']['infection_entities']:
         pred_sites[ent] = {'umls_names': ent._.kb_ents,
                            'needs_normalization': True,
                            'canonical_name': ent.text}

    ############################################################################

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