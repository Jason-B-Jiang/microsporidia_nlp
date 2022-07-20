# -----------------------------------------------------------------------------
#
# Predict microsporidia localities from paper titles + abstracts: V2
#
# Jason Jiang - Created: 2022/05/25
#               Last edited: 2022/07/20
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

from tkinter import N
from flashgeotext.geotext import GeoText
import spacy
import geocoder
import pickle
from os.path import exists
from typing import Dict, List, Tuple, Optional
import pandas as pd
import re
import copy
from pathlib import Path

################################################################################

## Initialize language models
nlp = spacy.load('en_core_web_md')
geotext = GeoText()

################################################################################

## Global variables
SPACY_LOCALITIES = ['FAC', 'GPE', 'LOC']  # spacy entities corresponding to localities

if exists('geonames_cache.pickle'):
    # load in cached geonames search results
    with open('geonames_cache.pickle', 'r') as f:
        GEONAMES_CACHE = pickle.load(f)
else:
    GEONAMES_CACHE = {}

USERNAME = 'jiangjas'  # fill in your own geonames username here

################################################################################

## Helper functions for predicting regions + subregions from texts

def remove_leading_determinant(span: spacy.tokens.span.Span) -> str:
    """Remove leading determinant from a spacy span for a location entity, and
    return the resulting text.
    Ex: the Weddell Sea -> Weddell Sea

    Input:
        span: spaCy span representing a location entity
    
    Return:
        string of the spaCy span, with leading determinant removed
    """
    if span[0].pos_ == 'DET':
        # remove leading determinant and return resulting text
        return span[1:].text
    
    return span.text


def get_spacy_preds(txt: str) -> List[spacy.tokens.span.Span]:
    """Return a list of spans corresponding to geographical location entities
    predicted by spaCy from a text.

    Input:
        txt: text to predict localities from.
    
    Output:
        List of spans corresponding to unique spaCy predicted locality entities
        Empty list if no spaCy predictions
    """
    # get spaCy locality predictions for txt
    spacy_preds = \
        [ent for ent in nlp(txt).ents if ent.label_ in SPACY_LOCALITIES]

    # keep only unique spaCy locality predictions
    loc_names = []
    for pred in spacy_preds:
        # remove leading determinant (ex: the Ocean of Weddel, remove 'the'),
        # in case the location is mentioned both with and without the determinant
        # in the text
        if remove_leading_determinant(pred) in loc_names:
            spacy_preds.remove(pred)
        else:
            loc_names.append(remove_leading_determinant(pred))
    
    return spacy_preds


def get_top_region_hit(location: str, geonames_result):
    """Top region hit from geonames search result is result with exact
    match to location in question, or the first search result if no
    exact name matches.
    """
    exact_matches = [res for res in geonames_result if \
        res.address.lower() == location.lower()]

    if exact_matches:
        return exact_matches[0]
    
    return geonames_result[0]


def get_most_likely_region(location, geonames_result, geo_preds_regions) -> \
    Optional[Tuple[str, str, bool]]:
    """Using the top 50 geonames results for a location and a list of countries
    already identified by flashgeotext, get the most likely region/country of
    origin for this location.

    If no geonames search results, return None
    """
    if not geonames_result:
        return None

    # likely regions are regions that have already been identified by flashgeotext
    # as countries in the text, and have also been identified by geonames as
    # potential countries for this location
    #
    # sort likely regions with i, so we can pick the first geonames
    # result that overlaps with flashgeotext predicted regions
    #
    # country = country of origin for location, address = "canonical" name for
    # location from geonames
    likely_regions = [(res.country, res.address, i) for res, i in \
        enumerate(geonames_result) if res.country in geo_preds_regions].sort(key = lambda x: x[2])

    top_region_hit = get_top_region_hit(location, geonames_result)

    if not likely_regions:
        # countries detected by flashgeotext in text doesn't correspond to top
        # result found by geonames, so just return the country + normalized
        # name of the top result
        if not top_region_hit.country:
            # no country, so treat location as a region
            return top_region_hit.address, top_region_hit.address
        
        return top_region_hit.country, top_region_hit.address
    
    else:
        # return first search result overlapping with regions predicted
        # by flashgeotext
        return likely_regions[0][0], likely_regions[0][1] 


def set_as_region_or_subregion(loc: str, regions: dict) -> None:
    """Assign a location (loc) as a subregion or a region using geonames.
    """
    region_names = list(regions.keys())

    if not loc in GEONAMES_CACHE:
        # this particular location not looked up in geonames yet, get top 50
        # search results and store in cache
        geonames_result = geocoder.geonames(loc, key=USERNAME, maxRows=50)
        GEONAMES_CACHE[loc] = geonames_result
    else:
        # fetch cached results
        geonames_result = GEONAMES_CACHE[loc]

    most_likely_region, normalized_location_name = \
        get_most_likely_region(loc, geonames_result, region_names)
    
    if not most_likely_region:
        # treat location as a region, if no results from geonames
        regions[loc] = {'subregions': {}, 'found_as': [loc]}
    elif most_likely_region == normalized_location_name:
        # if normalized location name is same as most likely region, then
        # location was a region itself according to geonames
        regions[most_likely_region] = {'subregions': {}, 'found_as': [loc]}
    else:
        # assign location as a subregion to the most likely region
        if most_likely_region in regions:
            if normalized_location_name not in regions[most_likely_region]['subregions']:
                regions[most_likely_region]['subregions'][normalized_location_name] = \
                    [loc]
            else:
                regions[most_likely_region]['subregions'][normalized_location_name].append(loc)
        else:
            regions[most_likely_region] = {'subregions': {normalized_location_name: [loc]}, \
                                           'found_as': [most_likely_region]}


def get_regions_and_subregions(spacy_preds: List[spacy.tokens.span.Span]) -> \
    Dict[str, dict]:
    """For a list of spaCy entities corresponding to probable geographic
    locations, assign them to their likely origin regions or set them as
    their own region.
    
    Return a dictionary with keys as region names, and values as dictionaries
    of subregions for each region and a list of aliased names to this region.
    """
    regions = {}
    undetermined = []

    for pred in spacy_preds:
        geotext_pred = geotext.extract(pred.text)

        if geotext_pred['countries']:
            # geotext predicted this location as a country, so add it as
            # a region using the geotext normalized name as the region
            # key name in the dictionary
            regions[list(geotext_pred['countries'].keys())[0]] = \
                {'subregions': {}, 'found_as': [pred.text]}

        else:
            # for locations tagged as cities or not recognized at all by
            # flashgeotext, add these to 'undetermined' so we can link
            # them back to their respective regions or set them as
            # their own regions
            undetermined.append(pred.text)

    # start assigning undetermined regions back to their regions/countries,
    # or mark them as their own independent regions
    for loc in undetermined:
        set_as_region_or_subregion(loc, regions)

    return regions


def predict_localities(txt: str) ->  Dict[str, dict]:
    """For a text, use flashgeotext + spaCy to predict countries and their
    associated cities/internal regions/etc.

    Input:
        txt: text to predict localities from

    Return:
        Dictionary with keys as regions and values as subregions
    """
    # get spaCy entities from text that correspond to geographical locations
    spacy_preds = get_spacy_preds(txt)

    # remove any location predictions that are actually just taxonomic names,
    # as those tend to get tagged as location entities by spacy
    #
    # Idea: have host/microsporidia species prediction script run first, then
    # use results from that to filter out taxonomic entities, instead of waiting
    # an hour to redo the predictions again here.
    #
    # spacy_preds = remove_taxonomic_entities(spacy_preds)

    # return the predicted locations as a dictionary of regions and their
    # subregions (ex: region = Australia, subregions = [Melbourne, Sydney])
    #
    # TODO - write a function that converts this to a string
    return get_regions_and_subregions(spacy_preds)

################################################################################

## Helper functions for normalizing recorded localities, so they're consistent
## with predicted localities

def create_locality_dictionary(recorded_locs: str) -> dict:
    """For a string of recorded localities from a text, turn this string into
    the same format as returned by predict_localities.
    
    Ex: U.S. (Baltimore); Belarus ->
    
        {'U.S.': {'subregions': {'Baltimore': ['Baltimore']}, 
                  'found_as': ['U.S.']},
         'Belarus': {'subregions': {}, 'found_as': ['Belarus']}}
    """
    recorded_locs = [s.strip() for s in recorded_locs.split(';')]
    locs_dict = {}

    for loc in recorded_locs:
        # extract region name away from parenthesized subregion text
        region = loc.split('(')[0].strip()

        # extract parenthesized subregion text
        subregions = re.search('(?<=\().+(?=\))', loc)
        if subregions:
            # this nasty list comprehension separates all comma/'|' separated
            # elements in the parentheses and returns an unnested list
            subregions = [item for l in [s.split(', ') for s in \
                subregions.group(0).split(' | ')] for item in l]
        else:
            subregions = []
        
        locs_dict[region] = {'subregions': {s: [s] for s in subregions},
                             'found_as': [region]}

    return locs_dict


def flashgeotext_normalize_recorded_localities(recorded_locs_dict: dict,
                                               recorded_locs) -> \
                                                   Tuple[dict, dict]:
    """Use flashgeotext to replace localities names in recorded_locs with
    their canonical geonames names, if possible.
    
    Return a tuple of recorded_locs_dict but with all location names replaced
    by normalized names from flashgeotext, and a dictionary of all undetected
    regions and their subregions
    """
    flashgeotext_preds = geotext.extract(recorded_locs)
    unnormalized_locs = {}
    normalized_locs = copy.deepcopy(recorded_locs_dict)

    for region in recorded_locs_dict:
        normalized_region = ''
        for country in flashgeotext_preds['countries']:
            if region in flashgeotext_preds['countries'][country]['found_as']:
                # set region key name in deep copied dictionary for normalized
                # recorded locations to the canonical name from flashgeotext
                normalized_locs[country] = normalized_locs.pop(region)
                normalized_region = country
                break

        if not normalized_region:
            # region not detected by flashgeotext, so add to unnormalized_locs
            # dictionary for normalization by geonames later
            unnormalized_locs[region] = \
                {'subregions': [s for s in recorded_locs_dict[region]['subregions']],
                 'normalized_region': False}
        else:
            # still add region to unnormalized locs dictionary, so we can keep track
            # of any subregions to this region that can't be normalized
            unnormalized_locs[normalized_region] = \
                {'subregions': [s for s in recorded_locs_dict[region]['subregions']],
                 'normalized_region': True}

        for subregion in recorded_locs_dict[region]['subregions']:
            for city in flashgeotext_preds['cities']:
                # replace entry for subregion in normalized dictionary with
                # canonical name from flashgeotext
                if subregion in flashgeotext_preds['cities'][city]['found_as']:
                    if normalized_region:
                        normalized_locs[normalized_region]['subregions'][city] = \
                            normalized_locs[normalized_region]['subregions'].pop(subregion)

                        # remove subregion from unnormalized_locs as it has been normalized
                        unnormalized_locs[normalized_region]['subregions'].remove(subregion)
                    else:
                        normalized_locs[region]['subregions'][city] = \
                            normalized_locs[region]['subregions'].pop(subregion)
                        
                        unnormalized_locs[region]['subregions'].remove(subregion)
                    
                    break

    return normalized_locs, unnormalized_locs


def geonames_normalize_recorded_localities(unnormalized_locs: Dict[str, List[str]],
                                           recorded_locs_dict: dict) -> None:
    """Look up recorded regions/subregions not detected by flashgeotext in geonames,
    and try to replace their entries in recorded_locs with their geonames canonical
    names.
    
    Modifies recorded_locs in-place.
    """
    for region in unnormalized_locs:
        normalized_region = region  # replace with canonical geonames name, if possible

        if not unnormalized_locs[region]['normalized_region']:
            # get canonical region name from geonames since name wasn't normalized
            # by flashgeotext
            normalized_region = get_geonames_canonical_region(region)
            if normalized_region is None:
                normalized_region = region
            else:
                recorded_locs_dict[normalized_region] = \
                    recorded_locs_dict.pop(region)

                unnormalized_locs[normalized_region] = unnormalized_locs.pop(region)
                unnormalized_locs[normalized_region]['normalized_region'] = True
        
        for subregion in unnormalized_locs[normalized_region]['subregions']:
            normalized_subregion = \
                get_geonames_canonical_subregion(subregion, normalized_region)
            
            if normalized_subregion is not None:
                recorded_locs_dict[normalized_region]['subregions'][normalized_subregion] = \
                    recorded_locs_dict[normalized_region]['subregions'].pop(subregion)
                
                unnormalized_locs[normalized_region]['subregions'].remove(subregion)


def normalize_recorded_localities(recorded_locs: str) -> str:
    """For a semi-colon separated list of regions and subregions recorded from a
    text, make the recorded location names consistent with predicted location names
    using flashgeotext and geonames.
    
    Ex: U.S. (Baltimore | Miami, Florida); Belarus ->
        United States (Baltimore | Miami | Florida); Belarus
    """
    # Create a dictionary of regions to subregions, formatted in the same way
    # predict_localities returns
    recorded_locs_dict = create_locality_dictionary(recorded_locs)

    # Normalize region/subregion names with flashgeotext, if possible, returning
    # a dictionary of undetected regions to undetected subregions, to try and
    # normalize with geonames
    recorded_locs_dict, unnormalized_locs = \
        flashgeotext_normalize_recorded_localities(
            recorded_locs_dict, recorded_locs
            )

    # Normalize regions/subregions not detected by flashgeotext with geonames,
    # by searching them up in geonames and using the canonical name of the
    # best geonames hit
    geonames_normalize_recorded_localities(unnormalized_locs, recorded_locs_dict)

    # TODO - write function turning this into a string (same function can be
    # used for predicting and normalizing)
    return recorded_locs_dict

################################################################################

## Make locality predictions for each microsporidia species paper
## TODO - move all this code to a main() function

# Load in formatted dataframe of microsporidia species data, from
# src/1_format_data/3_misc_cleanup.R
microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Exclude species with >1 papers describing them (for now)
microsp_data = microsp_data[microsp_data['num_papers'] < 2]

# Make locality predictions using title_abstract column
microsp_data = microsp_data.assign(
    # pred = predicted
    pred_locality = lambda df: df['title_abstract'].map(
        lambda txt: predict_localities(txt)[1]
    )
)

microsp_data = microsp_data.assign(
    locality_normalized = lambda df: df['locality'].map(
        lambda locs: normalize_recorded_localities(locs)[1], na_action='ignore'
    )
)

################################################################################

## Write outputs
## TODO - move all this code to a main() function

# Write predictions to results folder
microsp_data[['species', 'title_abstract', 'locality', 'locality_normalized', 'pred_locality']].to_csv(
    Path('../../results/microsp_locality_predictions.csv')
    )

# Save cached geonames search results
# pickling isn't working for me, so I'm commenting out the code
# with open('geonames_cache.pickle', 'w') as f:
#     if GEONAMES_CACHE != {}:  # don't bother writing an empty cache
#         pickle.dump(GEONAMES_CACHE, f)