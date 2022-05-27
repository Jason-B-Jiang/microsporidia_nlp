from flashgeotext.geotext import GeoText
import spacy
import geocoder
import pickle
from os.path import exists
from typing import Dict, List
import pandas as pd
from pathlib import Path

################################################################################

## Initialize language models
geotext = GeoText()
nlp = spacy.load('en_core_web_md')

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

def get_flashgeotext_preds(txt: str) -> Dict[str, dict]:
    """Return flashgeotext predictions for countries and cities from a text.

    Rename the 'countries' and 'cities' keys to 'regions' and 'subregions', to
    allow for broader geographical locations
    (ex: region = ocean, subregion = landmark in ocean)

    Input:
        txt: text to predict localities from.
    
    Return:
        Dictionary of countries and cities, as returned by GeoText.extract.
        Empty dictionary if no predictions from flashgeotext.
    """
    geo_preds = geotext.extract(input_text=txt)

    # rename countries to regions, and cities to subregions
    # region = country or any independent location that can't be assigned to
    # a country (ex: Pacific Ocean)
    # subregion = anything within a region, like cities or other landmarks
    geo_preds['regions'] = geo_preds.pop('countries')
    geo_preds['subregions'] = geo_preds.pop('cities')

    # assign subregions to each region later
    # also remove 'span_info' and 'count' key for each region, as we won't use it
    for region in geo_preds['regions']:
        geo_preds['regions'][region]['subregions'] = []
        geo_preds['regions'][region].pop('span_info')
        geo_preds['regions'][region].pop('count')

    # also remove 'span_info' and 'count' keys from subregions
    for subregion in geo_preds['subregions']:
        geo_preds['subregions'][subregion].pop('span_info')
        geo_preds['subregions'][subregion].pop('count')

    return geo_preds


def get_spacy_preds(txt: str, geo_preds: Dict[str, dict]) -> List[spacy.tokens.span.Span]:
    """Return a list of spans corresponding to spaCy locality entities that weren't
    predicted by flashgeotext.

    Input:
        txt: text to predict localities from.
        geo_preds: flashgeotext localitity predictions for txt, from
                   get_flashgeotext_preds.
    
    Output:
        List of spans corresponding to spaCy predicted locality entities that
        weren't picked up by flashgeotext.
        Empty list if no spaCy predictions or all predictions accounted for by
        flashgeotext.
    """
    # get spaCy locality predictions for txt
    spacy_preds = \
        [ent for ent in nlp(txt).ents if ent.label_ in SPACY_LOCALITIES]

    # Set spaCy predictions that have been found by geotext to None
    for i in range(len(spacy_preds)):
        found = False

        # check for spacy prediction in geotext regions, then subregions
        for region in geo_preds['regions']:
            if spacy_preds[i].text in geo_preds['regions'][region]['found_as'] + [region]:
                spacy_preds[i] = None
                found = True
                break

        if found:  # move on to next spacy pred if found in regions
            continue

        for subregion in geo_preds['subregions']:
            if spacy_preds[i].text in geo_preds['subregions'][subregion]['found_as'] + [subregion]:
                spacy_preds[i] = None
                break

    # remove Nones from spacy_preds and return it
    return(list(filter(None, spacy_preds)))


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


def assign_to_region_or_subregion(locations: List[str], geo_preds) -> None:
    """For a list of locations, check if it is a region or subregion of some
    other region, and update geo_preds dictionary accordingly.

    Input:
        locations: List of strings for localities of undetermined region/subregion
        geo_preds: dictionary of flashgeotext locality predictions from a text,
                   generated by get_flashgeo_preds
    
    Return:
        None, mutates geo_preds in-place
    """
    for loc in locations:
        if not loc in GEONAMES_CACHE:
            geonames_result = geocoder.geonames(loc, key = USERNAME)

            # store geonames search result in cache for future access
            GEONAMES_CACHE[loc] = geonames_result
        else:
            # fetch cached geonames search result for location
            geonames_result = GEONAMES_CACHE[loc]

        loc_country = geonames_result.country
        if not loc_country:  # no origin country for location, treat as a region
            geo_preds['regions'][loc] = {'found_as': [loc], 'subregions': []}
        else:  # assign location to origin country as its region
            if loc_country not in geo_preds['regions']:
                # origin country of location from geonames not yet in geo_preds
                geo_preds['regions'][loc_country] = {'found_as': [loc_country],
                                                     'subregions': [loc]}
            else:
                geo_preds['regions'][loc_country]['subregions'].append(loc)


def format_localitity_string(locality_preds: Dict[str, List[str]]) -> str:
    """Format dictionary of region predictions + their subregions as a string,
    in the form of 'Region 1 (subregion A | subregion B); Region 2 (...); ...'
    """
    loc_strs = []
    for region in locality_preds:
        loc_strs.append(
            region + ' (' + ' | '.join(locality_preds[region]['subregions']) + ')'
        )

    return '; '.join(loc_strs)
        

def predict_localities(txt: str) ->  Dict[str, dict]:
    """For a text, use flashgeotext + spaCy to predict countries and their
    associated cities/internal regions/etc.

    Input:
        txt: text to predict localities from

    Return:
        Dictionary of region names, with same keys as geotext['countries']
        + additional subregions key
    """
    # get flashgeotext predictions for countries and cities
    geo_preds = get_flashgeotext_preds(txt)

    # get spaCy locality predictions that weren't predicted by flashgeotext
    spacy_preds = get_spacy_preds(txt, geo_preds)

    # update geo_preds with spaCy predictions, adding regions and subregions
    # to regions as appropriate
    assign_to_region_or_subregion([remove_leading_determinant(pred) for pred in spacy_preds],
                                  geo_preds)

    # update geo_preds with subregions from itself, adding the subregions to
    # its corresponding regions as appropriate
    assign_to_region_or_subregion(list(geo_preds['subregions'].keys()), geo_preds)

    return format_localitity_string(geo_preds['regions'])

################################################################################

## Make locality predictions for each microsporidia species paper

# Load in formatted dataframe of microsporidia species data, from
# src/1_format_data/3_misc_cleanup.R
microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Exclude species with >1 papers describing them (for now)
microsp_data = microsp_data[microsp_data['num_papers'] < 2]

# Make locality predictions using title_abstract column
microsp_data = microsp_data.assign(
    # pred = predicted
    pred_locality = lambda df: df['title_abstract'].map(
        lambda txt: predict_localities(txt)
    )
)

################################################################################

## Write outputs

# Write predictions to results folder
microsp_data[['species', 'title_abstract', 'locality', 'pred_locality']].to_csv(
    Path('../../results/microsp_locality_predictions.csv')
    )

# Save cached geonames search results
# pickling isn't working for me, so I'm commenting out the code
# with open('geonames_cache.pickle', 'w') as f:
#     if GEONAMES_CACHE != {}:  # don't bother writing an empty cache
#         pickle.dump(GEONAMES_CACHE, f)