# -----------------------------------------------------------------------------
#
# Predict microsporidia localities from paper titles + abstracts: V2
#
# Jason Jiang - Created: 2022/05/25
#               Last edited: 2022/07/06
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

from flashgeotext.geotext import GeoText
import spacy
import geocoder
import pickle
from os.path import exists
from typing import Dict, List, Tuple, Optional
import pandas as pd
import re
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
            regions[most_likely_region] = {'subregions': {normalized_location_name: [loc]}, 'found_as': [most_likely_region]}


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

    # start assigning cities back to their regions/countries, and either
    # assign undetermined locations to their respective regions, or denote
    # them as their own regions
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
    # spacy_preds = remove_taxonomic_entities(spacy_preds)

    # return the predicted locations as a dictionary of regions and their
    # subregions (ex: region = Australia, subregions = [Melbourne, Sydney])
    return get_regions_and_subregions(spacy_preds)

################################################################################

def get_locality_dict(locs: str) -> Dict[str, List[str]]:
    """Converted recorded localities into dictionary, with regions as keys and
    lists of associated subregions as values.
    """
    locs = [s.strip() for s in locs.split(';')]
    locs_dict = {}

    for loc in locs:
        region = re.search('((?<=\) ).+|.+?(?= \())', loc)
        if region:
            region = region.group(0)
        else:
            region = loc

        # check if recorded region is actually a subregion of some place
        # ex: Siberia is a subregion of Russia
        geonames_res = geocoder.geonames(region, key=USERNAME)
        if geonames_res and geonames_res[0].country and geonames_res[0].country != region:
            locs_dict[geonames_res[0].country] = [region]
            region = geonames_res[0].country
        else:
            locs_dict[region] = []

        subregions = re.search('(?<=\().+(?=\))', loc)
        if subregions:
            subregions = subregions.group(0).split(' | ')
        else:
            subregions = []
        
        locs_dict[region].extend(subregions)

    return locs_dict


def get_locs_str(locs) -> str:
    """Docstring goes here.
    """
    loc_str = []
    for loc in locs:
        subs = []
        for sub in locs[loc]:
            if isinstance(sub, list):
                subs.append(', '.join(sub))
            else:
                subs.append(sub)
        
        loc_str.append(loc + ' (' + ' | '.join(subs) + ')')
    
    return '; '.join(loc_str)


def normalize_recorded_localities(locs: str) -> dict:
    """For a manually recorded entry of localities for a microsporidia species,
    normalize region and subregion names to what is found by flashgeotext,
    to allow for direct comparisions between recorded and predicted localities.

    Input:
        locs: string of recorded localities ('Locality' column) from supp. table
              S1 from Murareanu et al.
    Return:
        string formatted in same way as locs, but regions + subregions are
        renamed with normalized names from flashgeotext
    """
    loc_normalized_names = predict_localities(locs)[0]
    locs = get_locality_dict(locs)

    # list of lists for each region and subregion, where each sublist is
    # all 'found_as' names for the region and subregion
    all_region_names =\
        [[region] + loc_normalized_names[region]['found_as'] for \
            region in loc_normalized_names]

    all_subregion_names = []
    for region in loc_normalized_names:
        all_subregion_names.extend(loc_normalized_names[region]['subregions'])

    # replace region + subregion names in locs w/ normalized names from
    # loc_normalized_names
    regions = list(locs.keys())
    for region in regions:
        region_matches = [i for i, names in enumerate(all_region_names) if \
            any([name for name in names if name in region])]

        if region_matches:
            region_normalized = all_region_names[region_matches[0]][0]
        else:
            # region wasn't predicted by spacy + flashgeonames
            region_normalized = region
        
        for i in range(len(locs[region])):
            # get indices in all_subregion_names where subregion is found
            subregion_matches = [j for j, names in enumerate(all_subregion_names) \
                if [name for name in names if name in locs[region][i]]]
            
            # replace subregions w/ normalized names if found
            # do nothing if no matches for subregion found
            if subregion_matches:
                # more than 1 subregion associated with this recorded subregion
                # ex: a city and its associated province was recorded
                # replace old entry with list
                if len(subregion_matches) > 1:
                    # take first entry in each sublist as "normalized" subregion
                    # name
                    locs[region][i] =\
                        [all_subregion_names[i][0] if \
                            isinstance(all_subregion_names[i], list) \
                                else all_subregion_names[i] for i in subregion_matches]
                else:
                    locs[region][i] = all_subregion_names[subregion_matches[0]][0]
        
        locs[region_normalized] = locs.pop(region)

    return locs, get_locs_str(locs)

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

# Write predictions to results folder
microsp_data[['species', 'title_abstract', 'locality', 'locality_normalized', 'pred_locality']].to_csv(
    Path('../../results/microsp_locality_predictions.csv')
    )

# Save cached geonames search results
# pickling isn't working for me, so I'm commenting out the code
# with open('geonames_cache.pickle', 'w') as f:
#     if GEONAMES_CACHE != {}:  # don't bother writing an empty cache
#         pickle.dump(GEONAMES_CACHE, f)