# -----------------------------------------------------------------------------
#
# Predict microsporidia localities from paper titles + abstracts
#
# Jason Jiang - Created: 2022/05/19
#               Last edited: 2022/05/19
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

import spacy
from spacy.matcher import PhraseMatcher
import pandas as pd
import pycountry  # for making list of country names
from pathlib import Path

################################################################################

##Initialize spaCy English model + PhraseMatcher for detecting country names
nlp = spacy.load('en_core_web_md')
country_matcher = PhraseMatcher(nlp.vocab)

################################################################################

## Construct list of country names for country_matcher
countries = []
for country in pycountry.countries:
    countries.append(country.name)

for country in pycountry.historic_countries:  # for historic countries like USSR
    countries.append(country.name)

# Manually add abbreviations for USA
countries.extend(['U.S.A.', 'U.S.A', 'U.S.', 'U.S', 'US'])

# Add countries to phrase matcher
country_matcher.add("country", [nlp.make_doc(c) for c in countries])

# Write function for extracting country names from texts
def get_countries(txt: str,
                  matcher: spacy.matcher.phrasematcher.PhraseMatcher) -> str:
    doc = nlp(txt)
    matches = matcher(doc)

    return '; '.join(
        # convert to set then back to list so each country appears only once
        # in list
        list(set([doc[start:end].text for id, start, end in matches]))
        )

################################################################################

## Predict microsporidia countries from titles + abstracts

# Load in formatted dataframe of microsporidia species data, from
# src/1_format_data/3_misc_cleanup.R
microsp_data = pd.read_csv('../../data/manually_format_multi_species_papers.csv')

# Exclude species with >1 papers describing them (for now)
microsp_data = microsp_data[microsp_data['num_papers'] < 2]

# Add column for predicted country names in paper titles + abstracts
microsp_data = microsp_data.assign(
    # pred = predicted
    pred_locality = lambda df: df['title_abstract'].map(
        lambda txt: get_countries(txt, country_matcher)
    )
)

################################################################################

# Write predictions to results folder
microsp_data[['species', 'title_abstract', 'locality', 'pred_locality']].to_csv(
    Path('../../results/microsp_locality_predictions.csv')
    )