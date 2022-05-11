# -----------------------------------------------------------------------------
#
# Predict microsporidia/host species from paper titles and abstracts
#
# Jason Jiang - Created: 2022/05/09
#               Last edited: 2022/05/09
#
# Mideo Lab - Microsporidia text mining
#
# Use Taxonerd to pull species names from microsporidia paper titles and
# abstracts.
#
#
# -----------------------------------------------------------------------------

from taxonerd import TaxoNERD
import pandas as pd
from pathlib import Path 

################################################################################

## Load in abstract + traits for microsporidia species
abstracts = pd.read_csv('../../data/abstracts_traits.csv')

## Replace missing values (NaNs) in first_paper_title and abstract columns of
## abstracts with empty strings (to prevent error with Taxonerd)
abstracts[['first_paper_title', 'abstract']] = abstracts[['first_paper_title', 'abstract']].fillna('')

## Create Taxonerd object for extracting species names from text
taxonerd = TaxoNERD(model="en_ner_eco_biobert",  # use more accurate model
                    prefer_gpu=False,
                    with_abbrev=False)

def predict_species(taxonerd, txt) -> str:
    '''Docstring goes here
    '''
    taxonerd_pred = taxonerd.find_in_text(txt)
    if not taxonerd_pred.empty:
        return ', '.join(taxonerd_pred['text'])

    return ''

## Add column in abstracts for taxonerd predicted species from paper titles and
## abstracts
abstracts = abstracts.assign(  # yeah this code is a mess rn
    title_species = lambda df: df['first_paper_title'].map(
        lambda title: predict_species(taxonerd, title)
    ),
    abstract_species = lambda df: df['abstract'].map(
        lambda abstract: predict_species(taxonerd, abstract)
    )
)

## Test file output
filepath = Path('../../results/taxonerd_test.csv')
abstracts.to_csv(filepath)