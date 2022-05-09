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

################################################################################

## Load in abstract + host data for microsporidia species
abstracts = pd.read_csv('../../data/abstracts_hosts_matched.csv')

## Create Taxonerd object for extracting species names from text
taxonerd = TaxoNERD(model="en_ner_eco_biobert",  # use more accurate model
                    prefer_gpu=False,
                    with_abbrev=False)