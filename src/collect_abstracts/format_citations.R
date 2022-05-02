# ------------------------------------------------------------------------------
# Format citations in Microsporidia species dataset for collecting abstracts
# ------------------------------------------------------------------------------

library(tidyverse)
library(rentrez)

################################################################################

microsp_data <- read_csv('../../data/microsporidia_species.csv')

# Filter microsporidia species to species from 1977 to 2021
# Species from before 1977 come from Sprague book
# We can manually filter out 'Sprague' species from our 1977 collection
microsp_data <- microsp_data %>%
  mutate(year_first_described = get_year_first_identified(`Date Identified (year)`)) %>%
  # filter out cases where year first described is unknown or ambiguous
  # (NA or '?' in the microsporidia dataset for Date Identified)
  filter(!is.na(year_first_described), year_first_described >= 1977)

################################################################################

### Helper functions

get_year_first_identified <- Vectorize(function(years) {
  # ----------------------------------------------------------------------------
  # Get the first year a Microsporidia species was described in, as an integer
  #
  # Input:
  #   years: entry from the 'Date Identified (year)' column from microsporidia
  #          species dataset
  # ----------------------------------------------------------------------------
  return(as.integer(str_remove(str_split(years, '; ')[[1]][1], ' \\(.+')))
})