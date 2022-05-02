# ------------------------------------------------------------------------------
# Format citations in Microsporidia species dataset for collecting abstracts
# ------------------------------------------------------------------------------

library(tidyverse)
library(rentrez)

################################################################################

microsp_data <- read_csv('../../data/microsporidia_species.csv')

# Exclude any species coming from these papers
# Write reason why later
excluded_papers <- readLines('../../data/excluded_papers.txt')

# Filter microsporidia species to species from 1977 to 2021
# Species from before 1977 come from Sprague book
# We can manually filter out 'Sprague' species from our 1977 collection
microsp_data <- microsp_data %>%
  mutate(year_first_described = get_year_first_identified(`Date Identified (year)`),
         first_paper = get_first_reference(References)) %>%
  # filter out cases where year first described is unknown or ambiguous
  # (NA or '?' in the microsporidia dataset for Date Identified)
  filter(!is.na(year_first_described), year_first_described >= 1977) %>%
  rowwise() %>%
  filter(!any(str_detect(References, excluded_papers)))

write_csv(select(microsp_data, `Species Name`, year_first_described, first_paper) %>%
            arrange(year_first_described),
          '../../data/manually_add_PMIDs.csv')

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


get_first_reference <- Vectorize(function(ref) {
  # ----------------------------------------------------------------------------
  # Get reference for the first paper describing a particular Microsporidia species
  #
  # Input:
  #   ref: entry from References column of Microsporidia dataset
  # ----------------------------------------------------------------------------
  return(str_remove(trimws(str_split(ref, '\n')[[1]][1]), '^\\d\\. '))
})