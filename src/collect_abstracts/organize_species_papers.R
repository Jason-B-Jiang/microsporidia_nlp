# -----------------------------------------------------------------------------
#
# Organize microsporidia species papers
#
# Jason Jiang - Created: 2022/05/02
#               Last edited: 2022/05/09
#
# Mideo Lab - Microsporidia text mining
#
# Select microsporidia species for analysis and extract first papers describing
# each species.
#
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(writexl)

################################################################################

# Supplemental table S1 from Murareanu et al. 2021 (Microsporidia species dataset)
microsp_data <- read_csv('../../data/microsporidia_species.csv')

# Exclude any species coming from these papers
# I'll talk more about this during lab meeting
excluded_papers <- readLines('../../data/excluded_papers.txt')

microsp_data <- microsp_data %>%
  rename(species = `Species Name`) %>%
  mutate(year_first_described = get_year_first_identified(`Date Identified (year)`),
         first_paper_reference = get_first_reference(References),
         first_paper_title = NA,
         abstract = NA,
         notes = NA,
         # Is the paper in a foreign language?
         foreign = NA) %>%
  # filter out cases where year first described is unknown or ambiguous
  # (NA or '?' in the microsporidia dataset for Date Identified)
  filter(!is.na(year_first_described),
         # Filter microsporidia species to species from 1977 to 2021
         # Species from before 1977 come from Sprague book
         year_first_described >= 1977) %>%
  rowwise() %>%
  filter(!any(str_detect(References, excluded_papers))) %>%
  select(species, year_first_described, first_paper_reference,
         first_paper_title, abstract, notes, foreign) %>%
  arrange(year_first_described)

# Save microsp_dataframe as xlsx and manually collect paper abstracts
# (I'll explain why during lab meeting)
write_xlsx(microsp_data, '../../data/manually_collect_abstracts_2.xlsx')

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