# -----------------------------------------------------------------------------
#
# Format manually collected paper title/abstract data for microsporidia
#
# Jason Jiang - Created: 2022/05/18
#               Last edited: 2022/05/18
#
# Mideo Lab - Microsporidia text mining
#
# Clean extracted paper titles + abstracts, and add back in microsporidia trait
# data
#
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)

################################################################################

## Helper functions
clean_text <- Vectorize(function(text) {
  # ---------------------------------------------------------------------------
  # Replace tabs and newlines with spaces, remove excess + trailing whitespacee
  # from a string.
  #
  # Input:
  #   text: paper title or abstract
  # ---------------------------------------------------------------------------
  trimws(str_replace_all(
    str_replace_all(str_replace_all(text, '\t', ' '), '\n', ' '),
    ' +',
    ' '))
})

################################################################################

# Load in dataframe of manually extracted microsporidia paper titles + abstracts,
# from select_species_papers.R
abstracts <- read_xlsx('../../data/manually_collect_abstracts.xlsx') %>%
  # fix weird issue with true/false being turned into dates by xlsx format
  mutate(foreign = ifelse(!is.na(foreign), TRUE, FALSE),
         paper_title = clean_text(paper_title),
         abstract = clean_text(abstract),
         # create new column combining paper title and abstract
         title_abstract = ifelse(!is.na(abstract),
                                 str_c(paper_title, ' ::: ', abstract),
                                 paper_title)) %>%
  filter(!is.na(paper_title)) %>%
  select(-notes)

# Load in supplemental table S1 from Murareanu et al., 2021
microsp_data <- read_csv('../../data/microsporidia_species.csv') %>%
  select(-Timestamp, -`Host Environment`, -`Transmission`,
         -`Calculated Volume (µm³) (see methods)`,
         -`Calculated Polar Tubule (μm) (see methods)`,
         -`18S Accession #`, -`Has the genome been sequenced?`,
         -`Important Remarks`, -`Date Identified (year)`) %>%
  rename(species = `Species Name`, hosts = `Natural Host(s)`,
         hosts_natural = `Natural Host(s)`,
         hosts_experimental = `Experimental Host(s)`,
         host_lifestage = `Host Life Stage during Infection`,
         infection_site = `Site of Infection`,
         spore_length_avg = `Spore Length Average (µm)`,
         spore_width_avg = `Spore Width Average (µm)`,
         spore_shape = `Spore Shape (Class; Condition)`,
         locality = `Locality`,
         nucleus = Nucleus,
         # pt = polar tube
         pt_max = `Measured Polar Tubule Length Max (μm)`,
         pt_min = `Measured Polar Tubule Length Min (μm)`,
         pt_avg = `Measured Polar Tubule (μm)`,
         pt_coils_range = `Polar Tubule Coils Range`,
         pt_coils_avg = `Polar Tubule Coils Average`,
         all_refs = References)

################################################################################

# Merge abstracts and microsp_data dataframes by species names, to associate
# microsporidia traits with their papers
abstracts_traits <- merge(x = abstracts, y = microsp_data,
                          by = 'species', all = TRUE) %>%
  filter(!is.na(paper_title))

write_csv(abstracts_traits, '../../data/abstracts_traits.csv')