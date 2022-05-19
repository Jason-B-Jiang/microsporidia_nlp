# -----------------------------------------------------------------------------
#
# Miscellaneous cleanup
#
# Jason Jiang - Created: 2022/05/18
#               Last edited: 2022/05/18
#
# Mideo Lab - Microsporidia text mining
#
# For clean-up stuff I just thought of
#
#
# -----------------------------------------------------------------------------

library(tidyverse)

################################################################################

abstracts_traits <- read_csv('../../data/abstracts_traits.csv')

abstracts_traits <- abstracts_traits %>%
  rowwise() %>%
  mutate(num_papers = length(str_split(all_refs, '\n')[[1]])) %>%
  separate_rows(all_refs, sep = '\n') %>%
  mutate(all_refs = trimws(all_refs))

################################################################################

COLUMNS_OF_INTEREST <- c('species', 'year_first_described', 'paper_ref', 'all_refs',
                         'num_papers')

COLUMNS_TO_CLEAR <-
  colnames(abstracts_traits)[!(colnames(abstracts_traits) %in% COLUMNS_OF_INTEREST)]

for (i in 1 : nrow(abstracts_traits)) {
  # if num_papers > 1:
  #   set columns of non-interest to NA
  #   add numbering to species w/ refs numbering
  #   add paper_ref
  if (abstracts_traits$num_papers[i] > 1) {
    abstracts_traits[i, 'species'] <-
      str_c(abstracts_traits$species[i],
            ' (',
            str_extract(abstracts_traits$all_refs[i], '^\\d+(?=\\.)'),
            ')') 
    
    abstracts_traits[i, 'paper_ref'] <-
      str_remove(abstracts_traits$all_refs[i], '^\\d+\\. +')
    
    # if not entry for first paper, clear all columns for COLUMNS_TO_CLEAR an
    if (!str_detect(abstracts_traits$all_refs[i], '^1\\.')) {
      abstracts_traits[i, COLUMNS_TO_CLEAR] <-
        as.list(rep(NA, length(COLUMNS_TO_CLEAR)))
    }
  }
}

# Save and manually clean up
write_csv(abstracts_traits, '../../data/manually_format_multi_species_papers.csv')