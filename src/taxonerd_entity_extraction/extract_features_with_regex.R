# -----------------------------------------------------------------------------
#
# Text
#
# Jason Jiang - Created: 2022/05/10
#               Last edited: 2022/05/10
#
# Mideo Lab - Microsporidia text mining
#
# Text
#
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(glue)  # For formatting strings like in Python

################################################################################

## Global variables for regex patterns to extract spore length/width data
NUMBER <- '\\d\\.?\\d*'  # ex: 4, 5.0, 4.05, etc

DIMENSION  <-  # ex: 4.0 (3.0 - 5.0) um or 4 - 5 (4.5)
  glue('{number}\\s*\\W?\\s*{number}?\\(?\\s*{number}\\s*\\W\\s*{number}\\s*\\)\\s*\\W?[[:alpha:]]*')

LENGTH_WIDTH <-  # ex: 4.0 (3.0 - 5.0) um x 4 - 5 (4.5) um
  glue('{dimension}\\s*\\W*\\s*{dimension}\\s*\\W?\\w*')

################################################################################

## Extract dimensions

extract_spore_dimensions <- Vectorize(function(abstract) {
  # ---------------------------------------------------------------------------
  # Blah
  # ---------------------------------------------------------------------------
  return(str_c(str_extract_all(abstract, LENGTH_WIDTH)[[1]],
               collapse = ' | '))
})

microsp_data <- read_csv('../../data/abstracts_traits.csv') %>%
  rename(avg_length = `Spore Length Average (µm)`,
         avg_width = 'Spore Width Average (µm)') %>%
  select(species, year_first_described, first_paper_title, abstract,
         avg_length, avg_width) %>%
  mutate(extracted_dimensions = extract_spore_dimensions(abstract))

################################################################################

## Extract polar tube coils

## Try extracting numeric data immediately before or after detecting polar tube
## coil term?

detect_coils <- Vectorize(function(abstract) {
  # coil or turn appears as individual words in the abstract, suggesting polar
  # tube coil/turn data
  return(str_detect(abstract, ' *coil[s]? *') | str_detect(abstract, ' *turn[s]? *'))
})

microsp_data <- read_csv('../../data/abstracts_traits.csv') %>%
  select(species, year_first_described, first_paper_title, abstract,
         `Polar Tubule Coils Range`, `Polar Tubule Coils Average`) %>%
  mutate(has_coil_data = detect_coils(abstract))

false_neg <-
  filter(microsp_data, !is.na(`Polar Tubule Coils Range`) | !is.na(`Polar Tubule Coils Average`), !has_coil_data)

false_pos <-
  filter(microsp_data, is.na(`Polar Tubule Coils Range`) & is.na(`Polar Tubule Coils Average`), has_coil_data)
