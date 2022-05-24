# -----------------------------------------------------------------------------
#
# Evaluate accuracy of predicted microsporidia traits
#
# Jason Jiang - Created: 2022/05/17
#               Last edited: 2022/05/23
#
# Mideo Lab - Microsporidia text mining
#
#
# -----------------------------------------------------------------------------

library(tidyverse)

################################################################################

## Evaluate polar tube coil predictions

format_coil_range <- function(pt_coils_range) {
  # ---------------------------------------------------------------------------
  # Turn ranges that Excel turned into dates, back into ranges
  # ---------------------------------------------------------------------------
  if (!str_detect(pt_coils_range, '2022-')) {
    return(pt_coils_range)
  }
  
  pt_coils_range <- str_c(
    as.integer(str_split(str_remove(pt_coils_range, '2022-'), '-')[[1]]),
    collapse = '-'
  )
  
  return(pt_coils_range)
}


pt_coil_preds <- read_csv('../../results/microsp_pt_predictions.csv') %>%
  rowwise() %>%
  mutate(pt_coils_range = ifelse(!is.na(pt_coils_range),
                                        format_coil_range(pt_coils_range),
                                        NA)) %>%
  filter(!is.na(pred_pt_coil) | !is.na(pt_coils_range) | !is.na(pt_coils_avg))

################################################################################

## Evaluate locality predictions
## ~36% precision and recall, if going by exact matches for countries only

format_locality <- Vectorize(function(locality) {
  str_c(str_remove(str_split(locality, '; ')[[1]], ' ?\\(.+\\) ?'),
        collapse = '; ')
})


get_true_pos_locality <- Vectorize(function(locality, pred_locality) {
  locality <- str_split(locality, '; ')[[1]]
  pred_locality <- str_split(pred_locality, '; ')[[1]]
  
  true_pos <- pred_locality[which(pred_locality %in% locality)]
  if (length(true_pos) == 0) {
    return(NA)
  } else {
    return(str_c(true_pos, collapse = '; ')[[1]])
  }
}, vectorize.args = c('locality', 'pred_locality'))


get_false_pos_locality <- Vectorize(function(locality, pred_locality) {
  locality <- str_split(locality, '; ')[[1]]
  pred_locality <- str_split(pred_locality, '; ')[[1]]
  
  false_pos <- pred_locality[which(!(pred_locality %in% locality))]
  if (length(false_pos) == 0) {
    return(NA)
  } else {
    return(str_c(false_pos, collapse = '; ')[[1]])
  }
}, vectorize.args = c('locality', 'pred_locality'))


get_false_neg_locality <- Vectorize(function(locality, pred_locality) {
  locality <- str_split(locality, '; ')[[1]]
  pred_locality <- str_split(pred_locality, '; ')[[1]]
  
  false_neg <- locality[which(!locality %in% pred_locality)]
  if (length(false_neg) == 0) {
    return(NA)
  } else {
    return(str_c(false_neg, collapse = '; ')[[1]])
  }
}, vectorize.args = c('locality', 'pred_locality'))


get_precision <- Vectorize(function(true_pos, false_pos) {
  pos <- sapply(c(true_pos, false_pos),
                function(x) {ifelse(is.na(x), 0, length(str_split(x, '; ')[[1]]))})
  
  if (identical(pos, c(0, 0))) {
    return(0)
  } else {
    # precision = TP / (TP + FP)
    return(as.numeric(pos[1] / (pos[1] + pos[2]))) 
  }
}, vectorize.args = c('true_pos', 'false_pos'))


get_recall <- Vectorize(function(true_pos, false_neg) {
  pos <- sapply(c(true_pos, false_neg),
                function(x) {ifelse(is.na(x), 0, length(str_split(x, '; ')[[1]]))})
  
  if (identical(pos, c(0, 0))) {
    return(0)
  } else {
    # precision = TP / (TP + FN)
    return(as.numeric(pos[1] / (pos[1] + pos[2]))) 
  }
}, vectorize.args = c('true_pos', 'false_neg'))


locality_preds <- read_csv('../../results/microsp_locality_predictions.csv') %>%
  filter(!is.na(locality)) %>%
  mutate(locality_formatted = format_locality(locality),
         true_pos = get_true_pos_locality(locality_formatted, pred_locality),
         false_pos = get_false_pos_locality(locality_formatted, pred_locality),
         false_neg = get_false_neg_locality(locality_formatted, pred_locality),
         precision = get_precision(true_pos, false_pos),
         # NaNs being returned instead of 0 for some reason, temporary fix
         precision = ifelse(is.nan(precision), 0, precision),
         recall = get_recall(true_pos, false_neg))

################################################################################

## Evaluate microsporidia species names + hosts predictions
## ~58% of microsporidia names successfully extracted

format_species_name <- Vectorize(function(species) {
  # TODO - allow abbreviated name matches
  species <- str_remove(species, '( \\d+| from .+)')
  
  if (str_detect(species, '\\(')) {
    alt_name <- str_remove_all(str_extract(species, '\\(=?.+\\)'), '(\\(|\\)|= )')
    return(str_c(str_remove(species, ' \\(.+'), '; ', alt_name))
  }
  
  return(species)
})


check_microsp_prediction <- function(species, pred_species) {
  if (is.na(pred_species)) {
    return(FALSE)
  }
  
  # make matches case insensitive
  species <- tolower(str_split(species, '; ')[[1]])
  pred_species <- tolower(str_split(pred_species, '; ')[[1]])
  
  for (sp in species) {
    if (any(str_detect(pred_species, sp))) {
      return(TRUE)
    }
  }
  
  return(FALSE)
}


microsp_host_preds <- read_csv('../../results/microsp_and_host_predictions.csv') %>%
  mutate(species_formatted = format_species_name(species)) %>%
  rowwise() %>%
  mutate(microsp_predicted = check_microsp_prediction(species_formatted, pred_species))