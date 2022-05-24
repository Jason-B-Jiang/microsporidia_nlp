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

WORDS_TO_DIGITS <- list('one' = 1, 'two' = 2, 'three' = 3, 'four' = 4,
                        'five' = 5, 'six' = 6, 'seven' = 7, 'eight' = 8,
                        'nine' = 9, 'ten' = 10)

PT_COIL_RANGE <- '\\d{1,2}\\.?\\d?( to | or | and | ?- ?| ?– ?| ?— ?|\\/)?\\d*\\.?\\d?'

date_to_range <- function(pt_coils_range) {
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

turn_number_words_to_digits <- function(pt_coil_pred) {
  # ---------------------------------------------------------------------------
  # Replace number words in extracted polar tube coil data with digits
  # ex: 'five to six coils' -> '5 to 6 coils'
  # ---------------------------------------------------------------------------
  # get individual coil predictions
  pt_coil_pred <- str_split(pt_coil_pred, ' \\|\\|\\| ')[[1]]
  
  for (i in 1 : length(pt_coil_pred)) {
    pred <- str_split(pt_coil_pred[i], ' ')[[1]]  # split into words
    
    for (j in 1 : length(pred)) {
      if (tolower(pred[j]) %in% names(WORDS_TO_DIGITS)) {
        pred[j] <- WORDS_TO_DIGITS[[tolower(pred[j])]]
      }
    }
    
    pt_coil_pred <- replace(pt_coil_pred, i, str_c(pred, collapse = ' '))
  }
  
  return(str_c(pt_coil_pred, collapse = ' ||| '))
}


extract_pt_coil_range <- function(pt_coil_pred) {
  return(
    str_c(sapply(str_split(pt_coil_pred, ' \\|\\|\\| ')[[1]],
                         function(s) {str_extract(s, PT_COIL_RANGE)}),
                        collapse = ' ||| ')
  )
}


convert_pt_coil_ranges_to_medians <- function(pred_pt_coil_formatted) {
  pt_coil_preds <- str_split(pred_pt_coil_formatted, ' \\|\\|\\| ')[[1]]
  
  for (i in 1 : length(pt_coil_preds)) {
    # get average/median of a range of polar tube coils
    pt_coil_preds[i] <-
      mean(
        as.numeric(str_split(pt_coil_preds[i], '( to | or | and | ?- ?| ?– ?| ?— ?|\\/)')[[1]])
        )
  }
  
  return(str_c(pt_coil_preds, collapse = '; '))
}


get_pt_coil_precision <- function(pt_coils_avg, pred_pt_coil_avg) {
  pt_coils_avg <- str_split(pt_coils_avg, '; ')[[1]]
  pred_pt_coil_avg <- str_split(pred_pt_coil_avg, '; ')[[1]]
  
  true_pos <- length(pred_pt_coil_avg[pred_pt_coil_avg %in% pt_coils_avg])
  false_pos <- length(pred_pt_coil_avg[!(pred_pt_coil_avg %in% pt_coils_avg)])
  
  return(true_pos / (true_pos + false_pos))
}


get_pt_coil_recall <- function(pt_coils_avg, pred_pt_coil_avg) {
  pt_coils_avg <- str_split(pt_coils_avg, '; ')[[1]]
  pred_pt_coil_avg <- str_split(pred_pt_coil_avg, '; ')[[1]]
  
  true_pos <- length(pred_pt_coil_avg[pred_pt_coil_avg %in% pt_coils_avg])
  false_neg <- length(pt_coils_avg[!(pt_coils_avg %in% pred_pt_coil_avg)])
  
  return(true_pos / (true_pos + false_neg))
}


pt_coil_preds <- read_csv('../../results/microsp_pt_predictions.csv') %>%
  rowwise() %>%
  mutate(pt_coils_range = ifelse(!is.na(pt_coils_range),
                                        date_to_range(pt_coils_range),
                                        NA),
         pred_pt_coil_formatted = extract_pt_coil_range(turn_number_words_to_digits(pred_pt_coil)),
         pred_pt_coil_avg = convert_pt_coil_ranges_to_medians(pred_pt_coil_formatted)) %>%
  filter(!is.na(pred_pt_coil) | !is.na(pt_coils_range) | !is.na(pt_coils_avg)) %>%
  mutate(precision = get_pt_coil_precision(pt_coils_avg, pred_pt_coil_avg),
         recall = get_pt_coil_recall(pt_coils_avg, pred_pt_coil_avg))
  

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