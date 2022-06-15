# -----------------------------------------------------------------------------
#
# Evaluate accuracy of predicted microsporidia traits
#
# Jason Jiang - Created: 2022/05/17
#               Last edited: 2022/06/10
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

## Evaluate polar tube length predictions

extract_pt_length_value <- Vectorize(function(pred_pt) {
  # ---------------------------------------------------------------------------
  # Extract polar tube length from string where length was found, and convert
  # ranges of length to averages.
  # ---------------------------------------------------------------------------
  mean(as.numeric(
    str_split(str_extract(pred_pt, '\\d+\\.?\\d?(–|-| to )?\\d*\\.?\\d?'),
            '(–|-| to )')[[1]]
  ))
})


format_recorded_pt_length <- Vectorize(function(pt_avg) {
  as.numeric(str_split(pt_avg, ' ')[[1]][1])
})

pt_len_preds <- read_csv('../../results/microsp_pt_predictions.csv') %>%
  select(species, title_abstract, pred_pt, pt_max, pt_min, pt_avg) %>%
  filter(!is.na(pred_pt) | !is.na(pt_max) | !is.na(pt_min) | !is.na(pt_avg)) %>%
  mutate(pred_pt_formatted = ifelse(is.na(pred_pt), NA, extract_pt_length_value(pred_pt)),
         pt_avg_formatted = format_recorded_pt_length(pt_avg),
         tp = ifelse(!is.na(pred_pt_formatted),
                     as.numeric(pred_pt_formatted == pt_avg_formatted),
                     0),
         fp = as.numeric(!is.na(pred_pt_formatted) & pred_pt_formatted != pt_avg_formatted),
         fn = as.numeric(is.na(pred_pt_formatted) & !is.na(pt_avg_formatted)))

pt_len_precision <- sum(pt_len_preds$tp) / (sum(pt_len_preds$tp) + sum(pt_len_preds$fp))
pt_len_recall <- sum(pt_len_preds$tp) / (sum(pt_len_preds$tp) + sum(pt_len_preds$fn))

################################################################################

## Evaluate locality predictions

# In general, we want to see if predictions are a substring of the recorded
# data
# This is because we expect recorded regions to be a little more descriptive,
# while predicted regions should be more concise as we are extracting
# spaCy entities and dictionary matches for countries/subregions

# TP: number of regions + subregions that are actually in recorded data
  # for subregions, allow prediction to be substring of recorded region, as
  # recorded regions tend to be very descriptive

# FP: number of regions + subregions not in recorded data
  # i.e: prediction is not a substring of any of the recorded data

# FN: number of regions + subregions in recorded data but not in predictions
  # i.e: recorded region is not a substring of any predictions, and no predictions
  #      are substrings of the recorded region

get_false_pos_locality <- function(locality, pred_locality) {
  locality <- str_split(locality, '; ')[[1]]
  pred_locality <- str_split(pred_locality, '; ')[[1]]
  false_pos <- c()
  
  for (pred in pred_locality) {
    pred_region <- str_remove(pred, ' \\(.*\\)')
    pred_subregions <-
      str_split(str_extract(pred, '(?<=\\().*(?=\\))'), ' \\| ')[[1]]
    
    # 
  }
}

microsp_locality_preds <- read_csv('../../results/microsp_locality_predictions.csv') %>%
  mutate(false_pos = get_false_pos_locality(locality, pred_locality),
         true_pos = get_true_pos_locality(locality, pred_locality),
         false_neg = get_false_neg_locality(locality, pred_locality),
         precision = get_locality_precision(true_pos, false_pos),
         recall = get_locality_recall(true_pos, false_neg))

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

# Only 49% accuracy w/ genus species name extraction
microsp_host_preds <- read_csv('../../results/microsp_and_host_predictions.csv') %>%
  mutate(species_formatted = format_species_name(species)) %>%
  rowwise() %>%
  mutate(microsp_predicted = check_microsp_prediction(species_formatted, pred_species),
         microsp_predicted_2 = check_microsp_prediction(species_formatted, pred_species_2))

################################################################################

## Evaluate Microsporidia infection site predictions

get_infection_site_tp <- function(recorded, pred) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  if (is.na(pred)) {
    return(0)
  }
  
  if (is.na(recorded)) {
    recorded <- character(0)
  } else {
    recorded <- str_split(recorded, '; ')[[1]]
  }
  
  pred <- str_split(pred, '; ')[[1]]
  
  return(length(recorded[recorded %in% pred]))
}


get_infection_site_fp <- function(recorded, pred) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  if (is.na(pred)) {
    return(0)
  }
  
  if (is.na(recorded)) {
    recorded <- character(0)
  } else {
    recorded <- str_split(recorded, '; ')[[1]]
  }

  pred <- str_split(pred, '; ')[[1]]
  
  return(length(pred[!(pred %in% recorded)]))
}


get_infection_site_fn <- function(recorded, pred) {
  # ---------------------------------------------------------------------------
  # Docstring goes here
  # ---------------------------------------------------------------------------
  if (is.na(pred)) {
    return(0)
  }
  
  if (is.na(recorded)) {
    recorded <- character(0)
  } else {
    recorded <- str_split(recorded, '; ')[[1]]
  }
  
  pred <- str_split(pred, '; ')[[1]]
  
  return(length(recorded[!(recorded %in% pred)]))
}

infection_sites_in_abstract <- function(title_abstract, infection_site) {
  infection_site <- sapply(str_split(infection_site, '; ')[[1]],
                           function(x) {tolower(str_remove(x, ' ?\\(.+\\)'))})
  
  return(any(str_detect(tolower(title_abstract), infection_site)))
}


microsp_infection_preds <- read_csv('../../results/microsp_infection_site_predictions.csv') %>%
  filter(num_papers < 2) %>%  # look at species with only 1 paper for now
  rowwise() %>%
  mutate(tp = get_infection_site_tp(infection_site_normalized, pred_infection_site),
         fp = get_infection_site_fp(infection_site_normalized, pred_infection_site),
         fn = get_infection_site_fn(infection_site_normalized, pred_infection_site),
         infection_site_in_abstract = infection_sites_in_abstract(title_abstract, infection_site))

# 28% precision, 36% recall
# 40% precision, 50% recall if exclude cases where tissues not in abstract
infection_precision <-
  sum(microsp_infection_preds$tp) / (sum(microsp_infection_preds$tp) + sum(microsp_infection_preds$fp))

infection_recall <-
  sum(microsp_infection_preds$tp) / (sum(microsp_infection_preds$tp) + sum(microsp_infection_preds$fn))

################################################################################

## Evaluate predicted localities

get_region_subregions <- function(locality) {
  locality <- str_split(locality, '; ')[[1]]
  region_subregions <- list()
  
  for (loc in locality) {
    region <- str_extract(loc, '.+(?= \\()')[[1]]
    subregions <- str_extract(loc, '(?<=\\().+?(?=\\))')
    
    if (!is.na(subregions)) {
      region_subregions[[region]] <- str_split(subregions, ' \\| ')[[1]]
    } else {
      region_subregions[[region]] <- character(0)
    }
  }
  
  return(region_subregions)
}


get_num_regions_and_subregions <- function(locality) {
  if (is.na(locality)) {
    return(0)
  }
  
  locality <- get_region_subregions(locality)
  
  # number of regions + number of all subregions
  return(length(locality) + sum(unlist(lapply(locality, length))))
}


get_locality_tp <- function(pred, recorded) {
  if (is.na(pred) | is.na(recorded)) {
    return(0)
  }
  
  pred <- get_region_subregions(pred)
  recorded <- get_region_subregions(recorded)
  
  tp <- 0
  for (region in names(recorded)) {
    # true positive if recorded region among predicted regions
    tp <- tp + (region %in% names(pred))
    
    for (subregion in recorded[[region]]) {
      # true positive subregion if any of the predicted subregions for this
      # region are a substring of the recorded subregion
      if (!is.null(pred[[region]])) {
        tp <- tp + any(str_detect(subregion, pred[[region]]))
      }
    }
  }
  
  return(tp)
}


get_locality_fp <- function(pred, recorded) {
  if (is.na(pred)) {
    return(0)
  } else if (is.na(recorded)) {
    # if no recorded localities, then all predicted regions and subregions are
    # false pos
    return(get_num_regions_and_subregions(pred))
  }
  
  pred <- get_region_subregions(pred)
  recorded <- get_region_subregions(recorded)
  
  fp <- 0
  for (region in names(pred)) {
    if (!(region %in% names(recorded))) {
      # predicted region was not a recorded region, so count as false pos
      fp <- fp + 1
      
      # because predicted region is incorrect, all predicted subregions are
      # also false pos
      fp <- fp + length(pred[[region]])
      
    } else {
      # check if predicted subregions are correct for this correctly predicted
      # region
      for (subregion in pred[[region]]) {
        fp <- fp + !any(sapply(recorded[[region]],
                               function(x) {str_detect(x, subregion)}))
      }
    }
  }
  
  return(fp)
}


get_locality_fn <- function(pred, recorded) {
  if (is.na(pred)) {
    # if no predicted regions/subregions, then num false pos is all recorded
    # regions and subregions
    return(get_num_regions_and_subregions(recorded))
  } else if (is.na(recorded)) {
    return(0)
  }
  
  pred <- get_region_subregions(pred)
  recorded <- get_region_subregions(recorded)
  
  fn <- 0
  for (region in names(recorded)) {
    if (!(region %in% names(pred))) {
      # region wasn't in predictions, so count as false negative
      # all associated subregions are also false negatives
      fn <- fn + 1
      fn <- fn + length(recorded[[region]])
      
    } else {
      for (subregion in recorded[[region]]) {
        if (!(subregion %in% pred[[region]])) {
          # none of the predicted subregions for this region are a substring
          # of the recorded subregion, so count as false neg
          fn <- fn + !any(str_detect(subregion, pred[[region]]))
        }
      }
    }
  }
  
  return(fn)
}


check_if_regions_in_abstract <- function(title_abstract, locality) {
  if (is.na(locality)) {
    return(FALSE)
  }
  
  locality <- str_split(locality, '; ')[[1]]
  regions <- sapply(locality, function(x) {str_extract(x, '.+(?= \\()')[[1]]})
  
  return(any(str_detect(title_abstract, regions)))
}


microsp_locality_preds <- read_csv('../../results/microsp_locality_predictions.csv') %>%
  # manual fix for correctly associating Siberia with Russia
  mutate(locality_normalized = ifelse(locality_normalized == 'Siberia ()',
                                      'Russia (Siberia)',
                                      locality_normalized),
         locality_normalized = ifelse(locality_normalized == 'Western Siberia (Novosibirsk)',
                                      'Russia (Novosibirsk, Western Siberia)',
                                      locality_normalized),
         locality_normalized = ifelse(locality_normalized == '(pond (Tom River, Tomsk, Western Siberia)',
                                      'Russia (Tom River, Tomsk, Western Siberia)',
                                      locality_normalized),
         locality_normalized = ifelse(locality_normalized == '(pond (Krotovo Lake, Troitskoe, Novosibirsk, Western Siberia)',
                                      'Russia (Krotovo Lake, Troitskoe, Novosibirsk region, Western Siberia)',
                                      locality_normalized)) %>%
  rowwise() %>%
  mutate(tp = get_locality_tp(pred_locality, locality_normalized),
         fp = get_locality_fp(pred_locality, locality_normalized),
         fn = get_locality_fn(pred_locality, locality_normalized)) %>%
  mutate(regions_in_abstract = check_if_regions_in_abstract(title_abstract, locality_normalized))

# 18% precision, 45% recall
# if excluding cases where regions not in abstract, 24% precision and 77% recall
locality_precision <-
  sum(microsp_locality_preds$tp) / (sum(microsp_locality_preds$tp) + sum(microsp_locality_preds$fp))

locality_recall <-
  sum(microsp_locality_preds$tp) / (sum(microsp_locality_preds$tp) + sum(microsp_locality_preds$fn))