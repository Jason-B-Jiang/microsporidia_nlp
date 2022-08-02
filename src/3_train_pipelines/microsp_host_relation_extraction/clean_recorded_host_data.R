# -----------------------------------------------------------------------------
#
# Clean recorded microsporidia species + host names data
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/08/02
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

library(tidyverse)

################################################################################

main <- function() {
  species_hosts <- read_csv('../../../data/manually_format_multi_species_papers.csv') %>%
    select(species, num_papers, paper_title, title_abstract, hosts_natural,
           hosts_experimental) %>%
    mutate(hosts_natural = clean_species_names(hosts_natural),
           hosts_experimental = clean_species_names(hosts_experimental),
           species_cleaned = clean_species_names(species),
           all_hosts = combine_all_hosts(hosts_natural, hosts_experimental)) %>%
    rowwise() %>%
    mutate(hosts_not_in_text = get_species_not_in_text(all_hosts, title_abstract),
           microsp_not_in_text = get_species_not_in_text(species, title_abstract)) %>%
    ungroup() %>%
    mutate(hosts_not_in_text_corrected = NA,
           microsp_not_in_text_corrected = NA)
  
  to_manually_correct <- species_hosts %>%
    filter(!is.na(hosts_not_in_text) | !is.na(microsp_not_in_text))
  
  if (!file.exists('../../../data/manually_correct_recorded_hosts.csv')) {
    write_csv(to_manually_correct, '../../../data/manually_correct_recorded_hosts.csv')
    stop("Please manually annotate microsporidia and host names in csv outputted")
    
  } else {
    manually_corrected <- read_csv('../../../data/manually_correct_recorded_hosts.csv')
    
    if (all(is.na(manually_corrected$hosts_not_in_text_corrected)) &
        all(is.na(manually_corrected$microsp_not_in_text_corrected))) {
      stop("Microsporidia and host names in texts were not manually corrected")
      
    } else {
      tmp <- species_hosts %>%
        filter(!(species %in% manually_corrected$species)) %>%
        bind_rows(manually_corrected) %>%
        mutate(hosts_in_text = resolve_species(all_hosts, hosts_not_in_text,
                                               hosts_not_in_text_corrected),
               microsp_in_text = resolve_species(species_cleaned, microsp_not_in_text,
                                                 microsp_not_in_text_corrected)) %>%
        filter(num_papers == 1, !is.na(title_abstract),
               !is.na(hosts_in_text) | !is.na(microsp_in_text)) %>%
        select(species, title_abstract, species_cleaned, all_hosts, hosts_in_text,
               microsp_in_text) %>%
        mutate(all_hosts = ifelse(is.na(all_hosts), '', all_hosts),
               hosts_in_text = ifelse(is.na(hosts_in_text), '', hosts_in_text),
               microsp_in_text = ifelse(is.na(microsp_in_text), '', microsp_in_text)) %>%
        group_by(title_abstract) %>%
        mutate(species = str_c(species, collapse = ' || '),
               species_cleaned = str_c(species_cleaned, collapse = ' || '),
               all_hosts = str_c(all_hosts, collapse = ' || '),
               hosts_in_text = str_c(hosts_in_text, collapse = ' || '),
               microsp_in_text = str_c(microsp_in_text, collapse = ' || ')) %>%
        distinct(.keep_all = T)
      
      write_csv(tmp, '../../../data/formatted_species_host_data.csv')
    }
  }
}

################################################################################

## Helper functions

clean_species_names <- Vectorize(function(hosts) {
  hosts <- str_remove(str_split(hosts, '; ')[[1]], ' ?\\(.+\\)')
  return(str_c(str_remove(hosts, ' \\d+$'), collapse = '; '))
})


combine_all_hosts <- Vectorize(function(hosts_natural, hosts_experimental) {
  if (is.na(hosts_natural) & is.na(hosts_experimental)) {
    return(NA)
  } else if (is.na(hosts_natural)) {
    return(hosts_experimental)
  } else if (is.na(hosts_experimental)) {
    return(hosts_natural)
  }
  
  return(str_c(hosts_natural, hosts_experimental, sep = '; '))
}, vectorize.args = c('hosts_natural', 'hosts_experimental'))


get_abbreviated_species_name <- function(species) {
  species <- str_split(species, ' ')[[1]]
  
  if (length(species) == 1) {
    return(species)
  }
  
  # ugly string concat operation
  str_c(
    str_c(
      sapply(species[1 : length(species) - 1], function(s) {str_c(substr(s, 1, 1), '.')}),
      collapse = ' '),
    species[length(species)], sep = ' '
    )
}


get_species_not_in_text <- function(species, txt) {
  if (is.na(txt) | is.na(species)) {
    return(NA)
  }
  
  species <- sapply(str_split(species, '; ')[[1]], function(s) {clean_species_names(s)})
  species_not_in_text <- character()
  
  for (sp in species) {
    if (!str_detect(tolower(txt), tolower(fixed(sp))) & 
        !str_detect(tolower(txt), tolower(fixed(get_abbreviated_species_name(sp))))) {
      species_not_in_text <- c(species_not_in_text, sp)
    }
  }
  
  if (length(species_not_in_text) == 0) {
    return(NA)
  }
  
  return(str_c(species_not_in_text, collapse = '; '))
}


resolve_species <- Vectorize(function(all_species, not_in_text, corrected) {
  if (is.na(all_species)) {
    return(NA)
  }
  
  # in case we are looking at microsporidia species data
  all_species <- str_split(clean_species_names(all_species), '; ')[[1]]
  not_in_text <- str_split(clean_species_names(not_in_text), '; ')[[1]]
  corrected <- str_split(corrected, '; ')[[1]]
    
  if (all(is.na(not_in_text))) {
    return(str_c(all_species, collapse = '; '))
    
    } else if (all(is.na(corrected))) {
      tmp = str_c(all_species[!(all_species %in% not_in_text)], collapse = '; ')
    return(ifelse(tmp == '', NA, tmp))
  }
  
  return(
    str_c(c(all_species[!(all_species %in% not_in_text)], corrected), collapse = '; ')
  )
}, vectorize.args = c('all_species', 'not_in_text', 'corrected'))


get_match_spans <- function(subtext, text) {
  # note: string searches are case-insensitive
  matches = str_locate_all(tolower(text), fixed(tolower(subtext)))[[1]]
  match_spans <- c()
  
  if (nrow(matches) == 0) {
    return(character())
  }
  
  for (i in 1 : nrow(matches)) {
    match_spans <- c(match_spans, str_c(c(matches[i,][1] - 1, matches[i,][2]),
                                        collapse = '-'))
  }
  
  return(str_c(match_spans, collapse = '; '))
}


get_spans_from_text <- function(entities, text) {
  if (is.na(entities) | is.na(text)) {
    return(NA)
  }
  
  entities <- str_split(entities, '; ')[[1]]
  spans <- c()
  
  for (ent in entities) {
    abbrev <- get_abbreviated_species_name(ent)
    
    if (abbrev == ent) {
      spans <- c(spans, get_match_spans(ent, text))
    } else {
      spans <- c(spans, str_c(get_match_spans(ent, text), get_match_spans(abbrev, text), sep = '; '))
    }
  }
  
  return(str_c(spans, collapse = ' | '))
}

################################################################################

# Note to self:
# Microsporidia of the genus Amblyospora parasiting the adipose body of mosquito larvae of the genus Aedes and Culex has been studied
# example of host species names being implicitly given
#
# First report of microsporidian infections in solefishes from Senegal coast (West Africa):
# didn't record 2 species + host pairs actually mentioned in paper

main()