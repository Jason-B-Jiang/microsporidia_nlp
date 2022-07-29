# -----------------------------------------------------------------------------
#
# Clean recorded microsporidia host data
#
# Jason Jiang - Created: 2022/07/29
#               Last edited: 2022/07/29
#
# Mideo Lab - Microsporidia text mining
#
# -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)

################################################################################

main <- function() {
  microsp_df <- read_csv('../../data/manually_format_multi_species_papers.csv')
  
  # list mapping "canonical" species names back to original recorded
  # entries in dataset
  host_names <- make_host_names_list(
    read_xlsx('../../data/Table S2 Host Masterlist.xlsx')
  )
  
  to_manually_correct <- microsp_df %>%
    select(species, paper_title, title_abstract, hosts_natural, hosts_experimental) %>%
    mutate(hosts_natural = clean_host_names(hosts_natural),
           hosts_experimental = clean_host_names(hosts_experimental),
           hosts_not_in_text = get_hosts_not_in_text(hosts_natural,
                                                     hosts_experimental,
                                                     title_abstract),
           hosts_not_in_text_corrected = NA) %>%
    filter(!is.na(hosts_not_in_text))
  
  write_csv(to_manually_correct, '../../data/manually_correct_recorded_hosts.csv')
}

################################################################################

## Helper functions
make_host_names_list <- function(host_df) {
  host_names <- lapply(host_df$Host_original, function(x) {str_split(x, ',')[[1]]})
  names(host_names) <- host_df$Host_formatted
  
  return(host_names)
}


clean_host_names <- Vectorize(function(hosts) {
  hosts <- str_remove(str_split(hosts, '; ')[[1]], ' ?\\(.+\\)')
  return(str_c(hosts, collapse = '; '))
})


get_hosts_not_in_text <- Vectorize(function(hosts_natural, hosts_experimental, text) {
  hosts_natural <- str_split(hosts_natural, '; ')[[1]]
  hosts_experimental <- str_split(hosts_experimental, '; ')[[1]]
  
  not_in_text <- str_c(
    c(
    hosts_natural[which(!str_detect(text, fixed(hosts_natural)))],
    hosts_experimental[which(!str_detect(text, fixed(hosts_experimental)))]
    ),
  collapse = '; ', sep = '; '
  )
  
  if (not_in_text == '') {
    return(NA)
  }
  
  return(not_in_text)
}, vectorize.args = c('hosts_natural', 'hosts_experimental', 'text'))


hosts_not_in_host_names <- function(hosts_not_in_text, host_names) {
  hosts_not_in_text <- str_split(hosts_not_in_text, '; ')[[1]]
  
  missing <- c()
  for (h in hosts_not_in_text[which(!(hosts_not_in_text %in% names(host_names)))]) {
    if (!any(sapply(host_names, function(x) {h %in% x}))) {
      missing <- c(missing, h)
    }
  }
  
  if (length(missing) == 0) {
    return(NA)
  }
  
  return(str_c(missing, collapse = '; '))
}

################################################################################

# Note to self:
# Microsporidia of the genus Amblyospora parasiting the adipose body of mosquito larvae of the genus Aedes and Culex has been studied
# example of host species names being implicitly given

main()