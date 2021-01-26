library(dplyr)
library(ggsci)
library(ggplot2)
library(hash)
library(magrittr)
library(data.table)

### parse a name table
full_default <- read.delim('DIGE_full.csv', sep = '\t', header = T, stringsAsFactors = F) %>% mutate(V = coalesce(V, 0))

### correct the redundant annotations
full_default[full_default$Protein.ID == 'P23382', 'Protein.ID'] <- 'WP_023521456.1'
full_default[full_default$Protein.ID == 'WP_098301448.1', 'Protein.ID'] <- 'WP_023521456.1'
full_default[full_default$Protein.ID == 'WP_050842914.1', 'Protein.ID'] <- 'WP_023521456.1'
full_default[full_default$Protein.ID == 'WP_016117439.1', 'Protein.ID'] <- 'WP_023521456.1'
full_default[full_default$Protein.ID == 'AXE15630.1', 'Protein.ID'] <- 'P0A382'
full_default[full_default$Protein.ID == 'B8I1A8', 'Protein.ID'] <- 'WP_065703111.1'

full_default[full_default$Protein.ID == 'EEM29970.1', 'Protein.ID'] <- 'WP_011109926.1'
full_default[full_default$Protein.ID == 'E5L4Y0', 'Protein.ID'] <- 'WP_011109926.1'
full_default[full_default$Protein.ID == 'A0A243MFA1', 'Protein.ID'] <- 'WP_000156601.1'

### Remove the conspicuous, Cry-annotated spots from further analysis
full_default %<>% filter(!(Protein.ID %in% c('Q9ZIU5', 'CBX51930.1')))

### closure the duplicate entries
full <- full_default %>% aggregate(. ~ Protein.ID + Gel, ., any) %>% mutate(across(where(is.logical), as.numeric))

### write down the deduplicated name list
fwrite(list(full$Protein.ID) %>% unique, 'unique_DIGE_names.txt')

### parse an Emapper output
annot <- read.delim(file.path('annots', 'full', 'full.emapper.annotations'), sep = '\t', header = T, skip = 3, stringsAsFactors = F)

### add COG categories to the name table
full$COG <- sapply(c(1:nrow(full)), function(x) annot[annot$X.query_name == full$Protein.ID[x], 'COG.Functional.cat.'] %>% unlist()) %>%
  lapply(function(x) if(identical(x, character(0))) 'S' else x) %>% lapply(function(x) x %>% strsplit('') %>% unlist() %>% paste(collapse = ';')) %>% unlist()

### create a COG pivot in a sample- and strain-wise manner

subset_DIGE <- function(strain=NULL, gel=NULL, ful=full) {

    if (!is.null(strain)) {
    ful <- ful[ful[[strain]] == 1,]
  }
  
  if (!is.null(gel)) {
    ful <- ful[ful$Gel == gel,]
  }
  
  return(ful)
}

COG_table <- function(df) {
  cogs <- df[['COG']] %>% sapply(function(x) strsplit(x, ';') %>% unlist()) %>% unlist() %>% unname()
  return(cogs)
}

COG_table(subset_DIGE(strain = 'I')) %>% table()
COG_table(subset_DIGE(strain = 'T')) %>% table()
COG_table(subset_DIGE(strain = 'D')) %>% table()
COG_table(subset_DIGE(strain = 'V')) %>% table()
COG_table(subset_DIGE(strain = 'A')) %>% table()

# COG_table(subset_DIGE(strain = 'V', gel = 2)) %>% table()
# COG_table(subset_DIGE(strain = 'A', gel = 2)) %>% table()

### Map COG identifiers to their meanings

COG_names <- read.delim('COG_names.tsv', sep = '\t', header = F, stringsAsFactors = F)
COG_names <- hash(COG_names[,1], COG_names[,2])


### convert COG tables into DFs

COG_df <- function(df, strain) {
  df %<>% table() %>% data.frame() %>% mutate(Strain = strain)
  return(df)
}

### set a uniform color scheme throughout the paper

paper_colors <- c('#B772A1', '#A1B772', '#6B8E23', '#72A1B7', '#A60B0B', '#2980B9', '#727EB7', '#D95252', '#D6D626')

### Build COG distribution histograms

splitLabels <- function(x, aver=NA) {
  if (is.na(aver)) {
  aver <- mean(sapply(x, function(y) nchar(y))) %>% round()
  }
  splitAfter <- function(z, marg) {
    occs <- grep('( )|(, )', strsplit(z, '')[[1]])
    dists <- sapply(occs[occs <= marg], function(x) x - marg)
    margind <- occs[occs <= marg][which(dists == max(dists))]
    startind <- ifelse(substr(z, margind + 1, margind + 1) == ' ', 2, 1)
    new_str <- paste(c(substr(z, 1, margind), '\n', substr(z, margind + startind, nchar(z))), collapse = '')
    return(new_str)
  }
  outp <- sapply(x, function(y) ifelse(nchar(y) > aver, splitAfter(y, aver), y)) %>% unname()
  return(outp)
}

COG_plotter <- function(...) {
  df <- do.call(rbind, list(...)) %>% `colnames<-`(c('COG', 'Frequency', 'Strain'))
  df$COG <- factor(df$COG, levels = c(sort(levels(df$COG)[levels(df$COG) != 'S']), 'S'))
  output <- ggplot(df, aes(x = Strain, y = Frequency)) +
    geom_bar(aes(fill = reorder(COG, -Frequency)), position = position_dodge(preserve = 'single'), stat = 'identity') +
    theme_bw() +
    theme(axis.title.x = element_text(
                         size = 30,
                         color = 'black'
                         ),
          axis.text.x = element_text(
                        size = 28,
                        color = 'black',
                        angle = 45,
                        vjust = 0.4,
                        hjust = 0.3
                        ),
          axis.title.y = element_text(
            size = 30,
            color = 'black'
            ),
          axis.text.y = element_text(
            size = 28,
            color = 'black'
            ),
          panel.background = element_rect(),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(
            size = 0.3,
            linetype = c('28'), 
            colour = "black"
            ),
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 28)) +
    scale_fill_manual(name = 'COG Term' ,
                      labels = sapply(levels(reorder(df$COG, -df$Frequency)),
                                      function(x) values(COG_names, key = x)) %>% 
                               unname() %>% 
                               splitLabels(),
                      values = paper_colors) +
    scale_color_discrete(guide = FALSE)
  return(output)
}

### Wrap the previous functions into one master-function

# COG_analysis <- function(strain=NULL, gel=NULL, ful=full) {
#   
# }

  


COG_plotter(COG_df(COG_table(subset_DIGE(strain = 'I')), 'I') %>% mutate(Strain = rep('800/3 (S)', 4)),
            COG_df(COG_table(subset_DIGE(strain = 'A')), 'A') %>% mutate(Strain = rep('800/3-15 (S)', 4)),
            COG_df(COG_table(subset_DIGE(strain = 'V')), 'V') %>% mutate(Strain = rep('800/3 (V)', 7)))

COG_plotter(COG_df(COG_table(subset_DIGE(strain = 'I')), 'I') %>% mutate(Strain = rep('800/3', 4)),
            COG_df(COG_table(subset_DIGE(strain = 'D')), 'D') %>% mutate(Strain = rep('109/25', 3)),
            COG_df(COG_table(subset_DIGE(strain = 'T')), 'T') %>% mutate(Strain = rep('800/15', 4)))
