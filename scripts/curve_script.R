library(dplyr)
library(ggplot2)
library(lattice)

df <- read.delim('growth_curves_edited.csv', sep='\t', header=T, stringsAsFactors = F)
df %>% View()

strains <- df$Sample %>% unique() %>% .[matches('^((?!Blank))', perl = T, vars=.)]
media <- df$Group %>% unique() %>% .[matches('^((?!E))', perl = T, vars = .)]

timeConvert <- function(x) {
  x <- stringr::str_split(x, 'h') %>% unlist()
  h <- x[1] %>% as.numeric()
  if (x[2] == '') return(h) 
  else {
    m <- x[2] %>% stringr::str_split('min') %>% unlist() %>% dplyr::first() %>% as.numeric()
    return(h + m / 60)
  }
}

### Figure 1

I <- df %>% filter(Sample == 'Sample X3', Group == 'B')
Im <- I[, -c(1:3)] %>% as.matrix() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  purrrlyr::by_row(..f = function(x) {x <- sd(unlist(x)[2:5]) / sqrt(length(x[2:5]))}) %>%
  transmute(Time = rowname %>% regexPipes::gsub('X', '') %>% regexPipes::gsub('\\.', '') %>% sapply(timeConvert),
    Mean = purrr::pmap_dbl(list(V1, V2, V3, V4), mean), Se = .out %>% unlist()) %>% mutate(upper = Mean + Se, lower = Mean - Se)

A <- df %>% filter(Sample == 'Sample X1', Group == 'B')
Am <- A[, -c(1:3)] %>% as.matrix() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  purrrlyr::by_row(..f = function(x) {x <- sd(unlist(x)[2:5]) / sqrt(length(x[2:5]))}) %>%
  transmute(Time = rowname %>% regexPipes::gsub('X', '') %>% regexPipes::gsub('\\.', '') %>% sapply(timeConvert),
            Mean = purrr::pmap_dbl(list(V1, V2, V3, V4), mean), Se = .out %>% unlist()) %>% mutate(upper = Mean + Se, lower = Mean - Se)

V <- df %>% filter(Sample == 'Sample X3', Group == 'A')
Vm <- V[, -c(1:3)] %>% as.matrix() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  purrrlyr::by_row(..f = function(x) {x <- sd(unlist(x)[2:5]) / sqrt(length(x[2:5]))}) %>%
  transmute(Time = rowname %>% regexPipes::gsub('X', '') %>% regexPipes::gsub('\\.', '') %>% sapply(timeConvert),
            Mean = purrr::pmap_dbl(list(V1, V2, V3, V4), mean), Se = .out %>% unlist()) %>% mutate(upper = Mean + Se, lower = Mean - Se)

Fig1 <- rbind(Im, Am, Vm) %>% mutate(Culture = as.factor(c(rep('I', nrow(Im)), rep('A', nrow(Am)), rep('V', nrow(Vm)))))

ggplot(Fig1, aes(x = Time, y = Mean)) + 
  geom_line(size = 0.5, aes(color = Culture)) + 
  geom_point(size = 0.5, aes(color = Culture)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Culture,  width = 0.1), alpha = 0.75) +
  labs(x = 'Time, h', y = 'Mean absorbance') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30, color = 'black'),
        axis.text.x = element_text(size = 28, color = 'black'),
        axis.title.y = element_text(size = 30, color = 'black'),
        axis.text.y = element_text(size = 28, color = 'black'),
        panel.background = element_rect(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.2, linetype = c('28'),
                                        colour = "black"),
        legend.title = element_text(size = 30, color = 'black'),
        legend.text = element_text(size = 28, color = 'black')) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  scale_color_manual(name = 'Bacterial culture', 
                       labels = c('Strain 800/3-15, T3 medium',
                                  'Strain 800/3, T3 medium',
                                  'Strain 800/3, LB medium'),
                       values = c('#B772A1','#A60B0B','#2980B9'))

### Figure 2

D <- df %>% filter(Sample == 'Sample X2', Group == 'B')
Dm <- D[, -c(1:3)] %>% as.matrix() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  purrrlyr::by_row(..f = function(x) {x <- sd(unlist(x)[2:5]) / sqrt(length(x[2:5]))}) %>%
  transmute(Time = rowname %>% regexPipes::gsub('X', '') %>% regexPipes::gsub('\\.', '') %>% sapply(timeConvert),
            Mean = purrr::pmap_dbl(list(V1, V2, V3, V4), mean), Se = .out %>% unlist()) %>% mutate(upper = Mean + Se, lower = Mean - Se)

Tu <- df %>% filter(Sample == 'Sample X4', Group == 'B')
Tm <- Tu[, -c(1:3)] %>% as.matrix() %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% 
  purrrlyr::by_row(..f = function(x) {x <- sd(unlist(x)[2:5]) / sqrt(length(x[2:5]))}) %>%
  transmute(Time = rowname %>% regexPipes::gsub('X', '') %>% regexPipes::gsub('\\.', '') %>% sapply(timeConvert),
            Mean = purrr::pmap_dbl(list(V1, V2, V3, V4), mean), Se = .out %>% unlist()) %>% mutate(upper = Mean + Se, lower = Mean - Se)

Fig2 <- rbind(Im, Dm, Tm) %>% mutate(Strain = as.factor(c(rep('I', nrow(Im)), rep('D', nrow(Dm)), rep('T', nrow(Tm)))))

ggplot(Fig2, aes(x = Time, y = Mean)) + 
  geom_line(size = 0.5, aes(color = Strain)) + 
  geom_point(size = 0.5, aes(color = Strain)) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = Strain,  width = 0.1), alpha = 0.75) +
  labs(x = 'Time, h', y = 'Mean absorbance') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 30, color = 'black'),
        axis.text.x = element_text(size = 28, color = 'black'),
        axis.title.y = element_text(size = 30, color = 'black'),
        axis.text.y = element_text(size = 28, color = 'black'),
        panel.background = element_rect(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_line(size = 0.2, linetype = c('28'),
                                        colour = "black"),
        legend.title = element_text(size = 30, color = 'black'),
        legend.text = element_text(size = 28, color = 'black')) +
  guides(colour = guide_legend(override.aes = list(size = 1.5))) +
  scale_color_manual(name = 'Strain', 
                     labels = c('109/25',
                                '800/3',
                                '800/15'),
                     values = c('#B772A1','#A60B0B','#2980B9'))

