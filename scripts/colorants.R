library(readxl)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(tidyverse)
library(grid)
library(reshape2)
library(colorspace)
library(stats)

landscape_size <- 8



### AUX. FUNCTIONS

# binary representation of integer n (returns a string)
binary <- function(n) sapply(n, FUN = function(n) {
  binary_str <- ""
  while (n > 0) {
    if (n %% 2 == 1) {
      binary_str <- paste("1", binary_str, sep="")
    } else {
      binary_str <- paste("0", binary_str, sep="")
    }
    n <- bitwShiftR(n, 1)
  }
  if (binary_str == "") {
    binary_str <- "0"
  }
  return(binary_str)
})

# add leading zeros to binary number (as string)
add_leading_zeros <- function(n, length.out) as.character(sapply(n, FUN = function(n) {
  if (nchar(n) >= length.out) {
    return(n)
  } else {
    num_zeros <- length.out - nchar(n)
    return(paste(c(rep("0", num_zeros), substr(n, 1, length.out - num_zeros)), collapse = ''))
  }
}))

# binary to integer
binary_to_integer <- function(binary_str) sapply(binary_str, FUN = function(binary_str) { # returns n given its binary representation (as a string) -- this has a precision limit of 2^30 (will only work for landscapes of 30 species or less, but that protocol would require over 1e7 96-well plates anyway so we will never have experiments that big)
  n <- 0
  for (i in 1:nchar(binary_str)) {
    bit <- substr(binary_str, i, i)
    if (bit == "1") {
      n <- n + 2^(nchar(binary_str)-i)
    }
  }
  return(n)
})




### LOAD DATA

data <- read.table('../data/colorants.txt',
                   header = T,
                   colClasses = c('character', 'numeric', 'numeric'))

data$nspecies <- sapply(data$community,
                        FUN = function(x) sum(strsplit(x, split = '')[[1]] == 1))

# plot spectra of monocultures
ggplot(data[data$nspecies %in% c(0, 1), ], aes(x = wavelength, y = absorbance, group = community, color = community)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c('black', '#e76d57', '#86a8db', '#22b37b', '#a541cc', '#4446d4', '#acd152', '#e7da5c', '#b38900')) +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank())



### PLOTS

# plot layout
design <- matrix(NA, nrow = landscape_size + 1,
                 ncol = max(c(choose(landscape_size, ceiling(landscape_size/2)), choose(landscape_size, floor(landscape_size/2)))))
panel_index <- 1
for (i in 0:landscape_size) {
  n_panels <- choose(landscape_size, i)
  design[i + 1, 1:n_panels] <- panel_index:(panel_index + n_panels - 1)
  panel_index <- panel_index + n_panels
}

# sort communities by number of species & add aux. variables
data$community_id <- binary_to_integer(data$community)
data <- data[order(data$nspecies, data$community_id, decreasing = F), ]
data$community <- factor(data$community, levels = unique(data$community))
data$spec <- 'co-culture'


### PLOT

# base plot
myplot <-
  ggplot(data,
         aes(x = wavelength, y = absorbance, group = spec, color = spec)) +
  geom_blank() +
  # geom_line(linewidth = 0.25,
  #           color = 'black') +
  scale_x_continuous(breaks = c(500, 700),
                     name = 'Wavelength (nm)') +
  scale_y_continuous(breaks = pretty_breaks(n = 2),
                     name = 'Absorbance') +
  facet_manual(vars(community),
               design = design,
               widths = 1,
               heights = 0.3) +
  theme_bw() + 
  theme(aspect.ratio = 0.3,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

# plot monocultures

monocs_sum_all <- data.frame(wavelength = numeric(0),
                             community = character(0),
                             nspecies = numeric(0),
                             absorbance = numeric(0),
                             spec = character(0))

for (comm in unique(data$community)) {
  
  monocs <- which(strsplit(as.character(comm), split = '')[[1]] == '1')
  if(length(monocs) >= 1) {
    
    monocs <- sapply(1:length(monocs),
                     FUN = function(i) {
                       monoc_id <- rep('0', landscape_size)
                       monoc_id[monocs[i]] <- '1'
                       monoc_id <- paste(monoc_id, collapse = '')
                       return(monoc_id)
                     })
    
    monocs_spect <- do.call(rbind,
                            lapply(monocs,
                                   FUN = function(m) {
                                     df <- data[data$community == m, ]
                                     df$spec <- m
                                     return(df)
                                   }))
    monocs_spect$community <- comm
    monocs_spect$community <- factor(monocs_spect$community,
                                     levels = levels(data$community))
    
    myplot <- myplot +
      geom_line(data = monocs_spect,
                linewidth = 0.25,
                linetype = 'dashed')
    
    # sum of monoc. spectra
    if(length(monocs) > 1) {
      
      monocs_sum <- aggregate(absorbance ~ wavelength + community + nspecies,
                              data = monocs_spect,
                              FUN = sum)
      monocs_sum$spec <- 'sum'
      
      myplot <- myplot +
        geom_line(data = monocs_sum,
                  linewidth = 0.5,
                  linetype = 'dashed')
      
      monocs_sum_all <- rbind(monocs_sum_all,
                              monocs_sum)
      
    }
    
  }
  
}

# plot consortia spectra
for (comm in unique(data$community)) {
  if(sum(strsplit(comm, split = '')[[1]] == '1') != 1) {
    
    myplot <- myplot +
      geom_line(data = data[data$community == comm, ],
                linewidth = 0.5)
    
  }
}

mycolors <- c('#E76D57', '#ACD152', '#7A9DD2', '#F7AFBE', '#E5C91F', '#A88BC0', '#5CA985', '#BBAD75')

myplot <- myplot +
  scale_color_manual(name = 'Species',
                     values = c(mycolors, 'black', 'black'))

ggsave(myplot,
       filename = paste('../plots/colorants_spectra.pdf', sep = ''),
       width = 1500,
       height = 200,
       units = 'mm',
       limitsize = F)














# absolute and relative errors
plot_this <- merge(data[, c('community', 'wavelength', 'absorbance')],
                   monocs_sum_all[, c('community', 'wavelength', 'absorbance')],
                   by = c('community', 'wavelength'),
                   suffixes = c('.obs', '.sum'))


plot_this$absErr <- plot_this$absorbance.sum - plot_this$absorbance.obs
plot_this$relErr <- abs(plot_this$absErr)/plot_this$absorbance.obs

plot_this$xpos <- rnorm(nrow(plot_this), 0, 0.1)

# relative error
ggplot(plot_this[plot_this$absorbance.obs > 0.1, ], aes(x = xpos, y = 100*relErr)) +
  # geom_point(alpha = 0.01,
  #            cex = 0.5,
  #            color = 'gray') +
  geom_bin2d(binwidth = 0.5*c(1/80, 1)) +
  geom_boxplot(x = 0, outlier.shape = NA, width = 0.2) +
  scale_fill_gradient(low = 'gray90', high = 'black') +
  scale_y_continuous(limits = c(0, 100),
                     name = 'Relative deviation between\nexpected and empirical absorbances (%)') +
  #scale_x_continuous(limits = c(-0.4, 0.4)) +
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5) +
  guides(fill = 'none')

ggsave(filename = paste('../plots/colorants_relErr_all.pdf', sep = ''),
       width = 50,
       height = 100,
       units = 'mm',
       limitsize = F)

# same but only for consortium 00000111

ggplot(plot_this[plot_this$absorbance.obs > 0.1 & plot_this$community == '10000101', ],
       aes(x = xpos, y = 100*relErr)) +
  # geom_point(alpha = 0.01,
  #            cex = 0.5,
  #            color = 'gray') +
  geom_bin2d(binwidth = 0.5*c(1/80, 1)) +
  geom_boxplot(x = 0, outlier.shape = NA, width = 0.2) +
  scale_fill_gradient(low = 'gray90', high = 'black') +
  scale_y_continuous(limits = c(0, 100),
                     name = 'Relative deviation between\nexpected and empirical absorbances (%)') +
  #scale_x_continuous(limits = c(-0.4, 0.4)) +
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5) +
  guides(fill = 'none')

ggsave(filename = paste('../plots/colorants_relErr_sample.pdf', sep = ''),
       width = 50,
       height = 100,
       units = 'mm',
       limitsize = F)


# relative error vs consortium size
plot_this$nspecies <- sapply(as.character(plot_this$community),
                             FUN = function(x) sum(strsplit(x, split = '')[[1]] == '1'))
ggplot(plot_this[plot_this$absorbance.obs > 0.1 & plot_this$nspecies > 1, ], aes(x = nspecies + xpos/2, y = 100*relErr, group = nspecies)) +
  # geom_point(alpha = 0.01,
  #            cex = 0.5,
  #            color = 'gray') +
  geom_bin2d(binwidth = 0.5*c(1/80, 1)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) +
  scale_fill_gradient(low = 'gray90', high = 'black') +
  scale_y_continuous(name = 'Relative deviation between\nexpected and empirical\nabsorbances (%)') +
  scale_x_continuous(name = 'Number of colorants in consortium',
                     breaks = 2:8) +
  theme_bw() + 
  theme(aspect.ratio = 0.5,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        panel.border = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5) +
  guides(fill = 'none')

ggsave(filename = paste('../plots/colorants_relErr_vs_size.pdf', sep = ''),
       width = 150,
       height = 75,
       units = 'mm',
       limitsize = F)











# absolute error
ggplot(plot_this[plot_this$absorbance.obs > 0.1, ], aes(x = xpos, y = abs(absErr))) +
  # geom_point(alpha = 0.01,
  #            cex = 0.5,
  #            color = 'gray') +
  geom_bin2d(binwidth = 0.5*c(1/80, 1/300)) +
  geom_boxplot(x = 0, outlier.shape = NA, width = 0.2) +
  scale_fill_gradient(low = 'gray90', high = 'black') +
  scale_y_continuous(limits = c(0, 0.5),
                     name = 'Absolute deviation between\nexpected and empirical absorbances (A.U.)') +
  #scale_x_continuous(limits = c(-0.4, 0.4)) +
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5) +
  guides(fill = 'none')

ggsave(filename = paste('../plots/colorants_absErr_all.pdf', sep = ''),
       width = 50,
       height = 100,
       units = 'mm',
       limitsize = F)


if (F) {
  # max. absolute error across wavelengths
  plot_this <- aggregate(absErr ~ community,
                         data = plot_this,
                         FUN = function(x) max(abs(x)))
  plot_this$xpos <- rnorm(nrow(plot_this), 0, 0.1)
  
  ggplot(plot_this, aes(x = xpos, y = absErr)) +
    # geom_point(alpha = 0.01,
    #            cex = 0.5,
    #            color = 'gray') +
    #geom_bin2d(binwidth = 0.5*c(1/80, 1/300)) +
    geom_point(alpha = 0.25) +
    geom_boxplot(x = 0, outlier.shape = NA, width = 0.2) +
    scale_y_continuous(limits = c(0, 0.5),
                       name = 'Max. absolute deviation between\nexpected and empirical absorbances (A.U.)') +
    scale_x_continuous(limits = c(-0.4, 0.4)) +
    theme_bw() + 
    theme(aspect.ratio = 3,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(vjust = 0),
          panel.border = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 13)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5) +
    guides(fill = 'none')
  
  ggsave(filename = paste('../plots/colorants_maxAbsErr_all.pdf', sep = ''),
         width = 50,
         height = 100,
         units = 'mm',
         limitsize = F)
}




# is there an effect of wavelength on deviations?
p <- quantile(plot_this$wavelength, probs = seq(0, 1, length.out = 10))
plot_this$quantile <- sapply(plot_this$wavelength,
                             FUN = function(x) which.min(abs(x - p)))
median_wls <- aggregate(wavelength ~ quantile,
                        data = plot_this,
                        FUN = median)
plot_this$quantile_wl <- median_wls$wavelength[plot_this$quantile]
ggplot(plot_this[plot_this$absorbance.obs > 0.1, ], aes(x = wavelength, y = 100*relErr, group = quantile)) +
  # geom_point(alpha = 0.01,
  #            cex = 0.5,
  #            color = 'gray') +
  #geom_bin2d(binwidth = c(1, 1/2)) +
  #geom_point(alpha = 0.25) +
  geom_boxplot(aes(x = quantile_wl, y = 100*relErr),
               outlier.shape = NA, width = 20) +
  scale_fill_gradient(low = 'gray90', high = 'black') +
  scale_y_continuous(limits = c(0, 45),
                     name = 'Relative deviation between\nexpected and empirical absorbances (%)') +
  scale_x_continuous(name = 'Wavelength (nm)') +
  theme_bw() + 
  theme(aspect.ratio = 0.5,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(vjust = 0),
        panel.border = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5) +
  guides(fill = 'none')

ggsave(filename = paste('../plots/colorants_relErr_vs_wavelength.pdf', sep = ''),
       width = 150,
       height = 75,
       units = 'mm',
       limitsize = F)


