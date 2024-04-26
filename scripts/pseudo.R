library(readxl)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(tidyverse)
library(grid)
library(reshape2)
library(colorspace)

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

data <- read.table('../data/pseudo.txt', header = T, colClasses = c('character', 'character', 'numeric', 'numeric', 'numeric'))

data$nspecies <- sapply(data$community,
                        FUN = function(x) sum(strsplit(x, split = '')[[1]] == 1))

# filter out highest dil. factors and aggregate replicates
data <- data[data$dilution_factor == 0.0025, ]
data <- do.call(data.frame,
                aggregate(absorbance ~ community + nspecies + wavelength,
                          data,
                          FUN = function(x) c(mean = mean(x), sd = sd(x))))
colnames(data)[4:5] <- c('absorbance', 'sd')




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
         aes(x = wavelength, y = absorbance,
             ymin = absorbance - sd, ymax = absorbance + sd,
             group = spec, color = spec, fill = spec)) +
  geom_blank() +
  # geom_line(linewidth = 0.25,
  #           color = 'black') +
  scale_x_continuous(breaks = c(500, 700),
                     name = 'Wavelength (nm)') +
  # scale_y_continuous(breaks = pretty_breaks(n = 2),
  #                    name = 'Absorbance') +
  facet_manual(vars(community),
               design = design,
               widths = 1,
               heights = 0.3) +
  theme_bw() + 
  # theme(aspect.ratio = 0.3,
  #       panel.grid = element_blank(),
  #       strip.background = element_blank(),
  #       legend.position = 'none',
  #       panel.border = element_blank()) +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        strip.text = element_text(size = 8), #strip.text = element_blank(),
        #axis.ticks.length = unit(0.25, "mm"),
        #axis.ticks = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16))
  # annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  # annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)



# plot monocultures

monocs_sum_all <- data.frame(wavelength = numeric(0),
                             community = character(0),
                             nspecies = numeric(0),
                             absorbance = numeric(0),
                             spec = character(0))

for (comm in unique(data$community)) {
  
  monocs <- which(strsplit(as.character(comm), split = '')[[1]] == '1')
  if(length(monocs) > 1) {
    
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
      monocs_sum$sd <- NA
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
  if(sum(strsplit(as.character(comm), split = '')[[1]] == '1') != 1) {
    
    myplot <- myplot +
      geom_ribbon(data = data[data$community == comm, ],
                  color = NA,
                  alpha = 0.2) +
      geom_line(data = data[data$community == comm, ],
                linewidth = 0.5)
    
  }
}

# plot monocultures
singles <- data[data$nspecies == 1, ]
singles$spec <- singles$community

myplot <- myplot +
  geom_ribbon(data = singles,
              color = NA,
              alpha = 0.2) +
  geom_line(data = singles,
            linewidth = 0.5)

# adjust colors and scales for nice aesthetics
mycolors <- c('#E76D57', '#ACD152', '#7A9DD2', '#F7AFBE', '#E5C91F', '#A88BC0', '#5CA985', '#BBAD75')
mybreaks <- c(seq(0.01, 0.1, length.out = 10), seq(0.2, 1, length.out = 9), seq(2, 10, length.out = 9), seq(20, 100, length.out = 9))

myplot <- myplot +
  scale_color_manual(name = 'Species',
                     values = c(mycolors, 'black', 'black')) +
  scale_fill_manual(name = 'Species',
                     values = c(mycolors, 'black', 'black')) +
  scale_y_log10(name = 'Absorbance (A.U.)',
                breaks = mybreaks)

ggsave(myplot,
       filename = paste('../plots/pseudo_spectra_clean.pdf', sep = ''),
       width = 1500,
       height = 200,
       units = 'mm',
       limitsize = F)


# myplot <- myplot +
#   scale_y_log10(breaks = c(0.1, 1, 10),
#                 name = 'Absorbance')
# 
# ggsave(myplot,
#        filename = paste('../plots/pseudo_spectra_log.pdf', sep = ''),
#        width = 1500,
#        height = 200,
#        units = 'mm',
#        limitsize = F)


### small plot for main fig
# data <- data[order(data$community_id), ]
# data$community <- factor(data$community,
#                          levels = unique(data$community))
# myplot_main <- 
#   ggplot(data,
#          aes(x = wavelength, y = absorbance, ymin = absorbance - sd, ymax = absorbance + sd)) +
#     geom_ribbon(color = NA,
#                 fill = 'black',
#                 alpha = 0.2) +
#     geom_line(color = 'black') +
#     # facet_wrap(~ community,
#     #            nrow = 8,
#     #            dir = 'v') +
#     facet_wrap2(~ community,
#                 nrow = 8,
#                 dir = 'v',
#                 axes = 'all',
#                 remove_labels = 'all') +
#     # scale_y_continuous(name = 'Absorbance (A.U.)',
#     #                    breaks = pretty_breaks(n = 2)) +
#     scale_y_log10(name = 'Absorbance (A.U.)',
#                   breaks = seq(0.1, 1, length.out = 10)) +
#     scale_x_continuous(name = 'Wavelength (nm)',
#                        breaks = c(500, 700)) +
#     theme_bw() + 
#     theme(aspect.ratio = 1,
#           panel.grid = element_blank(),
#           strip.background = element_blank(),
#           legend.position = 'none',
#           panel.border = element_blank(),
#           axis.line.x = element_line(),
#           axis.line.y = element_line(),
#           axis.text = element_text(size = 13),
#           axis.title = element_text(size = 16),
#           #strip.text = element_text(size = 8),
#           axis.ticks.length = unit(0.25, "mm"),
#           #axis.ticks = element_blank(),
#           strip.text = element_blank())
#     # annotate("segment", x=-Inf, xend=Inf, y=0.09, yend=0.09, linewidth=0.5) +
#     # annotate("segment", x=-Inf, xend=-Inf, y=0.09, yend=Inf, linewidth=0.5)
# 
# ggsave(myplot_main,
#        filename = paste('../plots/pseudo_spectra_main_log_ticks.pdf', sep = ''),
#        width = 300,
#        height = 150,
#        units = 'mm',
#        limitsize = F)


### QUANTIFY INTERACTIONS

# # get interaction coefficients at all orders
# focal_wavelength <- 600
# data <- data[data$wavelength == focal_wavelength, c('community', 'absorbance')]
# df <- do.call(rbind,
#               lapply(1:nrow(data),
#                      FUN = function(i) as.data.frame(t(strsplit(as.character(data$community[i]), split = '')[[1]]))))
# df$y <- data$absorbance
# colnames(df) <- c(paste('x', 1:landscape_size), 'y')
# 
# # which base?
# base <- 'taylor' # one of 'taylor' or 'fourier'
# if (base == 'fourier') {
#   for (i in 1:landscape_size) df[df[, i] == 0, i] <- -1
# }
# 
# # model <- lm(y ~ . + .*. + .*.*. + .*.*.*. + .*.*.*.*. + .*.*.*.*.*. + .*.*.*.*.*.*. + .*.*.*.*.*.*.*.,
# #             data = df)
# model <- lm(y ~ .*.*.*.*.*.*.*.,
#             data = df)
# coef <- model$coefficients
# names(coef) <- gsub('`1', '', names(coef))
# names(coef) <- gsub('`x', '', names(coef))
# names(coef) <- gsub(' ', '', names(coef))
# names(coef) <- gsub(':', '', names(coef))
# names(coef)[1] <- ''
# 
# coef_df <- data.frame(coef_id = names(coef),
#                       coef_order = nchar(names(coef)),
#                       value = as.numeric(coef))
# ggplot(coef_df, aes(x = coef_order, y = value, group = coef_order)) +
#   geom_abline(slope = 0,
#               intercept = 0,
#               color = 'gray') +
#   geom_boxplot(data = coef_df[coef_df$coef_order > 0 & coef_df$coef_order < 8, ],
#                fill = 'white',
#                outlier.shape = NA) +
#   geom_jitter(width = 0.15,
#               alpha = 0.25) +
#   scale_x_continuous(name = 'Order',
#                      breaks = 0:8) +
#   scale_y_continuous(name = expression(paste('Coefficient value, ', italic(beta), ' (A.U.)', sep = ''))) +
#   theme_bw() + 
#   theme(aspect.ratio = 0.6,
#         panel.grid = element_blank(),
#         strip.background = element_blank(),
#         legend.position = 'none',
#         panel.border = element_blank(),
#         axis.text = element_text(size = 13),
#         axis.title = element_text(size = 16)) +
#   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
#   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)
# 
# ggsave(filename = paste('../plots/pseudo_coef_', base, '.pdf', sep = ''),
#        width = 100,
#        height = 100,
#        units = 'mm',
#        limitsize = F)
# 
# model_y <- predict(model, df[, 1:landscape_size])





