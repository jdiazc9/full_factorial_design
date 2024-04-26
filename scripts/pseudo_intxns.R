library(readxl)
library(tidyr)
library(ggplot2)
library(ggh4x)
library(scales)
library(tidyverse)
library(grid)
library(reshape2)
library(colorspace)
library(combinat)
library(pbapply)
library(igraph)

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

# filter out highest dil. factors (we only use 1:400)
data <- data[data$dilution_factor == 0.0025, ]
data <- do.call(data.frame,
                aggregate(absorbance ~ community + nspecies + wavelength,
                          data,
                          FUN = function(x) c(mean = mean(x), sd = sd(x))))
colnames(data)[4:5] <- c('absorbance', 'sd')






### FUNCTIONAL LANDSCAPE
focal_wavelength <- 600

df <- data[data$wavelength == focal_wavelength, ]

# wrapper function: plot functional landscape
plotFitnessGraph <- function(df) {
  
  # make edges of fitness graph
  genots <- df$community
  
  # number of mutations in genotype
  nMut <- function(genot) sapply(genot,
                                 FUN = function(genot_i) sum(strsplit(genot_i, split = '')[[1]] == '1'))
  
  # check if a genotype is descendant from another
  isDescendant <- function(this_genot, of_this_genot) {
    
    this_genot <- which(strsplit(this_genot, split = '')[[1]] == '1')
    of_this_genot <- which(strsplit(of_this_genot, split = '')[[1]] == '1')
    
    return(all(of_this_genot %in% this_genot))
    
  }
  
  #make edges
  makeEdges <- function(genots) {
    
    edges <- data.frame(source = character(0),
                        target = character(0),
                        source.nmut = numeric(0),
                        target.nmut = numeric(0))
    
    for(s in genots) {
      
      t <- genots[sapply(genots,
                         isDescendant,
                         of_this_genot = s) & nMut(genots) == nMut(s)+1]
      if(length(t)) {
        edges <- rbind(edges,
                       data.frame(source = s,
                                  target = as.character(t),
                                  source.nmut = as.numeric(nMut(s)),
                                  target.nmut = as.numeric(nMut(s)) + 1))
      }
      
    }
    
    edges <- cbind(edge_id = paste('edge_', 1:nrow(edges), sep = ''),
                   edges)
    
    return(edges)
    
  }
  
  # plot landscape
  plotGraph <- function(df, save.plot = F) {
    
    mycolors <- c('#939598', '#d68f28', '#415ba9', '#a96cad')
    
    n_mut <- max(nMut(df$community))
    
    landscape <- df[, c('community', 'absorbance')]
    colnames(landscape) <- c('genot', 'f')
    
    edges <- makeEdges(landscape$genot)
    
    df <- cbind(edges,
                source.f = setNames(landscape$f, landscape$genot)[edges$source],
                target.f = setNames(landscape$f, landscape$genot)[edges$target])
    df$source.f[is.na(df$source.f)] <- landscape$f[landscape$genot == '']
    
    if ('color' %in% colnames(landscape)) {
      df <- merge(df, landscape[, c('genot', 'color')], by.x = 'target', by.y = 'genot')
    } else {
      df$color <- 'A'
    }
    df <- df[, c('edge_id', 'source', 'target', 'source.nmut', 'target.nmut', 'source.f', 'target.f', 'color')]
    
    dfx <- gather(df[, c(1, 4, 5)], position, nmut, source.nmut:target.nmut)
    dfx$position <- setNames(c('source', 'target'), c('source.nmut', 'target.nmut'))[dfx$position]
    
    dfy <- gather(df[, c(1, 6, 7)], position, f, source.f:target.f)
    dfy$position <- setNames(c('source', 'target'), c('source.f', 'target.f'))[dfy$position]
    
    dfxy <- merge(dfx, dfy, by = c('edge_id', 'position'))
    
    df <- merge(dfxy, df[, c('edge_id', 'color')], by = 'edge_id')
    
    dy <- min(c(max(landscape$f) - landscape$f[1], landscape$f[1] - min(landscape$f)))
    dy <- round(dy/0.1)*0.1
    ybreaks <- seq(landscape$f[1] - 10*dy, landscape$f[1] + 10*dy, by = dy)
    
    myplot <-
      ggplot(df, aes(x = nmut, y = f, group = edge_id, color = color)) +
      # geom_abline(slope = 0,
      #             intercept = landscape$f[landscape$genot == paste(rep(0, n_mut), collapse = '')],
      #             color = '#d1d3d4') +
      geom_line() +
      scale_x_continuous(name = '# of species',
                         breaks = 0:n_mut,
                         labels = as.character(0:n_mut)) +
      scale_y_continuous(name = 'Function [a.u.]',
                         breaks = pretty_breaks(n = 3),
                         expand = c(0.05, 0.05)) +
      scale_color_manual(values = setNames(mycolors, LETTERS[1:length(mycolors)])) +
      theme_bw() +
      theme(aspect.ratio = 0.6,
            panel.grid = element_blank(),
            panel.border = element_blank(),
            legend.position = 'none',
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 13)) +
      annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
      annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)
    
    if (save.plot != F) {
      ggsave(myplot,
             file = paste('../plots/', save.plot, '.pdf', sep = ''),
             dpi = 600,
             width = 100,
             height = 80,
             units = 'mm')
    }
    
    return(myplot)
    
  }
  
  return(plotGraph(df))
  
}

myplot <- plotFitnessGraph(df[, c('community', 'absorbance')])
myplot <- myplot +
  scale_y_continuous(name = expression(paste(italic(Abs)[600], ' (A.U.)', sep = '')))
# myplot <- myplot +
#   geom_smooth(data = df[df$nspecies > 0, ],
#                                aes(x = nspecies, y = absorbance, group = NA, color = NA),
#                                method = 'lm',
#                                formula = y~x,
#                                se = F,
#                                color = 'black')
focal_edges <-c('11001000', '11001001', '11001010', '11001011')
myplot <- myplot +
  geom_line(data = df[df$community %in% focal_edges[c(1, 2, 4)], ],
            aes(x = nspecies, y = absorbance, group = NA, color = NA),
            color = 'red') +
  geom_line(data = df[df$community %in% focal_edges[c(1, 3, 4)], ],
            aes(x = nspecies, y = absorbance, group = NA, color = NA),
            color = 'red')

print(myplot)

ggsave(myplot,
       filename = '../plots/pseudo_landscape.pdf',
       width = 100,
       height = 100,
       units = 'mm',
       limitsize = F)



# fitness graph for best community (00001101)
focal_comms <- c('00000000',
                 '00000001',
                 '00000100',
                 '00001000',
                 '00000101',
                 '00001001',
                 '00001100',
                 '00001101')
myplot <- plotFitnessGraph(df[df$community %in% focal_comms, c('community', 'absorbance')])
myplot <- myplot +
  scale_y_continuous(name = expression(paste(italic(Abs)[600], ' (A.U.)', sep = ''))) +
  theme(aspect.ratio = 1)

print(myplot)

ggsave(myplot,
       filename = '../plots/pseudo_landscape_3.pdf',
       width = 150,
       height = 90,
       units = 'mm',
       limitsize = F)




### DIVERSITY-FUNCTION

ggplot(do.call(data.frame,
               aggregate(absorbance ~ nspecies,
                         data = df[df$nspecies > 0, ],
                         FUN = function(x) c(mean = mean(x), sd = sd(x)))),
       aes(x = nspecies, y = absorbance.mean,
           ymin = absorbance.mean - absorbance.sd, ymax = absorbance.mean + absorbance.sd)) +
  geom_smooth(method = 'lm',
              formula = y ~ x,
              color = 'gray',
              fullrange = T) +
  geom_errorbar(width = 0) +
  geom_point() +
  scale_x_continuous(limits = c(0.2, 8.8),
                     breaks = 1:8,
                     name = '# of species') +
  scale_y_continuous(name = expression(paste(italic(Abs)[600], ' (A.U.)', sep = ''))) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none',
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 13)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf,size=0.5)

ggsave(filename = '../plots/pseudo_divfun.pdf',
       width = 100,
       height = 100,
       units = 'mm',
       limitsize = F)








### GLOBAL EPISTASIS PATTERNS

data <- data[data$wavelength == focal_wavelength, c('community', 'absorbance')]
df <- do.call(rbind,
              lapply(1:nrow(data),
                     FUN = function(i) as.data.frame(t(strsplit(as.character(data$community[i]), split = '')[[1]]))))
df$y <- data$absorbance
colnames(df) <- c(paste('sp', rev(1:landscape_size), sep = ''), 'y')

gedf <- do.call(rbind,
                lapply(1:8,
                       FUN = function(i) {
                         gedf_i <- merge(df[df[, i] == 0, -i],
                                         df[df[, i] == 1, -i],
                                         by = colnames(df)[!(colnames(df) %in% c(colnames(df)[i], 'y'))],
                                         suffixes = c('_bg', '_knockin'))[, c('y_bg', 'y_knockin')]
                         return(data.frame(species = 1 + (8-i),
                                           background_f = gedf_i$y_bg,
                                           knockin_f = gedf_i$y_knockin,
                                           delta_f = gedf_i$y_knockin - gedf_i$y_bg))
                       }))

tls_regressions <- do.call(rbind,
                           lapply(1:8,
                                  FUN = function(i) {
                                    
                                    mytls <- prcomp(gedf[gedf$species == i, c('background_f', 'knockin_f')])$rotation
                                    slope <- mytls[2, 1]/mytls[1, 1]
                                    intercept <- mean(gedf$knockin_f[gedf$species == i]) - slope*mean(gedf$background_f[gedf$species == i])
                                    
                                    return(data.frame(species = i,
                                                      slope = slope,
                                                      intercept = intercept))
                                    
                                  }))

gedf$species <- paste('Species', gedf$species)
tls_regressions$species <- paste('Species', tls_regressions$species)

# F-vs-F plot
ggplot(gedf,
       aes(x = background_f, y = knockin_f)) +
  geom_abline(slope = 1, intercept = 0, color = 'gray') +
  geom_point() +
  geom_blank(aes(y = background_f, x = knockin_f)) +
  geom_abline(data = tls_regressions,
              aes(slope = slope, intercept = intercept),
              color = 'red') +
  facet_wrap(~ species, nrow = 2) +
  scale_x_continuous(name = 'Absorbance without focal sp. (A.U.)') +
  scale_y_continuous(name = 'Absorbance with focal sp. (A.U.)') +
  theme_bw() +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14))

# dF-vs-F plot
ggplot(gedf,
       aes(x = background_f, y = delta_f)) +
  geom_abline(slope = 0, intercept = 0, color = 'gray') +
  geom_point(alpha = 0.25) +
  geom_smooth(method = 'lm',
              se = F,
              formula = y ~ x,
              color = 'black') +
  facet_wrap(~ species, nrow = 2) +
  scale_x_continuous(name = expression(Background~italic(Abs)[600]~(A.U.))) +
  scale_y_continuous(name = expression(Delta*italic(Abs)[600]~(A.U.)),
                     breaks = c(0, 1),
                     labels = c('0.0', '1.0')) +
  theme_bw() +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14),
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/fees.pdf',
       width = 140,
       height = 70,
       units = 'mm',
       limitsize = F)

# dF-vs-F plot (long)
ggplot(gedf,
       aes(x = background_f, y = delta_f)) +
  geom_abline(slope = 0, intercept = 0, color = 'gray') +
  geom_point(alpha = 0.25) +
  geom_smooth(method = 'lm',
              se = F,
              formula = y ~ x,
              color = 'black') +
  facet_wrap(~ species, nrow = 1) +
  scale_x_continuous(name = expression(Background~italic(Abs)[600]~(A.U.))) +
  scale_y_continuous(name = expression(Delta*italic(Abs)[600]~(A.U.)),
                     breaks = c(0, 1),
                     labels = c('0.0', '1.0')) +
  theme_bw() +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14),
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/fees_long.pdf',
       width = 290,
       height = 70,
       units = 'mm',
       limitsize = F)



# RANDOMIZED LANDSCAPE
n_rand <- 100
rand_linmod <- do.call(rbind,
                     lapply(1:n_rand,
                            FUN = function(r) {
                              
                              df_rand <- cbind(df[, 1:8], y = sample(df$y))
                              
                              gedf_rand <- do.call(rbind,
                                                   lapply(1:8,
                                                          FUN = function(i) {
                                                            gedf_i <- merge(df_rand[df_rand[, i] == 0, -i],
                                                                            df_rand[df_rand[, i] == 1, -i],
                                                                            by = colnames(df_rand)[!(colnames(df_rand) %in% c(colnames(df_rand)[i], 'y'))],
                                                                            suffixes = c('_bg', '_knockin'))[, c('y_bg', 'y_knockin')]
                                                            return(data.frame(species = 1 + (8-i),
                                                                              background_f = gedf_i$y_bg,
                                                                              knockin_f = gedf_i$y_knockin,
                                                                              delta_f = gedf_i$y_knockin - gedf_i$y_bg))
                                                          }))
                              
                              linmod_params <- do.call(rbind,
                                                       lapply(1:8,
                                                              FUN = function(sp) {
                                                                
                                                                mylm_i <- lm(delta_f ~ background_f,
                                                                             data = gedf_rand[gedf_rand$species == sp, ])
                                                                return(data.frame(species = paste('Species', sp),
                                                                                  slope = as.numeric(mylm_i$coefficients[2]),
                                                                                  intercept = as.numeric(mylm_i$coefficients[1]),
                                                                                  R2 = summary(mylm_i)$r.squared))
                                                                
                                                              }))
                              
                              return(cbind(type = 'rnd', run = r, linmod_params))
                              
                            }))

emp_linmod <- cbind(type = 'emp', run = NA,
                    do.call(rbind,
                            lapply(1:8,
                                   FUN = function(sp) {
                                     
                                     mylm_i <- lm(delta_f ~ background_f,
                                                  data = gedf[gedf$species == paste('Species', sp), ])
                                     return(data.frame(species = paste('Species', sp),
                                                       slope = as.numeric(mylm_i$coefficients[2]),
                                                       intercept = as.numeric(mylm_i$coefficients[1]),
                                                       R2 = summary(mylm_i)$r.squared))
                                     
                                   })))

plot_this <- rbind(emp_linmod, rand_linmod)

# all randomizations
ggplot(gedf,
       aes(x = background_f, y = delta_f)) +
  geom_abline(slope = 0, intercept = 0, color = 'gray') +
  geom_point(alpha = 0.25) +
  # geom_smooth(method = 'lm',
  #             se = F,
  #             formula = y ~ x,
  #             color = 'black') +
  geom_abline(data = plot_this[plot_this$type == 'rnd', ],
              aes(slope = slope, intercept = intercept, group = run),
              alpha = 0.05,
              color = 'deepskyblue') +
  geom_abline(data = plot_this[plot_this$type == 'emp', ],
              aes(slope = slope, intercept = intercept),
              alpha = 1,
              color = 'black') +
  facet_wrap(~ species, nrow = 2) +
  scale_x_continuous(name = expression(Background~italic(Abs)[600]~(A.U.))) +
  scale_y_continuous(name = expression(Delta*italic(Abs)[600]~(A.U.)),
                     breaks = c(0, 1),
                     labels = c('0.0', '1.0')) +
  theme_bw() +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 14),
        panel.border = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/fees_vs_randomLanscape.pdf',
       width = 140,
       height = 70,
       units = 'mm',
       limitsize = F)

# one randomization
df_rand <- cbind(df[, 1:8], y = sample(df$y))
gedf_rand <- do.call(rbind,
                     lapply(1:8,
                            FUN = function(i) {
                              gedf_i <- merge(df_rand[df_rand[, i] == 0, -i],
                                              df_rand[df_rand[, i] == 1, -i],
                                              by = colnames(df_rand)[!(colnames(df_rand) %in% c(colnames(df_rand)[i], 'y'))],
                                              suffixes = c('_bg', '_knockin'))[, c('y_bg', 'y_knockin')]
                              return(data.frame(species = paste('Species', 1 + (8-i)),
                                                background_f = gedf_i$y_bg,
                                                knockin_f = gedf_i$y_knockin,
                                                delta_f = gedf_i$y_knockin - gedf_i$y_bg))
                            }))
ggplot(gedf_rand,
       aes(x = background_f, y = delta_f)) +
  geom_abline(slope = 0, intercept = 0, color = 'gray') +
  geom_point(alpha = 0.25) +
  geom_smooth(method = 'lm',
              se = F,
              formula = y ~ x,
              color = 'deepskyblue',
              fullrange = T) +
  facet_wrap(~ species, nrow = 2) +
  scale_x_continuous(name = expression(Background~italic(Abs)[600]~(A.U.))) +
  scale_y_continuous(name = expression(Delta*italic(Abs)[600]~(A.U.)),
                     breaks = c(0, 1),
                     labels = c('0.0', '1.0')) +
  theme_bw() +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        panel.border = element_blank(),
        axis.ticks = element_blank()) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/fees_sampleRandomLanscape.pdf',
       width = 100,
       height = 50,
       units = 'mm',
       limitsize = F)





### DISTRIBUTIONS OF FUNCTIONAL INTERACTIONS

# pairwise interactions
pairs <- t(combn(8, 2))
dfi <- do.call(rbind,
               lapply(1:nrow(pairs),
                      FUN = function(pair_i) {
                        
                        p <- pairs[pair_i, ]
                        p_cols <- 1 + (8 - p)
                        
                        bg <- df[df[, p_cols[1]] == 0 & df[, p_cols[2]] == 0, ]
                        bg_i <- df[df[, p_cols[1]] == 1 & df[, p_cols[2]] == 0, ]
                        bg_j <- df[df[, p_cols[1]] == 0 & df[, p_cols[2]] == 1, ]
                        bg_ij <- df[df[, p_cols[1]] == 1 & df[, p_cols[2]] == 1, ]
                        
                        out <- merge(merge(merge(bg[, -p_cols], bg_i[, -p_cols],
                                                 by = colnames(df[, 1:8])[-p_cols],
                                                 suffixes = c('', '_i')),
                                           bg_j[, -p_cols],
                                           by = colnames(df[, 1:8])[-p_cols],
                                           suffixes = c('', '_j')),
                                     bg_ij[, -p_cols],
                                     by = colnames(df[, 1:8])[-p_cols],
                                     suffixes = c('', '_ij'))
                        out$delta_i <- out$y_i - out$y
                        out$delta_j <- out$y_j - out$y
                        out$eps_ij <- out$y_ij - out$y_i - out$y_j + out$y
                        
                        out <- cbind(out, pi = 0, pj = 0)
                        colnames(out)[(ncol(out) - 1):ncol(out)] <- paste('sp', p, sep = '')
                        bg <- sapply(1:nrow(out),
                                     FUN = function(i) paste(out[i, paste('sp', 8:1, sep = '')], collapse = ''))
                        
                        out <- data.frame(species_i = paste('Species', p[1]),
                                          species_j = paste('Species', p[2]),
                                          background = bg,
                                          delta_i = out$delta_i,
                                          delta_j = out$delta_j,
                                          eps_ij = out$eps_ij)
                        
                        return(out)
                        
                      }))

ggplot(dfi, aes(x = 0, y = eps_ij)) +
  geom_hline(yintercept = 0,
             color = 'gray') +
  geom_violin(alpha = 0.25,
              color = NA,
              fill = 'black') +
  geom_point(cex = 0.5) +
  facet_grid2(species_i ~ species_j,
              render_empty = F,
              switch = 'y') +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(name = expression(paste('Functional interaction between species ', italic (i), ' and ', italic(j), ', ', italic(epsilon)[italic(ij)], sep ='')),
                     position = 'right',
                     breaks = pretty_breaks(n = 2)) +
  theme_bw() +
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13,
                                  angle = 0),
        strip.text.y.left = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

ggsave(filename = '../plots/dfi.pdf',
       width = 150,
       height = 150,
       units = 'mm',
       limitsize = F)

# functional effects
dfe <- do.call(rbind,
               lapply(1:8,
                      FUN = function(sp_i) {
                        
                        sp_col <- 1 + (8 - sp_i)
                        
                        bg <- df[df[, sp_col] == 0, ]
                        bg_i <- df[df[, sp_col] == 1, ]
                        
                        out <- merge(bg[, -sp_col], bg_i[, -sp_col],
                                     by = colnames(df[, 1:8])[-sp_col],
                                     suffixes = c('', '_i'))
                        out$delta_i <- out$y_i - out$y
                        out <- cbind(out, pi = 0)
                        colnames(out)[ncol(out)] <- paste('sp', sp_i, sep = '')
                        bg <- sapply(1:nrow(out),
                                     FUN = function(i) paste(out[i, paste('sp', 8:1, sep = '')], collapse = ''))
                        
                        out <- data.frame(species = paste('Species', sp_i),
                                          background = bg,
                                          delta_i = out$delta_i)
                        
                        return(out)
                        
                      }))

ggplot(dfe, aes(x = species, y = delta_i, fill = species)) +
  geom_hline(yintercept = 0,
             color = 'gray') +
  geom_violin(alpha = 0.25,
              color = NA,
              width = 1) +
  geom_point(cex = 0.5) +
  scale_y_continuous(name = expression(paste('Functional effect, ', Delta*italic(Abs), ' (A.U.)', sep = '')),
                     breaks = pretty_breaks(n = 3)) +
  scale_fill_manual(values = c('#E76D57',
                               '#ACD152',
                               '#7A9DD2',
                               '#F7AFBE',
                               '#E5C91F',
                               '#A88BC0',
                               '#5CA985',
                               '#BBAD75')) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 13,
                                  angle = 0),
        strip.text.y.left = element_text(angle = 0),
        strip.text.x = element_text(angle = 90),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        legend.position = 'none') +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/dfe.pdf',
       width = 250,
       height = 100,
       units = 'mm',
       limitsize = F)


if(FALSE) {
  
  
  ### BEYOND-PAIRWISE DOMINANCE MODEL
  
  singles <- data[data$nspecies == 1, ]
  
  tst <- sapply(data$community[data$nspecies > 1],
                FUN = function(comm) {
                  
                  singles_in_comm <- sapply(which(strsplit(comm, split = '')[[1]] == 1),
                                            FUN = function(sp) {
                                              out <- rep(0, 8)
                                              out[sp] <- 1
                                              out <- paste(out, collapse = '')
                                              return(out)
                                            })
                  singles_in_comm <- singles[singles$community %in% singles_in_comm, ]
                  
                  comm_spect <- data[data$community == comm, c('wavelength', 'absorbance', 'sd')]
                  comm_spect <- merge(singles_in_comm, comm_spect,
                                      by = 'wavelength',
                                      all = T,
                                      suffixes = c('_singles', '_coculture'))
                  comm_spect$diff <- comm_spect$absorbance_coculture - comm_spect$absorbance_singles
                  comm_spect <- aggregate(diff ~ community,
                                          data = comm_spect,
                                          FUN = function(x) sum(x^2))
                  dom_monoc <- comm_spect$community[which.min(comm_spect$diff)]
                  
                  return(dom_monoc)
                  
                })
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  ### COMPLETE EXCLUSION MODEL
  # species are ranked from most to least competitive
  # the spectrum of a consortium c is modeled as the spectrum of the most competitive species in c
  
  data <- data_full
  monoc <- data[data$nspecies == 1 & data$wavelength == focal_wavelength, ]
  consortia <- data[data$nspecies > 1 & data$wavelength == focal_wavelength, ]
  unique_consortia <- unique(consortia$community)
  
  set.seed(0)
  ranks <- sample(1:8)
  
  exclusionModel <- function(target_consortia, ranks = 1:8) {
    
    do.call(rbind,
            lapply(target_consortia,
                   FUN = function(target_consortium) {
                     
                     sp_in_target <- which(strsplit(target_consortium, split = '')[[1]] == '1')
                     most_competitive_sp <- sp_in_target[which.min(ranks[sp_in_target])]
                     
                     pred <- rep(0, 8)
                     pred[most_competitive_sp] <- 1
                     pred <- paste(pred, collapse = '')
                     
                     predicted_spect <- monoc[monoc$community == pred, c('community', 'wavelength', 'absorbance', 'sd')]
                     colnames(predicted_spect)[3:4] <- paste(colnames(predicted_spect)[3:4], '_predicted', sep = '')
                     predicted_spect$community <- target_consortium
                     
                     return(predicted_spect)
                     
                   }))
    
  }
  
  # find optimal rank
  ranks_to_test <- permn(8)
  st <- Sys.time()
  ranks_R2 <- unlist(pblapply(ranks_to_test,
                              FUN = function(rank_i) {
                                
                                predictedSpect <- exclusionModel(unique_consortia, ranks = rank_i)
                                predictedSpect <- merge(predictedSpect, consortia, by = c('community', 'wavelength'), all = T)
                                
                                R2 <- summary(lm(absorbance ~ absorbance_predicted,
                                                 data = predictedSpect))$r.squared
                                
                                return(R2)
                                
                              }))
  end <- Sys.time()
  end - st
  
  
  
  
  
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
  
  # adjust colors for nice aesthetics
  mycolors <- c('#E76D57', '#ACD152', '#7A9DD2', '#F7AFBE', '#E5C91F', '#A88BC0', '#5CA985', '#BBAD75')
  
  myplot <- myplot +
    scale_color_manual(name = 'Species',
                       values = c(mycolors, 'black', 'black')) +
    scale_fill_manual(name = 'Species',
                      values = c(mycolors, 'black', 'black'))
  
  ggsave(myplot,
         filename = paste('../plots/pseudo_spectra.pdf', sep = ''),
         width = 1500,
         height = 200,
         units = 'mm',
         limitsize = F)
  
  
  myplot <- myplot +
    scale_y_log10(breaks = c(0.1, 1, 10),
                  name = 'Absorbance')
  
  ggsave(myplot,
         filename = paste('../plots/pseudo_spectra_log.pdf', sep = ''),
         width = 1500,
         height = 200,
         units = 'mm',
         limitsize = F)
  
  
  # small plot for main fig
  data <- data[order(data$community_id), ]
  data$community <- factor(data$community,
                           levels = unique(data$community))
  myplot_main <- 
    ggplot(data,
           aes(x = wavelength, y = absorbance, ymin = absorbance - sd, ymax = absorbance + sd)) +
    geom_ribbon(color = NA,
                fill = 'black',
                alpha = 0.2) +
    geom_line(color = 'black') +
    # facet_wrap(~ community,
    #            nrow = 8,
    #            dir = 'v') +
    facet_wrap2(~ community,
                nrow = 8,
                dir = 'v',
                axes = 'all',
                remove_labels = 'all') +
    # scale_y_continuous(name = 'Absorbance (A.U.)',
    #                    breaks = pretty_breaks(n = 2)) +
    scale_y_log10(name = 'Absorbance (A.U.)',
                  breaks = seq(0.1, 1, length.out = 10)) +
    scale_x_continuous(name = 'Wavelength (nm)',
                       breaks = c(500, 700)) +
    theme_bw() + 
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = 'none',
          panel.border = element_blank(),
          axis.line.x = element_line(),
          axis.line.y = element_line(),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          #strip.text = element_text(size = 8),
          axis.ticks.length = unit(0.25, "mm"),
          #axis.ticks = element_blank(),
          strip.text = element_blank())
  # annotate("segment", x=-Inf, xend=Inf, y=0.09, yend=0.09, linewidth=0.5) +
  # annotate("segment", x=-Inf, xend=-Inf, y=0.09, yend=Inf, linewidth=0.5)
  
  ggsave(myplot_main,
         filename = paste('../plots/pseudo_spectra_main_log_ticks.pdf', sep = ''),
         width = 300,
         height = 150,
         units = 'mm',
         limitsize = F)
  
  
  ### QUANTIFY INTERACTIONS
  
  # get interaction coefficients at all orders
  focal_wavelength <- 600
  data <- data[data$wavelength == focal_wavelength, c('community', 'absorbance')]
  df <- do.call(rbind,
                lapply(1:nrow(data),
                       FUN = function(i) as.data.frame(t(strsplit(as.character(data$community[i]), split = '')[[1]]))))
  df$y <- data$absorbance
  colnames(df) <- c(paste('x', 1:landscape_size), 'y')
  
  # which base?
  base <- 'taylor' # one of 'taylor' or 'fourier'
  if (base == 'fourier') {
    for (i in 1:landscape_size) df[df[, i] == 0, i] <- -1
  }
  
  # model <- lm(y ~ . + .*. + .*.*. + .*.*.*. + .*.*.*.*. + .*.*.*.*.*. + .*.*.*.*.*.*. + .*.*.*.*.*.*.*.,
  #             data = df)
  model <- lm(y ~ .*.*.*.*.*.*.*.,
              data = df)
  coef <- model$coefficients
  names(coef) <- gsub('`1', '', names(coef))
  names(coef) <- gsub('`x', '', names(coef))
  names(coef) <- gsub(' ', '', names(coef))
  names(coef) <- gsub(':', '', names(coef))
  names(coef)[1] <- ''
  
  coef_df <- data.frame(coef_id = names(coef),
                        coef_order = nchar(names(coef)),
                        value = as.numeric(coef))
  ggplot(coef_df, aes(x = coef_order, y = value, group = coef_order)) +
    geom_abline(slope = 0,
                intercept = 0,
                color = 'gray') +
    geom_boxplot(data = coef_df[coef_df$coef_order > 0 & coef_df$coef_order < 8, ],
                 fill = 'white',
                 outlier.shape = NA) +
    geom_jitter(width = 0.15,
                alpha = 0.25) +
    scale_x_continuous(name = 'Order',
                       breaks = 0:8) +
    scale_y_continuous(name = expression(paste('Coefficient value, ', italic(beta), ' (A.U.)', sep = ''))) +
    theme_bw() + 
    theme(aspect.ratio = 0.6,
          panel.grid = element_blank(),
          strip.background = element_blank(),
          legend.position = 'none',
          panel.border = element_blank(),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 16)) +
    annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
    annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)
  
  ggsave(filename = paste('../plots/pseudo_coef_', base, '.pdf', sep = ''),
         width = 100,
         height = 100,
         units = 'mm',
         limitsize = F)
  
  model_y <- predict(model, df[, 1:landscape_size])
  
  
}







