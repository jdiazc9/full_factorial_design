library(ggplot2)
library(ggbreak)

# consider N = 8 species
N <- 8

### LOAD DATA
df <- read.table('../data/pseudo.txt', header = T, colClasses = c('character', 'character', 'numeric', 'numeric', 'numeric'))

# filter out highest dil. factors and aggregate replicates
df <- df[df$dilution_factor == 0.0025 & df$wavelength == 600, ]
df <- do.call(data.frame,
              aggregate(absorbance ~ community,
                        df,
                        FUN = mean))

communities <- do.call(rbind,
                       lapply(1:nrow(df),
                              FUN = function (i) as.numeric(strsplit(df$community[i], split = '')[[1]])))
colnames(communities) <- paste('sp', N:1, sep = '')
df <- as.data.frame(cbind(communities, y = df$absorbance))
fun <- df$y




# Taylor base
# F = a_0 + a_1*x_1 + a_2*x_2 + ... + a_12*x_1*x_2 + ... + a_123*x_1*x_2*x_3 + ...; with x_i = 0,1; {a} are the coefficients
mydf <- as.data.frame(cbind(communities, fun))
taylor_lm <- lm(fun ~ .*.*.*.*.*.*.*., # eight dots here because there are N=8 species and we want to fit the full model (if we wanted to truncate at e.g. third order this would be fun ~ .*.*.)
                data = mydf)

# Fourier base
# F = b_0 + b_1*x_1 + b_2*x_2 + ... + b_12*x_1*x_2 + ... + b_123*x_1*x_2*x_3 + ...; with x_i = -1,1; {b} are the coefficients
communities[communities == 0] <- -1
mydf <- as.data.frame(cbind(communities, fun))
fourier_lm <- lm(fun ~ .*.*.*.*.*.*.*.,
                 data = mydf)

# make data frame to plot
plot_this <- data.frame(coef_id = c(names(taylor_lm$coefficients), names(fourier_lm$coefficients)),
                        base = rep(c('taylor', 'fourier'), each = 2^N),
                        coef_value = c(taylor_lm$coefficients, fourier_lm$coefficients))
plot_this$order <- sapply(plot_this$coef_id,
                          function(x) 1 + sum(strsplit(x, split = '')[[1]] == ':'))
plot_this$order[grepl('\\(', plot_this$coef_id)] <- 0
plot_this$base <- factor(plot_this$base,
                         levels = c('taylor', 'fourier'))

# fourier plot
plot_this$coef_value[plot_this$base == 'fourier' & plot_this$coef_id == '(Intercept)'] <- plot_this$coef_value[plot_this$base == 'fourier' & plot_this$coef_id == '(Intercept)'] - 0.8
ggplot(plot_this[plot_this$base == 'fourier', ], aes(x = order, y = coef_value, group = order)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = 'gray') +
  geom_boxplot(data = plot_this[plot_this$base == 'fourier' & plot_this$order > 0 & plot_this$order < 8, ],
               fill = 'white',
               outlier.shape = NA) +
  geom_jitter(width = 0.15,
              alpha = 0.25) +
  scale_x_continuous(name = 'Order',
                     breaks = 0:8,
                     limits = c(-0.4, 8.4),
                     expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste('Coefficient value, ', italic(beta), ' (A.U.)', sep = '')),
                     limits = c(-0.1, 0.2),
                     breaks = seq(-0.1, 0.2, by = 0.1),
                     labels = c('-0.1', '0.0', '0.1', '1.0')) +
  coord_cartesian(clip = 'off') +
  theme_bw() + 
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=0.1475, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=0.1525, yend=Inf, linewidth=0.5) +
  annotate("segment", x=-0.4, xend=-0.3, y=0.1525, yend=0.1525, linewidth=0.5) +
  annotate("segment", x=-0.4, xend=-0.3, y=0.1475, yend=0.1475, linewidth=0.5)

ggsave(filename = '../plots/pseudo_coef_fourier.pdf',
       width = 100,
       height = 100,
       units = 'mm',
       limitsize = F)

# taylor plot
ggplot(plot_this[plot_this$base == 'taylor', ], aes(x = order, y = coef_value, group = order)) +
  geom_abline(slope = 0,
              intercept = 0,
              color = 'gray') +
  geom_boxplot(data = plot_this[plot_this$base == 'taylor' & plot_this$order > 0 & plot_this$order < 8, ],
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

ggsave(filename = '../plots/pseudo_coef_taylor.pdf',
       width = 100,
       height = 100,
       units = 'mm',
       limitsize = F)

# variance explained by each order:
# the variance explained by interactions of order k is quantified as
# the sum of squares of the Fourier coefficients of order k (k >= 1)
var_explained <- sapply(1:N,
                        function(i) sum(plot_this$coef_value[plot_this$base == 'fourier' & plot_this$order == i]^2))
frac_var_explained <- var_explained / mean((fun - mean(fun))^2)

ggplot(data.frame(x = 1:N, y = frac_var_explained),
       aes(x = x, y = y)) +
  geom_line() +
  geom_line(aes(x = 1:N, y = cumsum(frac_var_explained)),
            linetype = 'dashed') +
  scale_x_continuous(name = 'Order') +
  scale_y_continuous(name = 'Fraction of\nvariance explained',
                     limits = c(0, 1.01)) +
  theme_bw() + 
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/pseudo_fracVariance.pdf',
       width = 65,
       height = 65,
       units = 'mm',
       limitsize = F)

# larger plot (for main figure)
ggplot(data.frame(x = 1:N, y = frac_var_explained),
       aes(x = x, y = y)) +
  geom_line() +
  geom_line(aes(x = 1:N, y = cumsum(frac_var_explained)),
            linetype = 'dashed') +
  scale_x_continuous(name = 'Order',
                     breaks = 1:8) +
  scale_y_continuous(name = 'Fraction of\nvariance explained',
                     limits = c(0, 1.01)) +
  theme_bw() + 
  theme(aspect.ratio = 0.7,
        panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 16)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, linewidth=0.5) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, linewidth=0.5)

ggsave(filename = '../plots/pseudo_fracVariance_large.pdf',
       width = 100,
       height = 100,
       units = 'mm',
       limitsize = F)

