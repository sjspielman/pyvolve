# Post-processing and plotting in R:
require(ggplot2)
require(cowplot)
require(dplyr)
dat <- read.csv("compare_scaling_results.csv")
dat %>% group_by(scaling, omega) %>% summarize(mean_tl = mean(treelen), sd_tl = sd(treelen)) -> dat.sum
limits <- aes(ymax = mean_tl + sd_tl, ymin = mean_tl - sd_tl)

theme_set(theme_cowplot() + theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11.5)))
ggplot(dat.sum, aes(x = omega, y = mean_tl, shape = scaling, color = scaling)) + 
geom_point(size = 3) + geom_line(size = 0.45, alpha = 0.75) + 
geom_errorbar(limits, width = 0.02, alpha = 0.75) + 
scale_shape_manual(values = c(18, 20), name = "Matrix Scaling", labels = c("Neutral", "Yang")) + 
scale_color_manual(values = c("dodgerblue3", "darkgoldenrod2"), name = "Matrix Scaling", labels = c("Neutral", "Yang")) + 
xlab("Simulated dN/dS") + ylab("Mean Tree Length") -> p

ggsave("yang_neutral_scaling.pdf", p, width = 5.25, height = 3.75)
    
    
    
    
    
    