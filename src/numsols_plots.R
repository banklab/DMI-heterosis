###############################################################################
# numsols_plots.r
# Julio Ayala Lopez
# 08.2024
#
# Script to generate heatmap plots of allele frequency and linkage 
# disequilibrium using numerical solutions of a two locus DMIs/heterosis ODE. 
###############################################################################


library(ggplot2)
library(dplyr)
library(rstudioapi)
library(viridis)

setwd(dirname(getActiveDocumentContext()$path)) # Works only if executing code manually, not sourcing



dfcodominant01 <- read.csv("../Data/codominant_sols_freqs_20250911_a=0.01.csv", header = F, sep = ",", col.names = c("ga", "r", "m","pA", "LD"))
dfcodominant001 <- read.csv("../Data/codominant_sols_freqs_20250911_a=0.001.csv", header = F, sep = ",", col.names = c("ga", "r", "m","pA", "LD"))
dfreces02 <- read.csv("../Data/recessive_sols_ext.csv", header = F, sep = ",", col.names = c("ga", "r", "m","pA", "LD"))
dfreces <- read.csv("../Data/recessive_sols_ext_001.csv", header = F, sep = ",", col.names = c("ga", "r", "m","pA", "LD"))


dfcodominant01$r2 <- (dfcodominant01$LD**2)/((dfcodominant01$pA**2)*((1-dfcodominant01$pA)**2))
dfcodominant01$rs <- dfcodominant01$r/0.01
dfcodominant001$r2 <- (dfcodominant001$LD**2)/((dfcodominant001$pA**2)*((1-dfcodominant001$pA)**2))
dfcodominant001$rs <- dfcodominant001$r/0.001
dfreces$r2 <- (dfreces$LD**2)/((dfreces$pA**2)*((1-dfreces$pA)**2))
dfreces$rs <- dfreces$r/0.001


dfreces <- dfreces %>% filter(r>0)

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# dfcodominant01 %>% filter(r<=0.49 & m==0 & ga<=8) %>%
#   ggplot(aes(x = ga, y=r/0.01, fill=pA)) + 
#   geom_tile() +
#   #scale_fill_gradientn("Freq. A", colors = viridis(n=256,end = 0.7,direction =-1),
#   scale_fill_gradientn(expression(P[A]), colors = tail(hcl.colors(100, "ag_GrnYl", rev=T),95),
#                        limits = c(0.24,0.51),
#                        breaks = c(0.25,0.375,0.5)) +
#   geom_contour(aes(z= pA), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
#   scale_y_continuous(trans='log10') +
#   #facet_grid(vars(m)) +
#   labs(title="Allele frequency  (a=0.01)") + 
#   xlab(expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")"))) + 
#   #xlab(expression(paste("line1 \nline2 a"))) +
#   ylab("Recombination rate") +
#   labs(fill=expression(P[A])) +
#   #guides(y.sec = guide_none("migration rate")) +
#   theme_bw() +
#   theme(text=element_text(size=18), 
#         axis.text=element_text(size=16), 
#         axis.title=element_text(size=18), 
#         plot.title=element_text(size=18), 
#         legend.text=element_text(size=16), 
#         legend.title=element_text(size=16, vjust=0.85),
#         legend.position = "bottom") +
#   xlim(0,4.25) +
#   guides(fill = guide_colourbar(barwidth = 31, barheight = 1))
# 
# dfcodominant001 %>% filter(r<=0.49 & m==0 & ga<=8) %>%
#   ggplot(aes(x = ga, y=r, fill=pA)) + 
#   geom_tile() +
#   #scale_fill_gradientn("Freq. A", colors = viridis(n=256,end = 0.7,direction =-1),
#   scale_fill_gradientn(expression(P[A]), colors = tail(hcl.colors(100, "ag_GrnYl", rev=T),95),
#                        limits = c(0.24,0.51),
#                        breaks = c(0.25,0.375,0.5)) +
#   geom_contour(aes(z= pA), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
#   scale_y_continuous(trans='log10') +
#   #facet_grid(vars(m)) +
#   labs(title="Allele frequency (a=0.001)") + 
#   xlab(expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")"))) + 
#   #xlab(expression(paste("line1 \nline2 a"))) +
#   ylab("Recombination rate") +
#   labs(fill=expression(P[A])) +
#   #guides(y.sec = guide_none("migration rate")) +
#   theme_bw() +
#   theme(text=element_text(size=18), 
#         axis.text=element_text(size=16), 
#         axis.title=element_text(size=18), 
#         plot.title=element_text(size=18), 
#         legend.text=element_text(size=16), 
#         legend.title=element_text(size=16, vjust=0.85),
#         legend.position = "bottom") +
#   xlim(0,4.25) +
#   guides(fill = guide_colourbar(barwidth = 31, barheight = 1))


## paper
dfcodominant001 %>%
  filter(r<=0.49 & m==0 & ga<=8 & r>=0.00000001 & rs<=100 & rs>=0.0001) %>%
  ggplot(aes(x = ga, y=rs, fill=pA)) + 
  geom_tile() +
  scale_fill_gradientn(expression(P[A]),
                       colors = tail(hcl.colors(100, "ag_GrnYl", rev=T),95),
                       limits = c(0.24,0.51),
                       breaks = c(0.25,0.375,0.5)) +
  geom_contour(aes(z= pA), colour="#FFF", binwidth = 0.02,
               linetype = "dashed", alpha=0.6) +
  # y-axis with automatic powers of 10
  scale_y_continuous(
    trans = 'log10',
    breaks = function(lims) {
      rng <- range(lims, na.rm = TRUE)
      pmin <- floor(log10(rng[1]))
      pmax <- ceiling(log10(rng[2]))
      10^seq(pmin, pmax)
    },
    labels = function(b) {
      exps <- log10(b)
      as.expression(lapply(exps, function(e) bquote(10^.(e))))
    }
  ) +
  # x-axis with float labels
  scale_x_continuous(
    limits = c(0, 4.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    #title=expression("Allele frequency  ("*alpha*"=0.001)"),
    x=expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")")),
    y=expression(r/alpha),
    fill=expression(P[A])
  ) +
  theme_linedraw() +
  theme(
    text=element_text(family = "Latin Modern Math",size=20), 
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20), 
    plot.title=element_text(size=20), 
    legend.text=element_text(size=18), 
    legend.title=element_text(size=20, vjust=0.85),
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 31, barheight = 1))

ggsave("../Fig/Fig_1A.png", width = 23, height = 20, units = "cm")

###############


## Updated with larger font
# dfcodominant01 %>% filter(r<=0.49 & m==0 & ga<=8 & r>=0.000001) %>%
#   ggplot(aes(x = ga, y=rs, fill=pA)) + 
#   geom_tile() +
#   #scale_fill_gradientn("Freq. A", colors = viridis(n=256,end = 0.7,direction =-1),
#   scale_fill_gradientn(expression(P[A]), colors = tail(hcl.colors(100, "ag_GrnYl", rev=T),95),
#                        limits = c(0.24,0.51),
#                        breaks = c(0.25,0.375,0.5)) +
#   geom_contour(aes(z= pA), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
#   scale_y_continuous(trans='log10') +
#   #facet_grid(vars(m)) +
#   labs(title="Allele frequency") + 
#   xlab(expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")"))) + 
#   #xlab(expression(paste("line1 \nline2 a"))) +
#   ylab("Recombination rate") +
#   labs(fill=expression(P[A])) +
#   #guides(y.sec = guide_none("migration rate")) +
#   theme_bw() +
#   theme(text=element_text(size=24), 
#         axis.text=element_text(size=20), 
#         axis.title=element_text(size=24), 
#         plot.title=element_text(size=24), 
#         legend.text=element_text(size=20), 
#         legend.title=element_text(size=18, vjust=0.85),
#         legend.position = "bottom") +
#   xlim(0,4.25) +
#   guides(fill = guide_colourbar(barwidth = 30, barheight = 1))
# ####
# 
# ggsave("../Fig/freq_heatmap_20x22.pdf", width = 20, height = 22, units = "cm")
# ggsave("../Fig/freq_heatmap_20x22_PPT.pdf", width = 20, height = 22, units = "cm")


# dfcodominant01 %>% filter(r<=0.49 & m==0 & ga<=8 & r>=0.000001) %>%
#   ggplot(aes(x = ga, y=r, fill=LD)) + 
#   geom_tile() +
#   scale_fill_gradientn(colors = rocket(n=10,end = 0.7,direction =1)) +
#   geom_contour(aes(z= LD), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
#   #scale_fill_gradientn(colors = hcl.colors(30, "BuPu", rev=T)) +
#   scale_y_continuous(trans='log10') +
#   #facet_grid(vars(m)) +
#   labs(title="Linkage disequilibrium") + 
#   xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
#   #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
#   ylab("Recombination rate") +
#   #guides(y.sec = guide_none("migration rate")) +
#   theme_bw() +
#   theme(text=element_text(size=18), 
#         axis.text=element_text(size=16), 
#         axis.title=element_text(size=18), 
#         plot.title=element_text(size=18), 
#         legend.text=element_text(size=16), 
#         legend.title=element_text(size=16, vjust=0.85),
#         legend.position = "bottom") +
#   xlim(0,4.25) +
#   guides(fill = guide_colourbar(barwidth = 30, barheight = 1))
# 
# 
# dfcodominant001 %>% filter(r<=0.49 & m==0 & ga<=8 & r>=0.000001) %>%
#   ggplot(aes(x = ga, y=r, fill=LD)) + 
#   geom_tile() +
#   scale_fill_gradientn(colors = rocket(n=10,end = 0.7,direction =1)) +
#   geom_contour(aes(z= LD), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
#   #scale_fill_gradientn(colors = hcl.colors(30, "BuPu", rev=T)) +
#   scale_y_continuous(trans='log10') +
#   #facet_grid(vars(m)) +
#   labs(title="Linkage disequilibrium") + 
#   xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
#   #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
#   ylab("Recombination rate") +
#   #guides(y.sec = guide_none("migration rate")) +
#   theme_bw() +
#   theme(text=element_text(size=18), 
#         axis.text=element_text(size=16), 
#         axis.title=element_text(size=18), 
#         plot.title=element_text(size=18), 
#         legend.text=element_text(size=16), 
#         legend.title=element_text(size=16, vjust=0.85),
#         legend.position = "bottom") +
#   xlim(0,4.25) +
#   guides(fill = guide_colourbar(barwidth = 30, barheight = 1))


## Paper 
dfcodominant001 %>%
  filter(r<=0.49 & m==0 & ga<=8 & r>=0.00000001 & rs<=100 & rs>=0.0001) %>%
  ggplot(aes(x = ga, y=rs, fill=LD)) + 
  geom_tile() +
  scale_fill_gradientn(colors = rocket(n=10,end = 0.7,direction =1)) +
  geom_contour(aes(z= LD), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
  # y-axis with automatic powers of 10
  scale_y_continuous(
    trans = 'log10',
    breaks = function(lims) {
      rng <- range(lims, na.rm = TRUE)
      pmin <- floor(log10(rng[1]))
      pmax <- ceiling(log10(rng[2]))
      10^seq(pmin, pmax)
    },
    labels = function(b) {
      exps <- log10(b)
      as.expression(lapply(exps, function(e) bquote(10^.(e))))
    }
  ) +
  # x-axis with float labels
  scale_x_continuous(
    limits = c(0, 4.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    #title=expression("Linkage disequilibrium  ("*alpha*"=0.001)"),
    x=expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")")),
    y=expression(r/alpha),
    fill=expression(D)
  ) +
  theme_linedraw() +
  theme(
    text=element_text(family = "Latin Modern Math",size=20), 
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20), 
    plot.title=element_text(size=20), 
    legend.text=element_text(size=18), 
    legend.title=element_text(size=20, vjust=0.85),
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 31, barheight = 1))

ggsave("../Fig/Fig_1B.png", width = 23, height = 20, units = "cm")


###########


## Plotting r2
dfcodominant001 %>% filter(r<=0.49 & m==0 & ga<=8 & r>=0.0000001) %>%
  ggplot(aes(x = ga, y=r, fill=r2)) + 
  geom_tile() +
  #scale_fill_gradientn(colors = hcl.colors(30, "BuPu", rev=T),
  scale_fill_gradientn(colors = rocket(n=10,end = 0.65,direction =-1),
                       breaks = c(0,0,0.5,1), 
                       limits = c(-0.05,1.05)) +
  scale_y_continuous(trans='log10') +
  #facet_grid(vars(m)) +
  labs(title=expression("Linkage disequilibrium ("*r^2*")"), fill=expression(r^2)) + 
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
  #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
  ylab("Recombination rate") +
  #guides(y.sec = guide_none("migration rate") +
  theme_bw() +
  theme(text=element_text(size=18), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        plot.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16))


ggsave("../Fig/LD_r2.pdf", width = 20, height = 22, units = "cm")


###### Plots with recessive DMIs

dfreces %>% filter(r<=0.49 & m==0 & r>=0.0000001 & rs>0.001 & rs<=50) %>%
  ggplot(aes(x = ga, y=rs, fill=pA)) + 
  geom_tile() +
  scale_fill_gradientn(expression(P[A]),
                       colors = tail(hcl.colors(100, "ag_GrnYl", rev=T),95),
                       limits = c(0.18,0.51),
                       breaks = c(0.2,0.3,0.4,0.5)) +
  geom_contour(aes(z= pA), colour="#FFF", binwidth = 0.02,
               linetype = "dashed", alpha=0.6) +
  # y-axis with automatic powers of 10
  scale_y_continuous(
    trans = 'log10',
    breaks = function(lims) {
      rng <- range(lims, na.rm = TRUE)
      pmin <- floor(log10(rng[1]))
      pmax <- ceiling(log10(rng[2]))
      10^seq(pmin, pmax)
    },
    labels = function(b) {
      exps <- log10(b)
      as.expression(lapply(exps, function(e) bquote(10^.(e))))
    }
  ) +
  # x-axis with float labels
  scale_x_continuous(
    limits = c(0, 16.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    #title=expression("Allele frequency  ("*alpha*"=0.001)"),
    x=expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")")),
    y=expression(r/alpha),
    fill=expression(P[A])
  ) +
  theme_linedraw() +
  theme(
    text=element_text(family = "Latin Modern Math",size=20), 
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20), 
    plot.title=element_text(size=20), 
    legend.text=element_text(size=18), 
    legend.title=element_text(size=20, vjust=0.85),
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 31, barheight = 1))


ggsave("../Fig/freq_recessive.png", width = 23, height = 20, units = "cm")



dfreces %>% filter(r<=0.49 & m==0 & r>=0.0000001 & rs>0.001) %>%
  ggplot(aes(x = ga, y=rs, fill=LD)) + 
  geom_tile() +
  scale_fill_gradientn(colors = rocket(n=10,end = 0.7,direction =1)) +
  geom_contour(aes(z= LD), colour="#FFF", binwidth = 0.02, linetype = "dashed", alpha=0.6)+
  # y-axis with automatic powers of 10
  scale_y_continuous(
    trans = 'log10',
    breaks = function(lims) {
      rng <- range(lims, na.rm = TRUE)
      pmin <- floor(log10(rng[1]))
      pmax <- ceiling(log10(rng[2]))
      10^seq(pmin, pmax)
    },
    labels = function(b) {
      exps <- log10(b)
      as.expression(lapply(exps, function(e) bquote(10^.(e))))
    }
  ) +
  # x-axis with float labels
  scale_x_continuous(
    limits = c(0, 16.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  labs(
    #title=expression("Linkage disequilibrium  ("*alpha*"=0.001)"),
    x=expression(paste("Relative strength of DMIs ("*gamma*"/"*alpha*")")),
    y=expression(r/alpha),
    fill=expression(D)
  ) +
  theme_linedraw() +
  theme(
    text=element_text(family = "Latin Modern Math",size=20), 
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20), 
    plot.title=element_text(size=20), 
    legend.text=element_text(size=18), 
    legend.title=element_text(size=20, vjust=0.85),
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 31, barheight = 1))

ggsave("../Fig/LD_recessive.png", width = 23, height = 20, units = "cm")



##

## r --> 0 
a_val <-0.001

dfcodominant01 %>% filter(r %in% c(1e-6, 1e-5, 1e-4)) %>%
  ggplot(aes(x=as.numeric(ga), y=LD, colour = as.factor(r))) +
  #geom_line(linetype="dashed") +
  geom_point(size = 2.5, alpha=0.75) +
  geom_line(data = \(d) d |> filter(r %in% c(1e-6, 1e-5, 1e-4)),
    aes(y = (-(4*a_val^2)/(4*a_val+ga*a_val)^2) + ((8*a_val + 3*ga*a_val)*8*r*(a_val^2))/(ga*a_val*(4*a_val+ga*a_val)^3), 
                group = r, colour = as.factor(r)), linetype = "solid", linewidth = 1) +
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
  ylab(expression("Linkage disequilibrium (D)")) + 
  labs(color="r") +
  # x-axis with float labels
  scale_x_continuous(
    limits = c(0, 4.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_colour_discrete(
    palette=scales::pal_brewer(palette = "Dark2"),
    name = expression(italic(r)),
    labels = function(x) parse(text = gsub("1e", "10^", x))
  ) +
  theme_linedraw() +
  theme(
    text=element_text(family = "Latin Modern Math",size=20), 
    axis.text=element_text(size=18), 
    axis.title=element_text(size=20), 
    plot.title=element_text(size=20), 
    legend.text=element_text(size=18), 
    legend.title=element_text(size=20),
    legend.position = c(0.9,0.05),
    legend.box.background = element_rect(colour = "black"),
    legend.direction = "horizontal",  
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

ggsave("../Fig/Fig_1C.png", width = 46, height = 15, units = "cm")

#pA <- ((2*a_val)/(4*a_val + ga*a_val)) - (4*r*a_val*(2*a_val + ga*a_val))/(ga*a_val*(4*a_val+ga*a_val)^2)
#D <- (-(4*a_val^2)/(4*a_val+ga*a_val)^2) + ((8*a_val + 3*ga*a_val)*8*r*(a_val^2))/(ga*a_val*(4*a_val+ga*a_val)^3)


## QLE
dfcodominant01 %>% filter(r %in% c(1e-1, 5e-1, 10^(-0.375), 10^(-0.5))) %>%
  ggplot(aes(x=ga, y=LD, colour = as.factor(r))) +
  #geom_line(linetype="dashed") +
  geom_point(size = 2) +
  geom_line(data = \(d) d |> filter(r %in% c(5e-1, 10^(-0.375), 10^(-0.5))),
            aes(y = -((a_val^2) * ga*a_val * (a_val + ga*a_val)^2) / (2*(2*a_val*ga*a_val)*(a_val*((ga*a_val)^3) + r*(2*a_val + ga*a_val)^3)), 
                group = r, colour = as.factor(r)), linetype = "dashed") +
  scale_colour_discrete(palette=scales::pal_brewer(palette = "Dark2")) +
  theme_bw() +
  labs(title="Linkage disequilibrium as a function of selection") + 
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
  theme(legend.position = "bottom")

# D <- 

library(scales)
a_val <- 0.001

dfcodominant001 %>%
  filter(r %in% c(1e-7, 1e-6, 1e-5)) %>%
  ggplot(aes(x = as.numeric(ga), y = LD, colour = as.factor(r))) +
  
  # Numerical (points)
  geom_point(aes(linetype = "Numerical", shape = "Numerical"), 
             size = 2.5, alpha = 0.75) +
  
  # Approximation (lines)
  geom_line(
    data = function(d) d %>% filter(r %in% c(1e-7, 1e-6, 1e-5)),
    aes(
      y = (-(4*a_val^2)/(4*a_val+ga*a_val)^2) + 
        ((8*a_val + 3*ga*a_val)*8*r*(a_val^2))/(ga*a_val*(4*a_val+ga*a_val)^3),
      group = r,
      linetype = "Approximation",
      shape = "Approximation"
    ),
    linewidth = 1
  ) +
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) +
  ylab(expression("Linkage disequilibrium (D)")) +
  scale_x_continuous(
    limits = c(0, 4.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_colour_brewer(
    palette = "Dark2",
    name = expression(r),
    labels = function(x) parse(text = gsub("1e", "10^", x))
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Numerical" = "blank", "Approximation" = "solid"),
    guide = "legend"
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Numerical" = 16, "Approximation" = NA)
  ) +
  
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2),
    shape = guide_legend(order = 2)
  ) +
  
  theme_linedraw() +
  theme(
    text = element_text(family = "Latin Modern Math", size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.position = c(0.98, 0.05),
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    legend.box.just = "right",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )


ggsave("../Fig/Fig_1C.png", width = 46, height = 15, units = "cm")

