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

setwd(dirname(getActiveDocumentContext()$path)) # Works only if executing code manually, not sourcing


#dfcodom <- read.csv("codominant_sols.csv", header = F, sep = ",", col.names = c("ga", "r", "m","nsols"))
dfcodominant <- read.csv("../Data/codominant_sols_freqs.csv", header = F, sep = ",", col.names = c("ga", "r", "m","pA", "LD"))
dfreces <- read.csv("../Data/recessive_sols.csv", header = F, sep = ",", col.names = c("ga", "r", "m","nsols"))

dfcodominant$r2 <- (dfcodominant$LD**2)/((dfcodominant$pA**2)*((1-dfcodominant$pA)**2))

dfreces <- dfreces %>% filter(r>0)
#df1 <- read.csv("solscodominant.csv", header = F, sep = ",", col.names = c("ga", "r", "nsols"))
#df2 <- read.csv("solscodominant1.csv", header = F, sep = ",", col.names = c("ga", "r", "nsols"))
#df3 <- read.csv("solscodominant2.csv", header = F, sep = ",", col.names = c("ga", "r", "nsols"))


dfcodominant %>% filter(r<=0.1 & m==0 & ga<=8 & r>=0.000001) %>%
  ggplot(aes(x = ga, y=r, fill=pA)) + 
  geom_tile() +
  scale_fill_gradientn("Freq. A",colors = hcl.colors(30, "ag_GrnYl", rev=T)) +
  scale_y_continuous(trans='log10') +
  #facet_grid(vars(m)) +
  labs(title="Allele frequency") + 
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
  #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
  ylab("Recombination rate") +
  labs(fill=expression(P[A])) +
  #guides(y.sec = guide_none("migration rate")) +
  theme_bw() +
  theme(text=element_text(size=18), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        plot.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16)) 

#dfcodominant %>% filter(r<=0.1 & m<=0.005 & ga<=8) %>%
dfcodominant %>% filter(r<=0.1 & m==0 & ga<=8 & r>=0.000001) %>%
  ggplot(aes(x = ga, y=r, fill=LD)) + 
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(30, "BuPu", rev=T)) +
  scale_y_continuous(trans='log10') +
  #facet_grid(vars(m)) +
  labs(title="Linkage disequilibrium") + 
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
  #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
  ylab("Recombination rate") +
  #guides(y.sec = guide_none("migration rate")) +
  theme_bw() +
  theme(text=element_text(size=18), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        plot.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16)) 


## Plotting r2
dfcodominant %>% filter(r<=0.1 & m==0 & ga<=8 & r>=0.000001) %>%
  ggplot(aes(x = ga, y=r, fill=r2)) + 
  geom_tile() +
  scale_fill_gradientn(colors = hcl.colors(30, "BuPu", rev=F)) +
  scale_y_continuous(trans='log10') +
  #facet_grid(vars(m)) +
  labs(title="r^2") + 
  xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
  #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
  ylab("Recombination rate") +
  #guides(y.sec = guide_none("migration rate")) +
  theme_bw() +
  theme(text=element_text(size=18), 
        axis.text=element_text(size=16), 
        axis.title=element_text(size=18), 
        plot.title=element_text(size=18), 
        legend.text=element_text(size=16), 
        legend.title=element_text(size=16)) 



# Old code
# dffull %>% filter(r<=0.1 & m<=0.005 & ga<=8) %>%
#   ggplot(aes(x = ga, y=r, fill=factor(nsols))) + 
#   geom_tile() +
#   scale_y_continuous(trans='log10') +
#   scale_fill_brewer('# of polymorphic\nequilibria', palette = "YlOrRd", labels = c("0 (No polymorphism)", "1 (Hybrids)", "2"))+
#   facet_grid(vars(m)) +
#   labs(title="Stable polymorphisms - codominant incompatibilities") + 
#   xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
#   #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
#   ylab("Recombination (log)") +
#   guides(y.sec = guide_none("migration rate")) +
#   theme_bw() +
#   theme(text=element_text(size=14), 
#         axis.text=element_text(size=12), 
#         axis.title=element_text(size=14), 
#         plot.title=element_text(size=14), 
#         legend.text=element_text(size=12), 
#         legend.title=element_text(size=12))  
# 
# 
# 
# dfcodom %>% filter(r<=0.1 & m==0) %>%
#   ggplot(aes(x = ga, y=r, fill=factor(nsols))) + 
#   geom_tile() +
#   scale_y_continuous(trans='log10') +
#   scale_fill_brewer('# of polymorphic\nequilibria', palette = "YlOrRd")+
#   labs(title="Stable polymorphisms - codominant incompatibilities (m=0)") + 
#   xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
#   #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
#   ylab("Recombination (log)") +
#   theme_bw() +
#   theme(text=element_text(size=20), 
#         axis.text=element_text(size=12), 
#         axis.title=element_text(size=20), 
#         plot.title=element_text(size=20), 
#         legend.text=element_text(size=12), 
#         legend.title=element_text(size=12))  



# dfreces %>% filter(r<=0.1 & m<=0.005) %>%
#   ggplot(aes(x = ga, y=r, fill=factor(nsols))) + 
#   geom_tile() +
#   scale_y_continuous(trans='log10') +
#   scale_fill_brewer('# of polymorphic equilibria', palette = "YlOrRd")+
#   facet_grid(vars(m)) +
#   labs(title="Stable polymorphisms - recessive incompatibilities") + 
#   xlab(expression("Relative strength of DMIs ("*gamma*"/"*alpha*")")) + 
#   #xlab(expression(gamma*"/"*alpha*" ("*alpha*"=0.02)")) + 
#   ylab("Recombination (log)") +
#   guides(y.sec = guide_none("migration rate")) +
#   theme_bw() 

