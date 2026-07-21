###############################################################################
# plot_neutral_invertedaxis.r
# 10.2025
#
# Script to generate plots of neutral heterozygosity linked to a DMI and overdominant loci
###############################################################################

library(ggplot2)
library(dplyr)
library(rstudioapi)
library(viridis)
library(grid)
library(scales)
library(cowplot)
library(stringr)

options(warn=1)


########### Loading dataframes
setwd(dirname(getActiveDocumentContext()$path)) # Change dir to file. OBS: Might not work if sourcing.


alldfs = data_frame()
alldfs_intra = data_frame()

for (filename in list.files(path="../Output_new/", pattern="^allele_ABN*")) {
  dfalleles = read.csv(paste0("../Output_new/",filename))
  splitname = strsplit(filename, "_")[[1]]
  dftype = ifelse(any(c("codom", "asym", "overdom", "recessive") %in% splitname), splitname[3], "codom+overdom")
  index_offset = as.integer(any(c("codom", "asym", "overdom", "recessive") %in% splitname))
  selcoeff = strsplit(splitname[3+index_offset], "=")[[1]]
  if (selcoeff[2]=="0a" & dftype=="codom+overdom"){
    g=0
    a=as.numeric(selcoeff[3])
    dftype = "overdom" # Change type to only overdominance
  }
  else{
    g = ifelse("overdom" %in% splitname, 0, as.numeric(selcoeff[3]))
    a = ifelse(any(c("codom", "overdom") %in% splitname), ifelse("overdom" %in% splitname, as.numeric(selcoeff[3]), 0), g/as.numeric(str_replace(selcoeff[2], "a", "")))  
  }
  a = ifelse(abs(a)==Inf,0,a)
  b = ifelse("asym" %in% splitname, 0, a)
  pA0=as.numeric(strsplit(splitname[4+index_offset], "=")[[1]][2])
  pB0=1-pA0
  dfalleles$g=g
  dfalleles$a=a
  dfalleles$b=b
  dfalleles$pA0=pA0
  dfalleles$pB0=pB0
  dfalleles$type = dftype
  dfalleles$filename = filename
  alldfs = rbind(alldfs,dfalleles)
}


for (filename in list.files(path="../Output_new/", pattern="^allele_ANB*")) {
  dfalleles = read.csv(paste0("../Output_new/",filename))
  splitname = strsplit(filename, "_")[[1]]
  dftype = ifelse(any(c("codom", "asym", "overdom", "recessive") %in% splitname), splitname[3], "codom+overdom")
  index_offset = as.integer(any(c("codom", "asym", "overdom", "recessive") %in% splitname))
  selcoeff = strsplit(splitname[3+index_offset], "=")[[1]]
  if (selcoeff[2]=="0a" & dftype=="codom+overdom"){
    g=0
    a=as.numeric(selcoeff[3])
    dftype="overdom"
  }
  else{
    g = ifelse("overdom" %in% splitname, 0, as.numeric(selcoeff[3]))
    a = ifelse(any(c("codom", "overdom") %in% splitname), ifelse("overdom" %in% splitname, as.numeric(selcoeff[3]), 0), g/as.numeric(str_replace(selcoeff[2], "a", "")))  
  }
  a = ifelse(abs(a)==Inf,0,a)
  b = ifelse("asym" %in% splitname, 0, a) # TODO: Change to arbitrary b in name
  pA0=as.numeric(strsplit(splitname[4+index_offset], "=")[[1]][2])
  pB0=1-pA0
  dfalleles$g=g
  dfalleles$a=a
  dfalleles$b=b
  dfalleles$pA0=pA0
  dfalleles$pB0=pB0
  dfalleles$type = dftype
  dfalleles$filename = filename
  alldfs_intra = rbind(alldfs_intra,dfalleles)
}


alldfs <- na.omit(alldfs)
alldfs_intra <- na.omit(alldfs_intra)

# Relative heterozygosity (Not used)
alldfs$rel_het <- alldfs$het/(1- (alldfs$pA0^2 + alldfs$pB0^2))
alldfs_intra$rel_het <- alldfs_intra$het/(1- (alldfs_intra$pA0^2 + alldfs_intra$pB0^2))



# log-like transform centered at 0
log_centered_trans <- function(x0 = 1e-4, base = 10) {
  trans_new(
    name = "log_centered",
    transform = function(x) sign(x) * log10(1 + abs(x)/x0),
    inverse = function(x) sign(x) * x0 * (10^abs(x) - 1),
    domain = c(-Inf, Inf)
  )
}


################ Plotting

# Setting params for plot
# pA_init <- 0.7
# pB_init <- 0.3
# a_value <- 0.02
# b_value <- 0.02
# g_values <- c(0.02,0.04,0.06) #c(0)
# r_ab_values <- c(1e-1,1e-2,1e-6)
# architecture_type="codom+overdom"

## params for plot with OD only
# pA_init <- 0.7
# pB_init <- 0.3
# a_value <- 0.02
# b_value <- 0.02
# g_values <- c(0)
# r_ab_values <- c(1e-1,1e-2,1e-6)
# architecture_type="overdom"

### Params for plot with DMIs only
pA_init <- 0.7
pB_init <- 0.3
a_value <- 0
b_value <- 0
g_values <- c(0.02,0.04,0.06) #c(0)
r_ab_values <- c(1e-1,1e-6,1e-2)
architecture_type="codom"

yaxis="het"
###########################

dfN <- alldfs %>%
  mutate(rn_signed = if_else(pA0 == pA_init, r_bc, -r_bc),
         freq = if_else(pA0==pA_init, freq, 1-freq)) %>%
  filter(pA0 %in% c(pA_init,pB_init), a==a_value, b==b_value, g %in% g_values, r_ab %in% r_ab_values, allele=="N",
         type==architecture_type)


# Get eq frequencies of selected alleles
# Filtered only for allele A, to get both sides when frequencies are equal.
ABvalues <- alldfs_intra %>%
  filter(allele %in% c("A"), a==a_value, b==b_value, g %in% g_values, r_ab %in% r_ab_values, pA0 %in% c(pA_init,pB_init), pos>=0.999 | pos<=0.001, type==architecture_type) %>%
  group_by(r_ab, allele, g, pA0) %>%
  summarise(het = mean(het, na.rm = TRUE),
            freq= mean(freq, na.rm = TRUE),
            rel_het = mean(rel_het, na.rm = TRUE),
            .groups="drop") %>%
  mutate(pos = if_else(pA0==pA_init, 0, 1),
         freq = if_else(pA0==pA_init, freq, 1-freq))

# Duplicate if pA_init=0.5
if (pA_init==0.5){
  dfN_dup <- dfN %>%
    mutate(rn_signed = -r_bc)
  dfN <- rbind(dfN, dfN_dup)
  ABvalues <- rbind(ABvalues, ABvalues%>%mutate(pos = 0))  
}


# ----------------------------
# First row
# ----------------------------
facet_vals_p1 <- g_values #unique(dfN$r_ab)

plot_list_p1 <- lapply(seq_along(facet_vals_p1), function(i) {
  df_sub <- dfN %>% filter(rn_signed < 0, g == facet_vals_p1[i])
  
  p <- ggplot(df_sub, aes_string(x = "r_bc", y = yaxis, colour = "as.factor(r_ab)")) +
    geom_point(size = 1) +
    geom_line(aes(group = interaction(r_ab,pA0)), alpha=0.4) +
    geom_hline(yintercept = 0, colour = "black", linetype = "solid", linewidth = 0.2) +
    scale_x_continuous(
      limits = c(min(dfN$r_bc), 0.5),
      breaks = c(1e-5, 1e-3, 1e-1),
      labels = c(expression(10^-5), expression(10^-3), expression(10^-1)),
      transform = c("log10", "reverse")
    ) +
    scale_color_manual(
      values = c(
        "#4b597d",
        "#59ced8",
        "#5df5f6"
      )) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0),
      limits = c(-0.05, case_when(
        yaxis=="het" ~ 0.5,
        yaxis=="freq" ~ 1.0,
        yaxis=="rel_het" ~ 2.0,
        TRUE ~ 1.0
      ))
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      axis.title.y = element_blank(),
      #axis.text.y = element_blank(),
      #axis.ticks.y = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.2),
      axis.ticks.length = unit(0.1, "cm"),
      # axis.title.x = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    coord_cartesian(clip = "off")
  
  # Hide x-axis for all but last facet
  if (i < length(facet_vals_p1)) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  } else {
    p <- p + geom_text(
      data = data.frame(x = max(dfN$rn_signed, na.rm = TRUE), y = 0),
      aes(x = x, y = y),
      label = "N ←",
      fontface = "bold",
      hjust = 0.5, vjust = 1.8,
      inherit.aes = FALSE
    ) +
      xlab(expression(r[an]))
    
  }
  
  p
})

# ----------------------------
# Third row
# ----------------------------
facet_vals_p3 <- g_values #unique(dfN$r_ab)

plot_list_p3 <- lapply(seq_along(facet_vals_p3), function(i) {
  df_sub <- dfN %>% filter(rn_signed > 0, g == facet_vals_p3[i])
  
  p <- ggplot(df_sub, aes_string(x = "r_bc", y = yaxis, colour = "as.factor(r_ab)")) +
    geom_point(size = 1) +
    geom_line(aes(group = interaction(r_ab,pA0)), alpha=0.4) +
    geom_hline(yintercept = 0, colour = "black", linetype = "solid", linewidth = 0.2) +
    scale_x_continuous(
      limits = c(min(dfN$r_bc), 0.65),
      breaks = c(1e-5, 1e-3, 1e-1),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      # labels = c(expression(10^-5), expression(10^-3), expression(10^-1)),
      trans = "log10"
    ) +
    scale_color_manual(
      values = c(
        "#4b597d",
        "#59ced8",
        "#5df5f6"
      )) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0),
      limits = c(-0.05, case_when(
        yaxis=="het" ~ 0.5,
        yaxis=="freq" ~ 1.0,
        yaxis=="rel_het" ~ 2.0,
        TRUE ~ 1.0
      )),
      sec.axis = sec_axis(~ ., name = expression(gamma), breaks = NULL, labels = NULL)
    ) +
    geom_text(
      data = data.frame(x = max(dfN$r_bc)+0.1, y = 0.4),
      aes(x = x, y = y),
      label = paste0("gamma == ",facet_vals_p3[i], " ~ ' ' ~alpha==", a_value),
      parse = T,
      inherit.aes = FALSE,
      hjust = 0, vjust = -1, angle=270
    ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.2),
      axis.ticks.length = unit(0.1, "cm"),
      #axis.title.x = element_blank(),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    ) +
    coord_cartesian(clip = "off")
  
  if (i < length(facet_vals_p3)) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  } else {
    p <- p + geom_text(
      data = data.frame(x = max(dfN$rn_signed, na.rm = TRUE), y = 0),
      aes(x = x, y = y),
      label = "→ N",
      fontface = "bold",
      hjust = 0.5, vjust = 1.8,
      inherit.aes = FALSE
    )  +
      xlab(expression(r[bn]))
  }
  
  p
})

# ----------------------------
# Center row
# ----------------------------
facet_vals_p2 <- g_values # unique(alldfs_intra$r_ab)

plot_list_p2 <- lapply(seq_along(facet_vals_p2), function(i) {
  df_sub <- alldfs_intra %>% filter(pA0 == pA_init, a == a_value, b==b_value, g %in% g_values,  r_ab %in% r_ab_values, allele == 'N', g == facet_vals_p2[i], type==architecture_type)
  
  p <- ggplot(df_sub, aes_string(x = "pos", y = yaxis, colour = "as.factor(r_ab)")) +
    geom_point(size = 1) +
    geom_line(aes(group = interaction(r_ab,pA0)), alpha=0.4) +
    geom_hline(yintercept = 0, colour = "black", linetype = "solid", linewidth = 0.2) +
    annotate("segment", x = 0, xend = 0, y = 0, yend = 0.5,
             colour = "black", linetype = "dotted", linewidth = 0.2) +
    annotate("segment", x = 1, xend = 1, y = 0, yend = 0.5,
             colour = "black", linetype = "dotted", linewidth = 0.2) +
    geom_point(
      data = ABvalues %>% filter(r_ab==facet_vals_p2[i]),
      aes_string(x = "pos", y = yaxis, colour = "as.factor(r_ab)"),
      inherit.aes = FALSE, size = 3
    ) +
    scale_color_manual(
      values = c(
        "#4b597d",
        "#59ced8",
        "#5df5f6"
      )) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0),
      limits = c(-0.05, case_when(
        yaxis=="het" ~ 0.5,
        yaxis=="freq" ~ 1.0,
        yaxis=="rel_het" ~ 2.0,
        TRUE ~ 1.0
      ))
    ) +
    # scale_x_continuous(
    #   breaks = c(0,0.5,1.0),
    #   labels = c(expression(0),expression(0.5*10^0), expression(10^0))
    # ) +
    theme_minimal(base_family = "Times New Roman") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks = element_line(colour = "black", linewidth = 0.2),
      axis.ticks.length = unit(0.1, "cm"),
      strip.text   = element_blank(),
      legend.position="none",
      # axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = unit(c(0, 0, 0.055, 0), "cm")
    ) +
    coord_cartesian(clip = "off")
  
  if (i < length(facet_vals_p2)) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  } else {
    p <- p + geom_text(
      data = data.frame(x = 1, y = 0),
      aes(x = x, y = y),
      label = "B",
      fontface = "bold",
      hjust = 0.5, vjust = 1.8,
      inherit.aes = FALSE
    ) +
      geom_text(
        data = data.frame(x = 0, y = 0),
        aes(x = x, y = y),
        label = "A",
        fontface = "bold",
        hjust = 0.5, vjust = 1.8,
        inherit.aes = FALSE
      )  +
      xlab("Relative position of N")
  }
  p
})


# Extract the legend from one plot (all plots share the same color scale)
legend <- get_legend(
  plot_list_p1[[1]] + 
    theme(legend.position = "bottom") +
    labs(colour = expression(r[ab]))
)

# Build each column (stacked vertically)
col1 <- plot_grid(plotlist = plot_list_p1, ncol = 1, align = "v")
col2 <- plot_grid(plotlist = plot_list_p2, ncol = 1, align = "v")
col3 <- plot_grid(plotlist = plot_list_p3, ncol = 1, align = "v")

# Combine columns horizontally
main_plot <- plot_grid(col1, col2, col3, ncol = 3, align = "h")


legend <- get_legend(
  plot_list_p2[[1]] + 
    theme(legend.position = "bottom") +
    labs(colour = expression(r[ab]))  # legend title
)

main_plot_shrunk <- plot_grid(main_plot, ncol = 1, scale = 0.93)

main_plot_labeled <- ggdraw(main_plot_shrunk) +
  draw_label("Heterozygosity (H)", x = 0.02, y = 0.5, angle = 90, vjust = 0.5, fontfamily = "Times New Roman")

final_plot <- plot_grid(main_plot_labeled, legend, ncol = 1, rel_heights = c(1, 0.05))

# Draw the final figure
print(final_plot)

#figname <- paste0("../Fig/",architecture_type,"a=",as.character(a_value),"pA=",as.character(pA_init),"_swappedlegend.pdf")
#save_plot(figname, final_plot, base_height = 6, base_width = 8)


#save_plot("../Fig/od_dmis_codominant_a=0.01.svg", final_plot, base_height = 6, base_width = 8)
#save_plot("../Fig/od_dmis_codominant_a=0.01.pdf", final_plot, base_height = 6, base_width = 8)
#save_plot("../Fig/od_dmis_codominant_a=0.02_pA=0.6.pdf", final_plot, base_height = 6, base_width = 8)
#save_plot("../Fig/od_dmis_codominant_a=0.02_pA=0.7.pdf", final_plot, base_height = 6, base_width = 8)

#save_plot("../Fig/od_dmis_recessive_a=0.01.svg", final_plot, base_height = 6, base_width = 8)
#save_plot("../Fig/od_dmis_recessive_a=0.01.pdf", final_plot, base_height = 6, base_width = 8)

#save_plot("../Fig/dmis_codominant.svg", final_plot, base_height = 6, base_width = 8)
#save_plot("../Fig/dmis_codominant.pdf", final_plot, base_height = 6, base_width = 8)



#########################
# This plot has to be generated in two steps, one for overdominance and one for DMIs only. 
###### overdominance plot
# Get plot for:
# alpha = 0.02; r_ab=0.01; g = c(0)
# row1 <- plot_grid(plotlist = c(plot_list_p1[[1]],plot_list_p2[[1]],plot_list_p3[[1]]), ncol = 3, align = "h")
# # row2 <- plot_grid(plotlist = c(plot_list_p1[[2]],plot_list_p2[[2]],plot_list_p3[[2]]), ncol = 3, align = "h")
# # row3 <- plot_grid(plotlist = c(plot_list_p1[[3]],plot_list_p2[[3]],plot_list_p3[[3]]), ncol = 3, align = "h")
# od_plot <- row1
# od_plot
###### recessive DMIs plot
dmis_plot <- main_plot
#
dmis_plot <- plot_grid(plotlist = c(plot_list_p1[[2]], plot_list_p2[[2]], plot_list_p3[[2]]), ncol = 3, align = "v")
row1 <- plot_grid(plotlist = c(plot_list_p1[[1]],plot_list_p2[[1]],plot_list_p3[[1]]), ncol = 3, align = "h")
row2 <- plot_grid(plotlist = c(plot_list_p1[[2]],plot_list_p2[[2]],plot_list_p3[[2]]), ncol = 3, align = "h")
row3 <- plot_grid(plotlist = c(plot_list_p1[[3]],plot_list_p2[[3]],plot_list_p3[[3]]), ncol = 3, align = "h")

dmis_plot <- row3

## Joint plot
od_dmis_plot <- plot_grid(od_plot, dmis_plot, rel_heights = c(0.8,1), ncol = 1)

od_dmis_plot_shrunk <- plot_grid(od_dmis_plot, ncol = 1, scale = 0.93)


od_dmis_plot_labeled <- ggdraw(od_dmis_plot_shrunk) +
   draw_label("Heterozygosity (H)", x = 0.02, y = 0.5, angle = 90, vjust = 0.5, fontfamily = "Times New Roman")


legend <- get_legend(
  plot_list_p2[[1]] +
    theme(legend.position = "bottom") +
    labs(colour = expression(r[ab]))  # legend title
)

od_dmis_plot_2 <- plot_grid(od_dmis_plot_labeled, legend, ncol = 1, rel_heights = c(1, 0.05))
#
print(od_dmis_plot_2)

#save_plot("../Fig/separate_od_DMIs_a=0.01.svg", od_dmis_plot_2, base_height = 4, base_width = 8)
#save_plot("../Fig/separate_od_DMIs_a=0.01.pdf", od_dmis_plot_2, base_height = 4, base_width = 8)

#save_plot("../Fig/separate_od_DMIs_a=0.02_pA=0.7.svg", od_dmis_plot_2, base_height = 4, base_width = 8)
save_plot("../Fig/separate_od_DMIs_a=0.02_pA=0.7_inv_r_ab.pdf", od_dmis_plot_2, base_height = 4, base_width = 8)


###########

# 
# row_plot_shrunk <- plot_grid(row_plot, ncol = 1, scale = 0.93)
# 
# row_plot_labeled <- ggdraw(row_plot_shrunk) +
#   draw_label("Heterozygosity (H)", x = 0.02, y = 0.5, angle = 90, vjust = 0.5)
# 
# row_plot <- plot_grid(row_plot_labeled, legend, ncol = 1, rel_heights = c(1, 0.05))
# 
# print(row_plot)
# # 
