###############################################################################
# plot_neutral.r
# 05.05.2026
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
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)

options(warn=1)

##############################################################################
##############################################################################
# Config and parameters
##############################################################################
##############################################################################

plot_type <- "combined"   # "combined" for DMIs+overdominance, or "separate" for OD and DMIs in separate panels
save_plot <- F
file_prefix <- str_replace_all(as.character(Sys.Date()), "-", "")

yaxis <- "het"
# Parameters per plot type

params_combined <- list(
  pA_init = 0.7,
  pB_init = 0.3,
  a_value = 0.02,
  b_value = 0.02,
  g_values  = c(0.02, 0.04, 0.06),
  r_ab_values = c(1e-1, 1e-2, 1e-6),
  architecture_type = "codom+overdom"
)

# Recessive DMIs + OD
# params_combined <- list(
#   pA_init = 0.7,
#   pB_init = 0.3,
#   a_value = 0.01,
#   b_value = 0.01,
#   g_values  = c(0.02, 0.04, 0.06),
#   r_ab_values = c(1e-1, 1e-2, 1e-6),
#   architecture_type = "recessive"
# )

# Only Codominant DMIs 
# params_combined <- list(
#   pA_init = 0.7,
#   pB_init = 0.3,
#   a_value = 0,
#   b_value = 0,
#   g_values  = c(0.02, 0.04, 0.06),
#   r_ab_values = c(1e-1, 1e-2, 1e-6),
#   architecture_type = "codom"
# )

# Codom+DMIs pA=0.5
# params_combined <- list(
#   pA_init = 0.5,
#   pB_init = 0.5,
#   a_value = 0.02,
#   b_value = 0.02,
#   g_values  = c(0.02, 0.04, 0.06),
#   r_ab_values = c(1e-1, 1e-2, 1e-6),
#   architecture_type = "codom+overdom"
# )

params_od <- list(
  pA_init = 0.7,
  pB_init = 0.3,
  a_value = 0.02,
  b_value = 0.02,
  g_values  = c(0),
  r_ab_values = c(1e-1, 1e-2, 1e-6),
  architecture_type = "overdom",
  row_idx = 2        # row from r_ab_values to plot
)

params_dmis <- list( # DMIs only
  pA_init = 0.7,
  pB_init = 0.3,
  a_value = 0,       
  b_value = 0,
  g_values  = c(0.02, 0.04, 0.06),
  r_ab_values = c(1e-1, 1e-2, 1e-6),
  architecture_type = "codom",
  row_idxs = c(1)  # which r_ab rows to stack
)


##############################################################################
##############################################################################


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

# build the three plot columns for one parameter set. Uses alldfs, and alldfs_intra loaded before

build_plot_lists <- function(params) {
  
  pA_init <- params$pA_init
  pB_init <- params$pB_init
  a_value <- params$a_value
  b_value <- params$b_value
  g_values <- params$g_values
  r_ab_values <- params$r_ab_values
  architecture_type <- params$architecture_type
  
  # dfN: neutral locus linked to A or B 
  dfN <- alldfs %>%
    mutate(
      rn_signed = if_else(pA0 == pA_init,  r_bc, -r_bc),
      freq = if_else(pA0 == pA_init, freq, 1 - freq)
    ) %>%
    filter(
      pA0 %in% c(pA_init, pB_init),
      a == a_value,
      b == b_value,
      g %in% g_values,
      r_ab %in% r_ab_values,
      allele == "N",
      type == architecture_type
    )
  

  # ABvalues - equilibrium frequencies of selected alleles
  ABvalues <- alldfs_intra %>%
    filter(
      allele %in% c("A"),
      a == a_value,
      b == b_value,
      g %in% g_values,
      r_ab %in% r_ab_values,
      pA0 %in% c(pA_init, pB_init),
      pos >= 0.999 | pos <= 0.001,
      type == architecture_type
    ) %>%
    group_by(r_ab, allele, g, pA0) %>%
    summarise(
      het = mean(het,     na.rm = TRUE),
      freq = mean(freq,    na.rm = TRUE),
      rel_het = mean(rel_het, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pos = if_else(pA0 == pA_init, 0, 1),
      freq = if_else(pA0 == pA_init, freq, 1 - freq)
    )
  
  # Duplicate both sides when pA_init == 0.5 (symmetric case)
  if (pA_init == 0.5) {
    dfN <- rbind(dfN, dfN %>% mutate(rn_signed = -r_bc))
    ABvalues <- rbind(ABvalues, ABvalues %>% mutate(pos = 0))
  }
  
  # warn if either side of dfN is empty (for pA0 mismatches)
  # n_left  <- sum(dfN$rn_signed < 0, na.rm = TRUE)
  # n_right <- sum(dfN$rn_signed > 0, na.rm = TRUE)
  # if (n_left == 0)  warning(sprintf("build_plot_lists [%s]: LEFT column (rn_signed<0) is empty. Check that pA0=%.1f rows exist in alldfs.", architecture_type, pB_init))
  # if (n_right == 0) warning(sprintf("build_plot_lists [%s]: RIGHT column (rn_signed>0) is empty. Check that pA0=%.1f rows exist in alldfs.", architecture_type, pA_init))
  
  # Pre-compute x ranges for each side separately (avoids label/scale errors)
  r_bc_left  <- dfN %>% filter(rn_signed < 0) %>% pull(r_bc)
  r_bc_right <- dfN %>% filter(rn_signed > 0) %>% pull(r_bc)
  x_min_left  <- if (length(r_bc_left)  > 0) min(r_bc_left)  else 1e-6
  x_min_right <- if (length(r_bc_right) > 0) min(r_bc_right) else 1e-6
  
  n_g <- length(g_values)
  colour_vals <- grDevices::colorRampPalette(c("dodgerblue3", "firebrick2"))(n_g)
  ylims <- c(-0.05, switch(yaxis, het = 0.5, freq = 1.0, rel_het = 2.0, 1.0))
  ybrks  <- c(0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
  n_facets <- length(r_ab_values)
  
  # Shared theme snippet to hide x axis on non-last rows
  hide_xaxis <- theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )
  
  
  # Left column (p1): N linked on the A side;  rn_signed < 0
  plot_list_p1 <- lapply(seq_len(n_facets), function(i) {
    df_sub <- dfN %>% filter(rn_signed < 0, r_ab == r_ab_values[i])
    
    p <- ggplot(df_sub, aes(x = r_bc, y = .data[[yaxis]], colour = as.factor(g))) +
    #p <- ggplot(df_sub, aes_string(x = "r_bc", y = yaxis, colour = "as.factor(g)")) +
      geom_point(size = 1) +
      geom_line(aes(group = interaction(g, pA0)), alpha = 0.4) +
      geom_hline(yintercept = 0, colour = "black", linetype = "solid", linewidth = 0.2) +
      scale_x_continuous(
        limits = c(x_min_left, 0.5),
        breaks = c(1e-5, 1e-3, 1e-1),
        labels = c(expression(10^-5), expression(10^-3), expression(10^-1)),
        transform = c("log10", "reverse")
      ) +
      scale_colour_manual(values = colour_vals) +
      scale_y_continuous(breaks = ybrks, limits = ylims) +
      theme_minimal(base_family = "Times New Roman") +
      theme(
        axis.title.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        strip.text = element_blank(),
        strip.background  = element_blank(),
        legend.position  = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        text=element_text(family="Times New Roman")
      ) +
      coord_cartesian(clip = "off")
    
    if (i < n_facets) {
      p <- p + hide_xaxis
    } else {
      p <- p +
        geom_text(
          data = data.frame(x = max(r_bc_left, na.rm = TRUE), y = 0),
          aes(x = x, y = y),
          label = "N \u2190",
          fontface = "bold",
          hjust = 0.5, vjust = 1.8,
          inherit.aes = FALSE
        ) +
        xlab(expression(r[an]))
    }
    p
  })
  
  # Right column (p3): N linked on the B side;  rn_signed > 0
  plot_list_p3 <- lapply(seq_len(n_facets), function(i) {
    df_sub <- dfN %>% filter(rn_signed > 0, r_ab == r_ab_values[i])
    
    #p <- ggplot(df_sub, aes_string(x = "r_bc", y = yaxis, colour = "as.factor(g)")) +
    p <- ggplot(df_sub, aes(x = r_bc, y = .data[[yaxis]], colour = as.factor(g))) +
      geom_point(size = 1) +
      geom_line(aes(group = interaction(g, pA0)), alpha = 0.4) +
      geom_hline(yintercept = 0, colour = "black", linetype = "solid", linewidth = 0.2) +
      scale_x_continuous(
        limits = c(x_min_right, 0.65),
        breaks = c(1e-5, 1e-3, 1e-1),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        trans  = "log10"
      ) +
      scale_colour_manual(values = colour_vals) +
      scale_y_continuous(
        breaks = ybrks,
        limits = ylims,
        sec.axis = sec_axis(~ ., name = expression(gamma), breaks = NULL, labels = NULL)
      ) +
      #fix here
      geom_text(
        data = data.frame(x = max(dfN$r_bc) + 0.1, y = 0.4),
        aes(x = x, y = y),
        #label = bquote(r[ab] == 10^.(log10(r_ab_values[i]))),
        label = paste0("r[ab] == 10^", log10(r_ab_values[i])), 
        parse = TRUE,
        inherit.aes = FALSE,
        hjust = 0, vjust = -1, angle = 270
      ) +
      theme_minimal(base_family = "Times New Roman") +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        strip.text = element_blank(),
        strip.background  = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm")
      ) +
      coord_cartesian(clip = "off")
    
    if (i < n_facets) {
      p <- p + hide_xaxis
    } else {
      p <- p +
        geom_text(
          data = data.frame(x = max(r_bc_right, na.rm = TRUE), y = 0),
          aes(x = x, y = y),
          label = "\u2192 N",
          fontface = "bold",
          hjust = 0.5, vjust = 1.8,
          inherit.aes = FALSE
        ) +
        xlab(expression(r[bn]))
    }
    p
  })
  

  # Centre column (p2): N positioned between A and B (alldfs_intra)
  plot_list_p2 <- lapply(seq_len(n_facets), function(i) {
    df_sub <- alldfs_intra %>%
      filter(
        pA0 == pA_init,
        a == a_value,
        b == b_value,
        g  %in% g_values,
        r_ab %in% r_ab_values,
        allele == "N",
        r_ab == r_ab_values[i],
        type == architecture_type
      )
    
    #p <- ggplot(df_sub, aes_string(x = "pos", y = yaxis, colour = "as.factor(g)")) +
    p <- ggplot(df_sub, aes(x = pos, y = .data[[yaxis]], colour = as.factor(g))) +
      geom_point(size = 1) +
      geom_line(aes(group = interaction(g, pA0)), alpha = 0.4) +
      geom_hline(yintercept = 0, colour = "black", linetype = "solid", linewidth = 0.2) +
      annotate("segment", x = 0, xend = 0, y = 0, yend = 0.5,
               colour = "black", linetype = "dotted", linewidth = 0.2) +
      annotate("segment", x = 1, xend = 1, y = 0, yend = 0.5,
               colour = "black", linetype = "dotted", linewidth = 0.2) +
      geom_point(
        data = ABvalues %>% filter(r_ab == r_ab_values[i]),
        aes(x = pos, y = .data[[yaxis]], colour = as.factor(g)),
        #aes_string(x = "pos", y = yaxis, colour = "as.factor(g)"),
        inherit.aes = FALSE,
        size = 3
      ) +
      scale_colour_manual(values = colour_vals) +
      scale_y_continuous(breaks = ybrks, limits = ylims) +
      theme_minimal(base_family = "Times New Roman") +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.2),
        axis.ticks.length = unit(0.1, "cm"),
        strip.text = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank(),
        plot.margin = unit(c(0, 0, 0.055, 0), "cm"),
        text=element_text(family="Times New Roman")
      ) +
      coord_cartesian(clip = "off")
    
    if (i < n_facets) {
      p <- p + hide_xaxis + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    } else {
      p <- p +
        geom_text(data = data.frame(x = 0, y = 0), aes(x = x, y = y),
                  label = "A", fontface = "bold", hjust = 0.5, vjust = 1.8,
                  inherit.aes = FALSE) +
        geom_text(data = data.frame(x = 1, y = 0), aes(x = x, y = y),
                  label = "B", fontface = "bold", hjust = 0.5, vjust = 1.8,
                  inherit.aes = FALSE) +
        xlab("Relative position of N")
    }
    p
  })
  
  list(p1 = plot_list_p1, p2 = plot_list_p2, p3 = plot_list_p3)
}


# Get horizontal row from plot lists at index i
make_row <- function(lists, i) {
  plot_grid(
    plotlist = list(lists$p1[[i]], lists$p2[[i]], lists$p3[[i]]),
    ncol = 3, align = "h"
  )
}

# Get 3-column grid (all rows stacked)
assemble_grid <- function(lists) {
  col1 <- plot_grid(plotlist = lists$p1, ncol = 1, align = "v")
  col2 <- plot_grid(plotlist = lists$p2, ncol = 1, align = "v", rel_heights = c(1,1,1))
  col3 <- plot_grid(plotlist = lists$p3, ncol = 1, align = "v")
  plot_grid(col1, col2, col3, ncol = 3, align = "h")
}

# extract gamma legend from centre column
make_legend <- function(plot_list_p2) {
  get_legend(
    plot_list_p2[[1]] +
      theme(legend.position = "bottom") +
      labs(colour = expression(gamma))
  )
}

# Wrap main plot with y-axis label and legend
wrap_with_ylabel <- function(main, legend) {
  shrunk  <- plot_grid(main, ncol = 1, scale = 0.93)
  labeled <- ggdraw(shrunk) +
    draw_label("Heterozygosity (H)", x = 0.02, y = 0.5, angle = 90, vjust = 0.5, fontfamily = "Times New Roman") 
  plot_grid(labeled, legend, ncol = 1, rel_heights = c(1, 0.05))
}


#### Executing plots

if (plot_type == "combined") { # For overdom+codom, etc.
  lists <- build_plot_lists(params_combined)
  main <- assemble_grid(lists)
  legend <- make_legend(lists$p2)
  final_plot <- wrap_with_ylabel(main, legend)
  print(final_plot)
  
  if (save_plot) {
    figname <- paste0(
      "../Fig/", file_prefix,
      params_combined$architecture_type,
      "_a=", params_combined$a_value,
      "pA=", params_combined$pA_init,
      ".pdf"
    )
    save_plot(figname, final_plot, base_height = 6, base_width = 8)
    message("Saved: ", figname)
  }
  
} else if (plot_type == "separate") { # For OD only and DMI only
  
  # Build plot lists independently 
  lists_od <- build_plot_lists(params_od)
  lists_dmis <- build_plot_lists(params_dmis)
  
  # OD plot
  od_plot <- make_row(lists_od, params_od$row_idx)
  
  # DMIs plot
  dmis_rows <- lapply(params_dmis$row_idxs, function(i) make_row(lists_dmis, i))
  dmis_plot <- plot_grid(plotlist = dmis_rows, ncol = 1)
  
  combined <- plot_grid(od_plot, dmis_plot, rel_heights = c(0.9, 1), ncol = 1)
  legend <- make_legend(lists_dmis$p2)
  final_plot <- wrap_with_ylabel(combined, legend)
  print(final_plot)
  
  if (save_plot) {
    figname <- paste0(
      "../Fig/", file_prefix,
      "separate_od_DMIs",
      "_a=", params_od$a_value,
      "_pA=", params_od$pA_init,
      "_added_r_ab.pdf"
    )
    save_plot(figname, final_plot, base_height = 6, base_width = 8)
    message("Saved: ", figname)
  }
}
