###############################################################################
# numsols_plots.r
# Julio Ayala Lopez
# 08.2024
#
# Script to generate heatmap plots of allele frequency and linkage
# disequilibrium using numerical solutions of a two locus BDMIs/heterosis ODE.
###############################################################################

library(ggplot2)
library(dplyr)
library(rstudioapi)
library(viridis)
library(ggiraph)
library(scales)
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)


setwd(dirname(getActiveDocumentContext()$path))

################################################################################
########################### Helper functions ###################################
################################################################################

load_eq_df <- function(filepath, alpha) {
  df <- read.csv(filepath, header = FALSE, sep = ",", col.names = c("ga", "r", "m", "pA", "LD"))
  df$r2 <- (df$LD^2) / ((df$pA^2) * ((1 - df$pA)^2))
  df$rs <- df$r / alpha
  return(df)
}

plot_dmi_heatmap <- function(df, r_lim = NULL, rs_lim = NULL, max_ga = 4.0,
                             fill_var = "pA", filename = "", width = 23, height = 20) {
  plot_df <- df %>% filter(m == 0 & ga <= max_ga)
  
  ## Bouund based on filters on r and r/alpha
  if (!is.null(r_lim)) {
    plot_df <- plot_df %>% filter(r >= r_lim[1] & r <= r_lim[2])
  }
  if (!is.null(rs_lim)) {
    plot_df <- plot_df %>% filter(rs >= rs_lim[1] & rs <= rs_lim[2])
  }
  
  if (fill_var == "pA") {
    fill_expr <- expression(P[A])
    is_recessive <- max(plot_df$ga, na.rm = TRUE) > 5
    # Set colour palette 
    if (is_recessive) {
      color_limits <- c(0.18, 0.51)
      color_breaks <- c(0.2, 0.3, 0.4, 0.5)
    } else {
      color_limits <- c(0.24, 0.51)
      color_breaks <- c(0.25, 0.375, 0.5)
    }
    
    fill_scale <- scale_fill_gradientn(
      fill_expr,
      colors = tail(hcl.colors(100, "ag_GrnYl", rev = TRUE), 95),
      limits = color_limits,
      breaks = color_breaks
    )
    contour_mapping <- aes(z = pA)
    
  } else if (fill_var == "LD") {
    fill_expr <- expression(D)
    fill_scale <- scale_fill_gradientn(colors = rocket(n = 10, end = 0.7, direction = 1))
    contour_mapping <- aes(z = LD)
  } else {
    stop("TODO: Add other variables if any (eg r2)")
  }
  
  p <- ggplot(plot_df, aes(x = ga, y = rs, fill = .data[[fill_var]])) +
    geom_tile() +
    fill_scale +
    geom_contour(
      contour_mapping,
      colour = "#FFF",
      binwidth = 0.02,
      linetype = "dashed",
      alpha = 0.6
    ) +
    scale_y_continuous(
      trans = "log10",
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
    scale_x_continuous(
      limits = c(0, max_ga),
      labels = function(x) sprintf("%.1f", x)
    ) +
    labs(
      x = expression(paste("Relative strength of BDMIs (" * gamma * "/" * alpha * ")")),
      y = expression(r / alpha),
      fill = fill_expr
    ) +
    theme_linedraw() +
    theme(
      text = element_text(family = "Times New Roman", size = 20),
      axis.text = element_text(size = 18),
      axis.title = element_text(size = 20),
      plot.title = element_text(size = 20),
      legend.text = element_text(size = 18),
      legend.title = element_text(size = 20, vjust = 0.85),
      legend.position = "bottom",
      plot.background = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    guides(fill = guide_colourbar(barwidth = 31, barheight = 1))
  
  if (filename != "") {
    ggsave(filename, plot = p, width = width, height = height, units = "cm")
  }
  
  return(p)
}

# Function to show exponential values in scientific format
scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

###############################################################################
####################### Data loading and plotting #############################
###############################################################################

dfcodominant01 <- load_eq_df("../Data/codominant_sols_freqs_20250911_a=0.01.csv", 0.01)
dfcodominant01 <- load_eq_df("../Data/codominant_sols_freqs_20260925_a=0.01.csv", 0.01) # New
dfcodominant001 <- load_eq_df("../Data/codominant_sols_freqs_20250911_a=0.001.csv", 0.001)
dfcodominant001 <- load_eq_df("../Data/codominant_sols_freqs_20260925_a=0.001.csv", 0.001) #New

dfreces005 <- load_eq_df("../Data/recessive_sols_ext_1.csv", 0.005) %>% filter(r > 0)
dfreces <- load_eq_df("../Data/recessive_sols_ext_001.csv", 0.001) %>% filter(r > 0)

# Fig_1A.pdf
plot_dmi_heatmap(
  df = dfcodominant001,
  r_lim = c(1e-7, 0.49),
  rs_lim = c(1e-4, 1e2),
  max_ga = 4.0,
  fill_var = "pA",  filename = "../Fig/Fig_2A.pdf"
)

# Fig_1B.pdf
plot_dmi_heatmap(df = dfcodominant001, r_lim = c(1e-7, 0.49),
                 rs_lim = c(1e-4, 1e2),  max_ga = 4.0, fill_var = "LD", filename = "../Fig/Fig_2B.pdf"
)

# freq_recessive005.pdf
plot_dmi_heatmap(df = dfreces005, r_lim = c(1e-7, 0.49), rs_lim = c(1e-4, 1e2),
                 max_ga = 16.0, fill_var = "pA",filename = "../Fig/freq_recessive_a=005.pdf"
)

# LD_recessive_a=0.005.pdf
plot_dmi_heatmap(df = dfreces005, r_lim = c(1e-7, 0.49), rs_lim = c(1e-4, 1e2),
                 max_ga = 16.0, fill_var = "LD", filename = "../Fig/LD_recessive_a=0.005.pdf"
)



###############################################################################
####################### LD approximations #####################################
###############################################################################

a_val <- 0.001

########### Fig 2C, Freq
dfcodominant001 %>%
  filter(
    r %in% c(0, 1e-4, 1e-3, 0.5),
    ga %in%
      sort(unique(dfcodominant001$ga))[seq(
        1,
        length(unique(dfcodominant001$ga)),
        2
      )]
  ) %>%
  ggplot(aes(x = as.numeric(ga), y = pA, colour = as.factor(r))) +
  #geom_point(shape=18, size=2.5)+
  # Exact lines
  #geom_hline(yintercept = c(0.25, min(dfcodominant001$pA)), linetype='dashed')+
  geom_line(
    data = dfcodominant001 %>% filter(as.numeric(ga) <= 2),
    aes(
      x = as.numeric(ga),
      y = ((a_val) / (2 * a_val + ga * a_val)),
      linetype = "LE"
    ),
    inherit.aes = FALSE,
    linewidth = 1,
    #colour = "seagreen2"
    colour = "#45f1f6"
  ) +
  geom_line(
    data = dfcodominant001,
    aes(
      x = as.numeric(ga),
      y = ((2 * a_val) / (4 * a_val + ga * a_val)),
      linetype = "No recombination"
    ),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "#47325f"
  ) +
  geom_point(shape = 19, size = 2.5) +
  xlab(expression("Relative strength of BDMIs (" * gamma * "/" * alpha * ")")) +
  ylab(expression("Equilibrium frequency (" * P[A] * ")")) +
  scale_x_continuous(limits = c(0, 4.0), labels = function(x) {
    sprintf("%.1f", x)
  }) +
  scale_y_continuous(limits = c(0, 0.5)) +
  # Manual colour gradient
  scale_color_manual(
    values = c("#47325f", "#54a7ba", "#59ced8", "#5df5f6"),
    labels = c(
      "0",
      expression(10^-4),
      expression(10^-3),
      "0.5"
    )
  ) +
  # Automatic colour gradient
  # scale_color_gradient(
  #   trans = pseudo_log_trans(sigma = 0.00001),
  #   breaks = c(1e-5, 1e-4, 1e-3, 0.5),
  #   labels = c("0", "1e-4", "1e-3", "5e-1"),
  #   guide = guide_coloursteps(
  #     barwidth = 16,
  #     barheight = 0.8,
  #     title.vjust = 1
  #   )
  # ) +
  guides(
    colour = guide_legend(
      title = "r",
      override.aes = list(size = 4)
    )
  ) +
  scale_linetype_manual(
    name = "Exact solutions",
    values = c("LE" = "solid", "No recombination" = "solid", "lowr" = "solid"),
    labels = c("LE", "No recombination", "low r approx (r=1e-5)")
  ) +
  guides(
    linetype = guide_legend(order = 2)
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.position = c(0.98, 0.10),
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    legend.box.just = "right",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )
ggsave("../Fig/Fig_2C.pdf", width = 23, height = 15, units = "cm")



########## Fig 2D, low r, LD plot
dfcodominant001 %>%
  filter(
    r %in% c(0, 1e-6, 1e-5, 1e-4),
    ga %in% sort(unique(dfcodominant001$ga))[seq(1, length(unique(dfcodominant001$ga)), 2)]) %>%
  ggplot(aes(x = as.numeric(ga), y = LD, colour = as.factor(r))) +
  # Numerical sols (points)
  geom_point(
    aes(linetype = "Numerical", shape = "Numerical"),
    size = 2.5,
    alpha = 1
  ) +
  # Approximation (lines)
  geom_line(
    data = function(d) d %>% filter(r %in% c(0, 1e-6, 1e-5, 1e-4)),
    aes(
      y = (-(4 * a_val^2) / (4 * a_val + ga * a_val)^2) +
        ((8 * a_val + 3 * ga * a_val) * 8 * r * (a_val^2)) /
        (ga * a_val * (4 * a_val + ga * a_val)^3),
      group = r,
      linetype = "Analytical",
      shape = "Analytical"
    ),
    linewidth = 0.75
  ) +
  # LE line
  # geom_line(
  #   data = function(d) d %>% filter(r %in% c(0.001, 0.5)),
  #   aes(
  #     y = -((a_val^2) * ga * a_val * (a_val + ga * a_val)) /
  #       (2 *
  #         (2 * a_val + ga * a_val) *
  #         (a_val * ((ga * a_val)^3) + r * (2 * a_val + ga * a_val)^3)),
  #     group = r,
  #     linetype = "Approximation",
  #     shape = "Approximation"
  #   ),
  #   linewidth = 0.75
  # ) +
  xlab(expression("Relative strength of BDMIs (" * gamma * "/" * alpha * ")")) +
  ylab(expression("Linkage disequilibrium (D)")) +
  scale_x_continuous(
    limits = c(0, 4.0),
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_continuous(
    limits = c(-0.25, 0)
  ) +
  scale_colour_brewer(
    palette = "Dark2",
    name = expression(r),
    labels = function(x) {
      parse(text = formatC(as.numeric(x), format = "e", digits = 1 ))
    }
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("Numerical" = "blank", "Analytical" = "solid"),
    guide = "legend"
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("Numerical" = 16, "Analytical" = NA)
  ) +
  scale_color_manual(
    name = "r",
    values = c("#47325f", "#4b597d", "#50809b", "#54a7ba", "#59ced8", "#5df5f6"),
    labels = c( "0", expression(10^-6), expression(10^-5), expression(10^-4))
  ) +
  guides(
    color = guide_legend(order = 2),
    linetype = guide_legend(order = 1),
    shape = guide_legend(order = 1)
  ) +
  
  theme_linedraw() +
  theme(
    text = element_text(family = "Times New Roman", size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.position = c(0.98, 0.10),
    legend.justification = c("right", "bottom"),
    legend.direction = "horizontal",
    legend.box.just = "right",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

ggsave("../Fig/Fig_2D.pdf", width = 23, height = 15, units = "cm")



################################################################################
################## Plotting mean fitness (unused) #############################
################################################################################

w_bar <- function(x1, x2, x3, x4, alpha, beta, gamma1, gamma2, gamma3, gamma4) {
  term1 <- 2 * x1 * (beta * x2 + alpha * x3 + x4 * (alpha + beta - gamma1))
  term2 <- 2 * x2 * (x3 * (alpha + beta - gamma1) + x4 * (alpha - gamma3))
  term3 <- -x4 * (-2 * beta * x3 + 2 * gamma2 * x3 + gamma4 * x4)
  
  term1 + term2 + term3
}

w_bar(0.5,0.5,0,0,0.001,0.001,0.001,0.002,0.002,0.004)


allele_to_haplo <- function(pA, pB, LD) {
  x1 <- 1 + LD - pA - pB + pA * pB
  x2 <- -LD + pB - pA * pB
  x3 <- -LD + pA - pA * pB
  x4 <- 1 - x1 - x2 - x3  # implied 4th haplotype frequency
  
  c(x1, x2, x3, x4)
}

scientific_superscript <- function(x) {
  exponent <- floor(log10(x))
  mantissa <- x / 10^exponent
  as.expression(mapply(function(m, e) {
    if (round(m, 6) == 1) {
      bquote(10^.(e))
    } else {
      bquote(.(m) %*% 10^.(e))
    }
  }, mantissa, exponent))
}


alpha <- 0.01
beta <- 0.01
gamma <- 0.02
gamma1 <- gamma/2
gamma2 <- gamma
gamma3 <- gamma
gamma4 <- gamma*2
x1 <- 0.2
x2 <- 0.4
x3 <- 0.4
x4 <- 0

alpha=0.001
# Modify df to have wbar
dfcodominant001 <- dfcodominant001 |> 
  mutate(
    pB = pA,
    x1 = 1 + LD - pA - pB + pA * pB,
    x2 = -LD + pB - pA * pB,
    x3 = -LD + pA - pA * pB,
    x4 = 1 - x1 - x2 - x3,
    wbar = w_bar(x1, x2, x3, x4, alpha, alpha, (ga*alpha)/2, (ga*alpha), (ga*alpha), (ga*alpha)*2)
  )


dfcodominant001 %>%
  filter(
    r <= 0.49 & m == 0 & ga <= 8 & r >= 0.00000001 & rs >= 0.0001 & rs <= 100
  ) %>%
  ggplot(aes(x = ga, y = rs, fill = wbar)) +
  geom_tile() +
  scale_fill_gradientn(
    # expression(P[A]),
    colors = tail(hcl.colors(100, "Heat", rev = T), 95),
    limits = c(0.0004992, 0.001),
    breaks = c(5e-4, 7.5e-4, 1e-3),
    labels = scientific_superscript
  ) +
  geom_contour(
    aes(z = wbar),
    colour = "#FFF",
    binwidth = 0.00002,
    linetype = "dashed",
    alpha = 0.6
  ) +
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
    x = expression(paste(
      "Relative strength of BDMIs (" * gamma * "/" * alpha * ")"
    )),
    y = expression(r / alpha),
    fill = expression(bar(w))
  ) +
  theme_linedraw() +
  theme(
    text = element_text(family = "Times New Roman", size = 20),
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    plot.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20, vjust = 0.85),
    legend.position = "bottom",
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  ) +
  guides(fill = guide_colourbar(barwidth = 31, barheight = 1))

