library(tidyverse)
library(dplyr)
library(minpack.lm)
library(ggplot2)
library(scico) #Testing diff colour palettes
library(extrafont)
loadfonts(device = "pdf", quiet = TRUE)

df <- alldfs %>%
  filter(pA0==0.7,a == 0, allele=="N", g>0, r_bc>0, type=="codom")

# Deterministic sweep
sweep_det <- function(r, H0, s) {
  # s here acts as the combined compound selection/recombination scale parameter (alpha * r / s)
  #H0 * (1 - exp(-r*(1 - s)/s)) # 
  ### From Charlesworth eq 8.5
  H0 * 1/(r+s)
  #H0 * (1 - exp(-r / s)) #If r close to 0, this should be the case?
}
# Deterministic bgs expectation
bgs_det <- function(r, H0, s) {
  H0 * exp(- s / r)
}

#### Wiehe and Stephen 1993 recurrent hitchhiking model; 
#### pi = theta*rho/(rho+k*alpha*v)
### Where k=0.075, v is the rate of adaptive substitution by site, and alpha is the strength of selection, 
### and theta=4*Ne*mu
rrh_det <- function(r, H0, s) {
  H0 * r/(r+s)
}

fit_results_bgs <- df %>%
  filter(r_bc > 0) %>%
  group_by(g) %>%
  do({
    sub <- .
    fit <- tryCatch(
      nlsLM(
        het ~ bgs_det(r_bc, H0, s),
        data = sub,
        start  = list(H0 = max(sub$het), s = unique(sub$g)),
        lower  = c(H0 = 0, s = 1e-10),
        upper  = c(H0 = 1, s = 1),
        control = nls.lm.control(maxiter = 500)
      ),
      error = function(e) { cat("ERROR g =", unique(sub$g), ":", e$message, "\n"); NULL }
    )
    if (!is.null(fit)) {
      data.frame(
        g      = unique(sub$g),
        r_ab  = unique(sub$r_ab),
        H0_fit = coef(fit)["H0"],
        s  = coef(fit)["s"],
        ratio  = coef(fit)["s"] / unique(sub$g),
        converged = TRUE
      )
    } else {
      data.frame(g = unique(sub$g), r_ab = unique(sub$r_ab),, H0_fit = NA, s = NA, ratio = NA, converged = FALSE)
    }
  }) %>%
  ungroup()

fit_results_sweep <- df %>%
  filter(r_bc > 0) %>%
  group_by(g, r_ab) %>%
  do({
    sub <- .
    fit <- tryCatch(
      nlsLM(
        het ~ sweep_det(r_bc, H0, s),
        data = sub,
        start  = list(H0 = max(sub$het), s = unique(sub$g)),
        lower  = c(H0 = 0, s = 1e-10),
        upper  = c(H0 = 1, s = 1),
        control = nls.lm.control(maxiter = 500)
      ),
      error = function(e) { cat("ERROR g =", unique(sub$g), ":", e$message, "\n"); NULL }
    )
    if (!is.null(fit)) {
      data.frame(
        g = unique(sub$g),
        r_ab  = unique(sub$r_ab),
        H0_fit = coef(fit)["H0"],
        s  = coef(fit)["s"],
        ratio  = coef(fit)["s"] / unique(sub$g),
        converged = TRUE
      )
    } else {
      data.frame(g = unique(sub$g), r_ab = unique(sub$r_ab),, H0_fit = NA, s = NA, ratio = NA, converged = FALSE)
    }
  }) %>%
  ungroup()


fit_results_rrh <- df %>%
  filter(r_bc > 0) %>%
  group_by(g, r_ab) %>%
  do({
    sub <- .
    fit <- tryCatch(
      nlsLM(
        het ~ rrh_det(r_bc, H0, s),
        data = sub,
        start  = list(H0 = max(sub$het), s = unique(sub$g)),
        lower  = c(H0 = 0, s = 1e-10),
        upper  = c(H0 = 1, s = 1),
        control = nls.lm.control(maxiter = 500)
      ),
      error = function(e) { cat("ERROR g =", unique(sub$g), ":", e$message, "\n"); NULL }
    )
    if (!is.null(fit)) {
      data.frame(
        g = unique(sub$g),
        r_ab  = unique(sub$r_ab),
        H0_fit = coef(fit)["H0"],
        s  = coef(fit)["s"],
        ratio  = coef(fit)["s"] / unique(sub$g),
        converged = TRUE
      )
    } else {
      data.frame(g = unique(sub$g), r_ab = unique(sub$r_ab),, H0_fit = NA, s = NA, ratio = NA, converged = FALSE)
    }
  }) %>%
  ungroup()

print(fit_results_sweep)
print(fit_results_bgs)
print(fit_results_rrh)


# Generate predicted curves
pred_df_bgs <- fit_results_bgs %>%
  filter(converged) %>%
  rowwise() %>%
  do({
    row <- .
    r_seq <- 10^seq(log10(min(df$r_bc[df$r_bc > 0])),log10(max(df$r_bc)),length.out = 300)
    data.frame(
      g = row$g,
      r_ab = row$r_ab,
      r_bc = r_seq,
      het_pred = bgs_det(r_seq, row$H0_fit, row$s)
    )
  })

pred_df_sweep <- fit_results_sweep %>%
  filter(converged) %>%
  rowwise() %>%
  do({
    row <- .
    r_seq <- 10^seq(log10(min(df$r_bc[df$r_bc > 0])),log10(max(df$r_bc)),length.out = 300)
    data.frame(
      g = row$g,
      r_ab = row$r_ab,
      r_bc = r_seq,
      het_pred = sweep_det(r_seq, row$H0_fit, row$s)
    )
  })

pred_df_rrh <- fit_results_rrh %>%
  filter(converged) %>%
  rowwise() %>%
  do({
    row <- .
    r_seq <- 10^seq(log10(min(df$r_bc[df$r_bc > 0])),log10(max(df$r_bc)),length.out = 300)
    data.frame(
      g = row$g,
      r_ab = row$r_ab,
      r_bc = r_seq,
      het_pred = rrh_det(r_seq, row$H0_fit, row$s),
      s = as.factor(row$s)
    )
  })


# Plot
rhh_plot <- df %>% filter(g > 0.005, r_ab >= 1e-2 | r_ab <= 1e-6) %>%
  ggplot() +
  geom_point(aes(x = r_bc, y = het, color = factor(g)), size = 2) +
  #geom_line(data = pred_df_bgs, aes(x = r_bc, y = het_pred, color = factor(g))) +
  #geom_line(data = pred_df_sweep, aes(x = r_bc, y = het_pred, color = factor(g)), linetype = "dashed") +
  geom_line(data = pred_df_rrh %>% filter(g > 0.005, r_ab >= 1e-2 | r_ab <= 1e-6), aes(x = r_bc, y = het_pred, color = factor(g), linetype = "")) +
  facet_wrap(.~r_ab, ncol = 1, strip.position = "right", labeller = label_bquote(r[ab] == 10^.(log10(r_ab)))) +
  scale_x_log10() +
  scale_colour_discrete(palette = grDevices::colorRampPalette(c("dodgerblue", "firebrick2"))(6)) +
  scale_linetype_manual(
    name = "Recurrent hitchhiking ",
    values = c("44", "solid", "dotted"),
    labels = c("", "", "")
  ) +
  labs(
    x = expression(r[bn]), 
    y = "Heterozygosity", 
    color = expression("BDMI " ~ (gamma))
  ) +
  guides(
    colour = guide_legend(nrow = 1, order = 2),
    linetype = guide_legend(order = 1)
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text.y = element_text(size = 10),
    legend.key.height = unit(0.01, 'cm'),
    text=element_text(family="Times New Roman"))

rhh_plot
ggsave("../Fig/Fig_s_rhh_A.pdf", plot = rhh_plot, 
       width = 23, height = 20, units = "cm", dpi=800)

rhh_plot_scoeff <-fit_results_rrh %>% 
  filter(g>0.005, r_ab>0) %>%
  ggplot(aes(x=r_ab, y=s, colour = as.factor(g))) +
  geom_point(size = 2) +
  geom_line() +
  #geom_smooth(se = 1) +
  labs(
    x = expression(r[ab]), 
    y = "Inferred selection coefficient (s)", 
    color = expression(gamma)
  ) +
  scale_x_continuous(trans='log10') +
  scale_colour_discrete(palette=grDevices::colorRampPalette(c("dodgerblue", "firebrick2"))(6))+
  theme_bw() +
  theme(legend.position = "bottom",
        text=element_text(family="Times New Roman")) +
  guides(colour = guide_legend(nrow = 1, order = 2))
rhh_plot_scoeff
ggsave("../Fig/Fig_s_rhh_B.pdf", plot = rhh_plot_scoeff, 
       width = 23, height = 20, units = "cm", dpi=800)

rhh_plot_heatmap <- fit_results_rrh %>% 
  filter(g>0.005, r_ab>0) %>%
  ggplot(aes(x=r_ab, y=g, fill = s)) +
  geom_point(size = 2) +
  geom_tile() +
  #geom_smooth(se = 1) +
  labs(
    x = expression(r[ab]), 
    fill = "Inferred selection coefficient (s)", 
    y = expression(gamma)
  ) +
  scale_fill_scico(palette = 'lajolla', direction = -1) +
  scale_x_continuous(trans='log10') +
  theme_bw() +
  theme(legend.position = "bottom",
        text=element_text(family="Times New Roman")) +
  guides(colour = guide_legend(nrow = 1, order = 2))

rhh_plot_heatmap
ggsave("../Fig/Fig_s_rhh_C.pdf", plot = rhh_plot_heatmap,
       width = 23, height = 20, units = "cm", dpi=800)
