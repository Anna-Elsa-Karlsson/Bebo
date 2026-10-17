## models_plots.R  -- figures for the trial-level models (PLOT step only)
## Loads models.rds from models_results.R and draws; no statistics here.
## Style as BEBO_behavior_glmer.R / BEBO_joint_glmer_plots.R: interact_plot
## for the age panels, plot(sim_slopes) for the interaction panels.
##
## Memory model figure (Fig5_memory_model):
##   A  Item similarity x age group
##   B  Alpha/beta power x age group
##   C  Item similarity x alpha/beta power
##   D  Item similarity x category similarity
##   E  Category similarity x alpha/beta power
## C-E use the saved simple-slopes object when the interaction was
## significant (follow-up computed in models_results.R), otherwise
## interact_plot (marked "n.s." in the panel title).
## Pair model figure (Fig_pair_model): similarity x power, per pair type.

suppressPackageStartupMessages({
  library(lme4); library(ggplot2); library(interactions); library(patchwork)
})

results_dir <- Sys.getenv("BEBO_RESULTS", "/Users/anna.karlsson/Documents/CL/CL/BEBO/Results_v2")
R <- readRDS(file.path(results_dir, "stats", "models.rds"))
fig_dir <- file.path(results_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

col_YA <- "#1E90FF"; col_OA <- "#505050"
theme_bebo <- theme_classic() +
  theme(text = element_text(size = 20), axis.text = element_text(size = 18),
        axis.title = element_text(size = 20), legend.text = element_text(size = 18),
        legend.title = element_text(size = 18), plot.title = element_text(size = 20, hjust = 0.5),
        plot.tag = element_text(size = 24, face = "bold"))
save_plot <- function(p, name, w, h) {
  for (ext in c("pdf", "png")) ggsave(file.path(fig_dir, paste0(name, ".", ext)), p, width = w, height = h, dpi = 300)
  cat("Saved:", name, "\n")
}

m_mem   <- R$m_mem
dat_mem <- R$dat_mem          # interact_plot looks up the model data by this name
pairs   <- R$pairs
age_vals <- c(-0.5, 0.5); age_labels <- c("Older", "Younger"); age_cols <- c(col_OA, col_YA)
y_lab <- "P('old')"

age_panel <- function(pred, xlab) {
  interact_plot(m_mem, pred = !!pred, modx = AgeNum, modx.values = age_vals, modx.labels = age_labels,
                interval = TRUE, int.type = "confidence", colors = age_cols, line.thickness = 0.6,
                x.label = xlab, y.label = y_lab, legend.main = "Age group") + theme_bebo
}
int_panel <- function(term, pred, modx, xlab, legend) {
  if (!is.null(R$ss[[term]])) {
    # plot(sim_slopes) shows the SLOPE of pred (logit) at each level of the moderator
    plot(R$ss[[term]]) + labs(x = paste0("Slope of ", tolower(sub(" \\(z\\)", "", xlab)), " (logit)"), y = legend) +
      theme_bebo
  } else {
    interact_plot(m_mem, pred = !!pred, modx = !!modx, interval = TRUE, int.type = "confidence",
                  x.label = xlab, y.label = y_lab, legend.main = legend) +
      ggtitle("(interaction n.s.)") + theme_bebo
  }
}

pA <- age_panel("ItemSimilarity", "Item similarity (z)")
pB <- age_panel("Power1st", "Alpha/beta power (z)")
pC <- int_panel("Power1st:ItemSimilarity", "ItemSimilarity", "Power1st", "Item similarity (z)", "Alpha/beta power")
pD <- int_panel("ItemSimilarity:CatSimilarity", "ItemSimilarity", "CatSimilarity", "Item similarity (z)", "Category similarity")
pE <- int_panel("Power1st:CatSimilarity", "CatSimilarity", "Power1st", "Category similarity (z)", "Alpha/beta power")
fig <- (pA | pB | pC) / (pD | pE | plot_spacer()) + plot_annotation(tag_levels = "A")
save_plot(fig, "Fig5_memory_model", 24, 13)

m_pair <- R$m_pair
pP <- interact_plot(m_pair, pred = PowerPair, modx = TypeNum, mod2 = AgeNum,
                    modx.values = c(-0.5, 0.5), modx.labels = c("Category similarity", "Item similarity"),
                    mod2.values = age_vals, mod2.labels = age_labels,
                    interval = TRUE, int.type = "confidence",
                    x.label = "Alpha/beta power of the pair (z)", y.label = "Similarity (Fisher z, corrected)",
                    legend.main = "Pair type") + theme_bebo
save_plot(pP, "Fig_pair_model", 14, 7)
