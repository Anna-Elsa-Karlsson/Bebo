## models_results.R  -- trial-level models (analysis plan section 2)
## RESULTS step only: fits the models, writes the statistics, saves the fitted
## models + simple-slopes objects. Figures: models_plots.R (no computation).
##
## Input  (from g3_extract.m):  <results>/tables/items.csv, pairs.csv
## Output: <results>/stats/models_stats.txt, <results>/stats/models.rds
##
## Coding (fixed in the plan):
##   AgeNum   YA +0.5, OA -0.5  (numeric; NEVER a factor with default
##            contrasts, which gives dummy coding with YA as reference)
##   catSize  4 exemplars +0.5, 2 exemplars -0.5
##   itemRep  4 presentations +0.5, 2 presentations -0.5
##   Power1st, PowerPair, ItemSimilarity, CatSimilarity: z-scored WITHIN person
## Set the results folder with the environment variable BEBO_RESULTS, or
## edit the default below.

suppressPackageStartupMessages({
  library(lme4); library(lmerTest); library(interactions); library(MuMIn); library(dplyr)
})

results_dir <- Sys.getenv("BEBO_RESULTS", "/Users/anna.karlsson/Documents/CL/CL/BEBO/Results_v2")
tab_dir  <- file.path(results_dir, "tables")
stat_dir <- file.path(results_dir, "stats")
dir.create(stat_dir, showWarnings = FALSE, recursive = TRUE)

ctrl_g <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
ctrl_l <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
zw <- function(x) as.numeric(scale(x))     # used inside group_by(ID): within-person z
muffle_r2 <- function(expr) withCallingHandlers(expr, warning = function(w)
  if (grepl("r-squared", conditionMessage(w))) invokeRestart("muffleWarning"))

sink(file.path(stat_dir, "models_stats.txt"), split = TRUE)
cat("BEBO trial-level models --", format(Sys.time()), "\nResults folder:", results_dir, "\n\n")

## ============================================================
## DATA
## ============================================================
items <- read.csv(file.path(tab_dir, "items.csv")) %>%
  mutate(ID = factor(ID), category = factor(category), item = factor(item),
         AgeNum = ifelse(YA == 1, 0.5, -0.5),
         catSizeC = ifelse(catSize == 4, 0.5, -0.5),
         itemRepC = ifelse(itemRep == 4, 0.5, -0.5)) %>%
  group_by(ID) %>%
  mutate(Power1st = zw(Power1st),
         ItemSimilarity = zw(Isim),
         CatSimilarity  = zw(Csim)) %>%   # scale() ignores NA: z over the items that have a Csim value
  ungroup() %>% as.data.frame()

pairs <- read.csv(file.path(tab_dir, "pairs.csv")) %>%
  mutate(ID = factor(ID), category = factor(category),
         AgeNum = ifelse(YA == 1, 0.5, -0.5),
         TypeNum = ifelse(Type == "Isim", 0.5, -0.5),
         catSizeC = ifelse(catSize == 4, 0.5, -0.5),
         itemRepC = ifelse(itemRep == 4, 0.5, -0.5)) %>%
  group_by(ID) %>% mutate(PowerPair = zw(PowerPair)) %>% ungroup() %>% as.data.frame()

cat("Items:", nrow(items), "| with a Csim value:", sum(!is.na(items$CatSimilarity)),
    "| persons:", nlevels(items$ID), "\n")
cat("Pairs:", nrow(pairs), "(Isim", sum(pairs$Type == "Isim"), ", Csim", sum(pairs$Type == "Csim"), ")\n\n")

## ============================================================
## 2a. CONDITION MODELS (descriptive check; do not decide 2b/2c)
## ============================================================
cat("==================== 2a. CONDITION MODELS ====================\n")
cond <- list()
for (dv in c("Power1st", "Isim", "Csim")) {
  cat("\n---", dv, "~ catSize * itemRep * Age ---\n")
  f <- as.formula(paste(dv, "~ catSizeC * itemRepC * AgeNum + (1|ID) + (1|category)"))
  cond[[dv]] <- lmer(f, data = items, control = ctrl_l)
  print(summary(cond[[dv]])$coefficients, digits = 3)
}

## ============================================================
## 2b. PAIR MODEL (Q4)
## ============================================================
cat("\n==================== 2b. PAIR MODEL ====================\n")
m_pair <- lmer(Sim ~ PowerPair * TypeNum * AgeNum + itemRepC + catSizeC + (1|ID) + (1|category),
               data = pairs, control = ctrl_l)
print(summary(m_pair), correlation = FALSE)
if (summary(m_pair)$coefficients["PowerPair:TypeNum", "Pr(>|t|)"] < .05) {
  cat("\n--- Follow-up: power slope per pair type ---\n")
  print(muffle_r2(sim_slopes(m_pair, pred = PowerPair, modx = TypeNum,
        modx.values = c(-0.5, 0.5), modx.labels = c("Csim", "Isim"))))
}

## ============================================================
## 2c. MEMORY MODEL (Q3)
## ============================================================
cat("\n==================== 2c. MEMORY MODEL ====================\n")
dat_mem <- items[!is.na(items$CatSimilarity), ]
m_mem <- glmer(Acc ~ Power1st * ItemSimilarity * CatSimilarity * AgeNum + catSizeC + itemRepC +
                 (1|ID) + (1 + AgeNum || category:item),
               data = dat_mem, family = binomial, control = ctrl_g)
print(summary(m_mem), correlation = FALSE)
cat("\nSingular fit:", isSingular(m_mem), "\n")
cat("R2 (Nakagawa, theoretical) marginal / conditional:\n")
print(round(MuMIn::r.squaredGLMM(m_mem)["theoretical", ], 3))

cat("\n--- Random-effects check: without the by-item age slope ---\n")
m_mem_noSlope <- glmer(Acc ~ Power1st * ItemSimilarity * CatSimilarity * AgeNum + catSizeC + itemRepC +
                         (1|ID) + (1|category:item),
                       data = dat_mem, family = binomial, control = ctrl_g)
cf1 <- summary(m_mem)$coefficients; cf0 <- summary(m_mem_noSlope)$coefficients
print(data.frame(term = rownames(cf1), est = round(cf1[,1],3), p = round(cf1[,4],4),
                 est_noSlope = round(cf0[,1],3), p_noSlope = round(cf0[,4],4), row.names = NULL))

## Follow-ups: ONLY for significant two-way interactions (plan 2c)
cat("\n--- Follow-ups (simple slopes at -1 SD / mean / +1 SD; Johnson-Neyman with FDR; logit) ---\n")
p_of <- function(term) { r <- rownames(cf1); i <- which(r == term)
  if (!length(i)) i <- which(r == paste(rev(strsplit(term, ":")[[1]]), collapse = ":"))
  if (length(i)) cf1[i, 4] else NA }
follow <- list(
  list(term = "Power1st:ItemSimilarity",      pred = "ItemSimilarity", modx = "Power1st"),
  list(term = "Power1st:CatSimilarity",       pred = "CatSimilarity",  modx = "Power1st"),
  list(term = "ItemSimilarity:CatSimilarity", pred = "ItemSimilarity", modx = "CatSimilarity"),
  list(term = "Power1st:AgeNum",              pred = "Power1st",       modx = "AgeNum"),
  list(term = "ItemSimilarity:AgeNum",        pred = "ItemSimilarity", modx = "AgeNum"),
  list(term = "CatSimilarity:AgeNum",         pred = "CatSimilarity",  modx = "AgeNum"))
ss <- list()
for (f in follow) {
  p <- p_of(f$term)
  cat(sprintf("\n%s: p = %.4f%s\n", f$term, p, ifelse(p < .05, "", " -> no follow-up")))
  if (is.na(p) || p >= .05) next
  args <- list(model = m_mem, pred = f$pred, modx = f$modx, jnplot = FALSE)
  if (f$modx == "AgeNum") {
    args$modx.values <- c(-0.5, 0.5); args$modx.labels <- c("Older", "Younger"); args$johnson_neyman <- FALSE
  } else { args$johnson_neyman <- TRUE; args$control.fdr <- TRUE }
  ss[[f$term]] <- muffle_r2(do.call(sim_slopes, args))
  print(ss[[f$term]])
}

## ============================================================
## SAVE
## ============================================================
saveRDS(list(cond = cond, m_pair = m_pair, m_mem = m_mem, m_mem_noSlope = m_mem_noSlope,
             ss = ss, dat_mem = dat_mem, pairs = pairs),
        file.path(stat_dir, "models.rds"))
cat("\nSaved models.rds\n")
sink()
