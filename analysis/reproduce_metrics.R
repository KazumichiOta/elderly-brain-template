#!/usr/bin/env Rscript
# Recompute numerical results only (no plots) from released CSV files.
# License: CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)
# You are free to share/adapt with attribution and ShareAlike.
#
# Inputs (repo root):
#   - IXI_region_dice.csv
#   - IXI_metrics(SSIM_CC_MSE).csv   [optional, if present]
#
# Outputs:
#   - outputs/mixed_effects_beta3_by_region.csv      # subcortical β3, 95%CI, p, FDR, crossing age
#   - outputs/subcortical_trend_slopes_tukey.csv     # Age slopes by template + Tukey-adjusted p
#   - outputs/tissue_interaction_and_trends.csv      # WM/GM/CSF β2, β3, crossing age, BH, Tukey
#   - outputs/metrics_interaction_and_trends.csv     # 3D-SSIM/CC/MSE interaction coef + CI + BH
#   - outputs/long_IXI_region_dice.csv               # (wide->long) normalized Dice table
#   - analysis/sessionInfo.txt                       # session info (always written here)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(lme4);  library(lmerTest); library(emmeans)
})

# -------------------- helpers --------------------
norm_key <- function(x) tolower(gsub("[^a-z0-9]", "", x))
find1 <- function(cands, want){
  nk <- norm_key(cands)
  for (w in want) {
    hit <- which(nk == norm_key(w))
    if (length(hit)) return(cands[hit[1]])
  }
  NA_character_
}
get_script_dir <- function(){
  args <- commandArgs(trailingOnly = FALSE)
  i <- grep("^--file=", args)
  if (length(i)) return(dirname(normalizePath(sub("^--file=", "", args[i[1]]))))
  if (!is.null(sys.frames()[[1]]$ofile)) return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  normalizePath(getwd())
}
# Age_cross = -b2 / b3, with delta-method CI
age_cross_delta <- function(beta_temp, beta3, var_tt, cov_t3, var_33){
  if (is.na(beta_temp) || is.na(beta3) || abs(beta3) < 1e-12) return(c(NA, NA, NA))
  est <- - beta_temp / beta3
  g1 <- -1 / beta3
  g2 <-  beta_temp / (beta3^2)
  var <- g1^2 * var_tt + 2*g1*g2 * cov_t3 + g2^2 * var_33
  se  <- sqrt(max(var, 0))
  c(est, est - 1.96*se, est + 1.96*se)
}

# -------------------- paths --------------------
script_dir <- get_script_dir()
root_dir   <- if (basename(script_dir) == "analysis") dirname(script_dir) else script_dir
out_dir    <- file.path(root_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ensure analysis/ exists for sessionInfo.txt
dir.create(file.path(root_dir, "analysis"), showWarnings = FALSE, recursive = TRUE)

# -------------------- load Dice CSV --------------------
csv_dice <- file.path(root_dir, "IXI_region_dice.csv")
if (!file.exists(csv_dice)) stop("CSV not found: ", csv_dice)
raw <- read_csv(csv_dice, show_col_types = FALSE)
nm  <- names(raw)

col_id       <- find1(nm, c("ID","Subject","SubjectID","subject_id"))
col_age      <- find1(nm, c("Age","AgeYears","Age (years)","age_years","age"))
col_template <- find1(nm, c("template","Template"))
if (any(is.na(c(col_id,col_age,col_template))))
  stop("Core columns not found. Available: ", paste(nm, collapse=", "))

# LONG or WIDE?
col_dice   <- find1(nm, c("Dice","DiceCoefficient","DiceValue","dice"))
col_region <- find1(nm, c("Region","ROI","Label","Structure","region"))
wide_idx   <- which(startsWith(tolower(nm), "dice"))  # e.g., Dice_brainstem, etc.

if (!is.na(col_dice) && !is.na(col_region)) {
  dice <- raw %>%
    rename(ID=!!col_id, Age=!!col_age, template=!!col_template,
           Region=!!col_region, Dice=!!col_dice)
} else if (length(wide_idx) > 0) {
  dice_cols <- nm[wide_idx]
  dice <- raw %>%
    rename(ID=!!col_id, Age=!!col_age, template=!!col_template) %>%
    pivot_longer(all_of(dice_cols), names_to="Region_raw", values_to="Dice") %>%
    mutate(Region = sub("^[Dd]ice[_-]*", "", Region_raw))
} else {
  stop("No Dice column(s) detected.")
}

# normalize values and labels
map_region <- c("amygdala"="amygdala","brainstem"="brainstem","caudate"="caudate",
                "hippocampus"="hippocampus","pallidum"="pallidum",
                "putamen"="putamen","thalamus"="thalamus",
                "wm"="WM","gm"="GM","csf"="CSF")
dice <- dice %>%
  mutate(
    ID   = factor(ID),
    Age  = suppressWarnings(as.numeric(Age)),
    Dice = suppressWarnings(as.numeric(Dice)),
    template = case_when(
      grepl("^mni",     tolower(template)) ~ "MNI152",
      grepl("^elderly", tolower(template)) ~ "Elderly",
      TRUE ~ as.character(template)
    ),
    template = factor(template, levels=c("MNI152","Elderly")),
    Region = tolower(Region),
    Region = ifelse(Region %in% names(map_region), unname(map_region[Region]), Region)
  ) %>%
  filter(!is.na(Age), !is.na(Dice), !is.na(template), !is.na(Region))

# save normalized long table
write_csv(dice, file.path(out_dir, "long_IXI_region_dice.csv"))

# Box–Cox (λ = 1.22) for Dice (as in manuscript)
lambda <- 1.22; eps <- 1e-6
bc <- function(x, lam=lambda) if (lam == 0) log(x + eps) else ((x + eps)^lam - 1) / lam
dice <- dice %>% mutate(Dice_bc = bc(Dice))

# -------------------- model extractor --------------------
extract_numbers <- function(df, label){
  fit <- lmer(Dice_bc ~ Age * template + (1|ID), data = df, REML=FALSE)
  s   <- summary(fit); co <- coef(s); rn <- rownames(co)

  r_int <- rn[grepl("Age:templateElderly|templateElderly:Age", rn)]
  if (length(r_int) != 1) r_int <- "Age:templateElderly"

  b_age <- co["Age","Estimate"]
  b_tmp <- if ("templateElderly" %in% rn) co["templateElderly","Estimate"] else NA_real_
  b_int <- if (r_int %in% rn) co[r_int,"Estimate"] else NA_real_
  se_i  <- if (r_int %in% rn) co[r_int,"Std. Error"] else NA_real_
  p_i   <- if (r_int %in% rn) co[r_int,"Pr(>|t|)"] else NA_real_
  ci_i  <- c(b_int - 1.96*se_i, b_int + 1.96*se_i)

  V <- as.matrix(vcov(fit))
  idx_t <- which(names(fixef(fit)) == "templateElderly")
  idx_i <- which(names(fixef(fit)) %in% c("Age:templateElderly","templateElderly:Age"))
  if (length(idx_t) && length(idx_i)){
    ac <- age_cross_delta(b_tmp, b_int, V[idx_t,idx_t], V[idx_t,idx_i], V[idx_i,idx_i])
  } else ac <- c(NA, NA, NA)

  et  <- suppressMessages(emtrends(fit, ~ template, var="Age"))
  est <- as.data.frame(summary(et))
  cmp <- as.data.frame(pairs(et, adjust="tukey"))
  elderly_slope <- est$Age.trend[est$template=="Elderly"]
  mni_slope     <- est$Age.trend[est$template=="MNI152"]
  delta_slope   <- if (nrow(cmp)) cmp$estimate[1] else NA_real_
  p_tukey       <- if (nrow(cmp)) cmp$p.value[1]  else NA_real_

  tibble(
    Label = label,
    Beta2 = b_age, Beta3 = b_int,
    InteractionCoef = b_int,
    Coef_CI_lower = ci_i[1], Coef_CI_upper = ci_i[2],
    p_value_raw = p_i,
    AgeCross_est = ac[1], AgeCross_Lower = ac[2], AgeCross_Upper = ac[3],
    Elderly_slope = elderly_slope, MNI_slope = mni_slope,
    Delta_slope = delta_slope, p_Tukey = p_tukey
  )
}

# -------------------- Subcortical (7 ROI) --------------------
subs <- c("amygdala","brainstem","caudate","hippocampus","pallidum","putamen","thalamus")
df_sub <- dice %>% filter(tolower(Region) %in% subs)
if (nrow(df_sub) > 0){
  num_sub <- df_sub %>%
    group_by(Region) %>% group_split() %>%
    lapply(function(d) extract_numbers(d, unique(d$Region))) %>%
    bind_rows() %>%
    mutate(p_FDR = p.adjust(p_value_raw, method="BH"))

  write_csv(num_sub %>% select(
    Region=Label, InteractionCoef, Coef_CI_lower, Coef_CI_upper,
    p_value_raw, AgeCross_est, AgeCross_Lower, AgeCross_Upper, p_FDR
  ), file.path(out_dir, "mixed_effects_beta3_by_region.csv"))

  write_csv(num_sub %>% select(
    Region=Label, Elderly_slope, MNI_slope, Delta_slope,
    `adj_p_Tukey` = p_Tukey
  ), file.path(out_dir, "subcortical_trend_slopes_tukey.csv"))
}

# -------------------- Tissue (WM / GM / CSF) --------------------
df_tissue <- dice %>% filter(Region %in% c("WM","GM","CSF"))
if (nrow(df_tissue) > 0){
  num_tissue <- df_tissue %>%
    group_by(Region) %>% group_split() %>%
    lapply(function(d) extract_numbers(d, unique(d$Region))) %>%
    bind_rows() %>%
    mutate(p_bh = p.adjust(p_value_raw, method="BH"))

  write_csv(num_tissue %>% transmute(
    Tissue = Label,
    Beta2, Beta3,
    CrossAge_est = AgeCross_est,
    CI_lower = AgeCross_Lower, CI_upper = AgeCross_Upper,
    p_interaction = p_value_raw, p_bh, p_Tukey
  ), file.path(out_dir, "tissue_interaction_and_trends.csv"))
}

# -------------------- Whole-brain metrics (optional) --------------------
csv_metrics <- list.files(root_dir, pattern="^IXI_metrics.*\\.csv$", full.names=TRUE)
if (length(csv_metrics) >= 1 && file.exists(csv_metrics[1])){
  wb <- read_csv(csv_metrics[1], show_col_types = FALSE)
  n2 <- names(wb)
  id2  <- find1(n2, c("ID","Subject","SubjectID","subject_id"))
  age2 <- find1(n2, c("Age","AgeYears","Age (years)","age_years","age"))
  tpl2 <- find1(n2, c("template","Template"))
  mtr  <- find1(n2, c("Metric","metric"))
  val  <- find1(n2, c("Value","value","Score","score"))
  if (!any(is.na(c(id2,age2,tpl2,mtr,val)))) {
    wb <- wb %>%
      rename(ID=!!id2, Age=!!age2, template=!!tpl2, Metric=!!mtr, Value=!!val) %>%
      mutate(
        template = case_when(
          grepl("^mni", tolower(template)) ~ "MNI152",
          grepl("^elderly", tolower(template)) ~ "Elderly",
          TRUE ~ as.character(template)
        ),
        template = factor(template, levels=c("MNI152","Elderly"))
      )
    tabs <- wb %>% split(.$Metric) %>% lapply(function(d){
      fit <- lmer(Value ~ Age*template + (1|ID), data=d, REML=FALSE)
      s   <- summary(fit); co <- coef(s)
      b   <- co["Age:templateElderly","Estimate"]
      se  <- co["Age:templateElderly","Std. Error"]
      p   <- co["Age:templateElderly","Pr(>|t|)"]
      tibble(Metric=unique(d$Metric),
             Coef=b, CI_lower=b-1.96*se, CI_upper=b+1.96*se, p_raw=p)
    }) %>% bind_rows() %>% mutate(p_BH = p.adjust(p_raw, method="BH"))
    write_csv(tabs, file.path(out_dir, "metrics_interaction_and_trends.csv"))
  } else {
    warning("Metrics CSV found but columns are not standard; skipping metrics table.")
  }
}

# -------------------- session info --------------------
writeLines(capture.output(sessionInfo()),
           file.path(root_dir, "analysis", "sessionInfo.txt"))

message("Done. Outputs saved to: ", out_dir)