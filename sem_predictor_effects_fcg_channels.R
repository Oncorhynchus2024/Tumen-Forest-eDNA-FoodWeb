# ============================================================
# STRICT SCI-Level Hierarchical Bayesian Piecewise SEM Pipeline
# FINAL OPTIMIZED VERSION (3 Mediators, 4 Shannon Layer)
# Added Controls: River (Random Intercept) + Width (Hidden Fixed Control)
# MEMORY SAFE VERSION: Fits sub-models separately to prevent C++ Stack Overflow
#
# FINAL VERSION:
# 1) adapt_delta and max_treedepth are truly passed into brm()
# 2) Added strict publication-level priors
# 3) River Width is still controlled in models, but hidden from figures
# 4) In AB figures, paths with |est| < 0.20 are hidden
# 5) AB path labels enlarged to >= 12 pt, closer to edges, and layout tightened
# ============================================================

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
set.seed(123)

# ------------------------------------------------------------
# Script overview
# ------------------------------------------------------------
# This script is part of the open research code archive for the Tumen River
# Basin study. It estimates hierarchical Bayesian piecewise SEMs linking the
# forest cover gradient (FCG), environmental mediators, multitrophic Shannon
# diversity, and food-web topology.
#
# Inputs expected in the working directory:
#   - abund.txt: genus-by-site abundance or presence table.
#   - genus_ffg.txt: mapping between genera and functional feeding groups.
#   - ffg_edges.txt: directed feeding links among functional feeding groups.
#   - env.txt: site-level environmental table containing Site, river, FCG
#              or the previous forest-gradient column name, river_width, TSM, and water-chemistry variables.
#   - spe_bacteria.txt, spe_protist.txt, spe_Metazoa.txt, spe_fish.txt:
#              ASV/OTU tables used to calculate Shannon diversity.
#
# Main outputs:
#   - Bayesian SEM diagnostics for each sub-model.
#   - Path summary tables with posterior medians, credible intervals, pd, and
#     significance markers.
#   - Publication-oriented SEM figures with river width controlled but hidden.
#
# Note on variable naming:
#   - FCG is used throughout this script as the formal manuscript term.
#   - A previous forest-gradient column name in env.txt is still accepted and internally renamed
#     to FCG for backward compatibility.
#
# This version focuses on predictor-centered cascading effects for panels C and D.

# -----------------------------
# 0) RUN MODE & MCMC SETTINGS
# -----------------------------
RUN_MODE <- "final"   # Use "final" for publication-quality figures; use "fast" for testing.

MODULARITY_N <- ifelse(RUN_MODE == "final", 1000, 30)
MCMC_ITER    <- ifelse(RUN_MODE == "final", 10000, 1000)
MCMC_WARMUP  <- ifelse(RUN_MODE == "final", 5000, 500)
MCMC_CHAINS  <- 4

CRI_LEVEL <- 0.95
PD_AUX_THRESHOLD <- 0.975

ADAPT_DELTA   <- ifelse(RUN_MODE == "final", 0.99, 0.95)
MAX_TREEDEPTH <- ifelse(RUN_MODE == "final", 15, 12)

# -----------------------------
# 1) Packages
# -----------------------------
pkgs <- c(
  "tidyverse","janitor","igraph","ggraph","tidygraph",
  "scales","ggrepel","patchwork","MASS","brms",
  "posterior","bayesplot","vegan","rstan"
)
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install) > 0){
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(MASS)
  library(brms)
  library(posterior)
  library(bayesplot)
  library(vegan)
  library(rstan)
})

options(mc.cores = max(1, parallel::detectCores() - 1))
rstan::rstan_options(auto_write = TRUE)

# -----------------------------
# 2) User settings
# -----------------------------
ABUND_FILE <- "abund.txt"
GENUS_FFG  <- "genus_ffg.txt"
FFG_EDGES  <- "ffg_edges.txt"
ENV_FILE   <- "env.txt"
SPE_BACTERIA_FILE <- "spe_bacteria.txt"
SPE_PROTIST_FILE  <- "spe_protist.txt"
SPE_METAZOA_FILE  <- "spe_Metazoa.txt"
SPE_FISH_FILE     <- "spe_fish.txt"

OUT_DIR <- file.path(getwd(), "SEM_Output_3Med_4Shan_Controlled_HideWidth_StrictPrior_Cut020_PredictorEffects")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_TOTAL_ABUND      <- 1
USE_PRESENCE_ONLY    <- TRUE
EDGE_DIRECTION       <- "prey_to_predator"   # or "predator_to_prey"

ALWAYS_INCLUDE_BASAL_RESOURCES <- TRUE
BASAL_ALWAYS_CAND <- c(
  "Plants","Fungi","Plankton","AquaticDetritus",
  "TerrestrialDetritus","Aquatic_Detritus",
  "Terrestrial_Detritus","Detritus","Algae"
)

NUTRIENT_VARS <- c("NH4","NO2","NO3","PO4")
PHYSCHEM_VARS <- c("T", "DO", "pH", "TDS", "Cond", "c", "Conductivity")  # Synchronized with the March 27 workflow; salinity is excluded.
TSM_VAR <- "TSM"

COL_POS <- "#2C7BB6"
COL_NEG <- "#D7191C"

# Only paths with |est| >= 0.20 are displayed in panels A and B.
EDGE_DISPLAY_THRESHOLD_ROW1 <- 0.20

# -----------------------------
# 2b) STRICT PRIORS
# -----------------------------
PRIORS_STRICT <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sigma"),
  prior(exponential(1), class = "sd")
)

# -----------------------------
# 3) Helpers
# -----------------------------
calc_skewness <- function(x){
  x <- x[is.finite(x)]
  n <- length(x)
  if(n < 3) return(NA_real_)
  s <- sd(x)
  if(!is.finite(s) || s == 0) return(NA_real_)
  sum((x - mean(x))^3 / n) / (s^3)
}

here_candidates <- c(
  getwd(),
  tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA)
)
here_candidates <- unique(na.omit(here_candidates))

find_file <- function(fname){
  for(d in here_candidates){
    f <- file.path(d, fname)
    if(file.exists(f)) return(f)
  }
  stop("Cannot find file: ", fname)
}

read_tsv_safe <- function(path){
  read.delim(
    path, sep = "\t", header = TRUE, check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
}

z <- function(x){
  x <- suppressWarnings(as.numeric(x))
  mu  <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if(!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

logit01 <- function(p){
  qlogis(pmin(pmax(suppressWarnings(as.numeric(p)), 0), 1) * 0.998 + 0.001)
}

save_pdf_tmp <- function(path, width, height, plot_fun){
  tmp <- paste0(path, ".tmp.pdf")
  if(file.exists(tmp)) suppressWarnings(file.remove(tmp))
  if(isTRUE(capabilities("cairo"))){
    grDevices::cairo_pdf(tmp, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::pdf(tmp, width = width, height = height, onefile = TRUE)
  }
  plot_fun()
  grDevices::dev.off()
  if(file.exists(path)) suppressWarnings(file.remove(path))
  ok <- file.rename(tmp, path)
  if(!isTRUE(ok)) warning("PDF overwrite failed: ", path)
}

nice_label <- function(x, label_map){
  y <- unname(label_map[x])
  ifelse(is.na(y) | y == "", x, y)
}

cr_interval_sig <- function(ci_lo, ci_hi){
  is.finite(ci_lo) && is.finite(ci_hi) && (ci_lo > 0 || ci_hi < 0)
}

stars_from_cri <- function(x){
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if(length(x) < 10) return("")
  
  q <- quantile(
    x,
    probs = c(0.0005, 0.005, 0.025, 0.975, 0.995, 0.9995),
    na.rm = TRUE, names = FALSE
  )
  
  if(q[1] > 0 || q[6] < 0) return("***")
  if(q[2] > 0 || q[5] < 0) return("**")
  if(q[3] > 0 || q[4] < 0) return("*")
  ""
}

summarize_posterior_vec <- function(x){
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  
  if(length(x) < 2){
    return(tibble(
      est = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      pd = NA_real_, sig = FALSE, pd_sig = FALSE, stars = ""
    ))
  }
  
  ci_lo <- quantile(x, 0.025, na.rm = TRUE, names = FALSE)
  ci_hi <- quantile(x, 0.975, na.rm = TRUE, names = FALSE)
  pd    <- max(mean(x > 0), mean(x < 0))
  
  sig    <- cr_interval_sig(ci_lo, ci_hi)
  pd_sig <- is.finite(pd) && pd >= PD_AUX_THRESHOLD
  stars  <- stars_from_cri(x)
  
  tibble(
    est   = median(x),
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    pd    = pd,
    sig   = sig,
    pd_sig = pd_sig,
    stars = stars
  )
}

do_pc1 <- function(df0, vars_raw, out_name){
  vars <- janitor::make_clean_names(vars_raw)
  mat <- df0[, vars[vars %in% names(df0)], drop = FALSE]
  mat <- mat[, sapply(mat, function(x) sum(is.finite(x)) >= 0.8 * nrow(mat)), drop = FALSE]
  
  if(ncol(mat) < 2){
    df0[[out_name]] <- 0
    return(df0)
  }
  
  for(nm in names(mat)){
    x <- suppressWarnings(as.numeric(mat[[nm]]))
    x[!is.finite(x)] <- median(x[is.finite(x)], na.rm = TRUE)
    mat[[nm]] <- x
  }
  
  pc1 <- as.numeric(prcomp(scale(as.matrix(mat)))$x[, 1])
  cc <- suppressWarnings(cor(pc1, rowSums(mat, na.rm = TRUE), use = "pairwise.complete.obs"))
  if(is.finite(cc) && cc < 0) pc1 <- -pc1
  
  df0[[out_name]] <- pc1
  df0
}

# ---- PhysChem PC1 helper functions synchronized with the updated March 27 workflow ----
replace_subscript_digits <- function(x){
    chartr("₀₁₂₃₄₅₆₇₈₉", "0123456789", x)
}

normalize_key <- function(x){
    x <- replace_subscript_digits(x)
    x <- tolower(trimws(x))
    x <- gsub("[^a-z0-9]+", "", x)
    x
}

find_col_relaxed <- function(df_names, candidates){
    keys <- normalize_key(df_names)
    cand_keys <- unique(normalize_key(candidates))
    for(ck in cand_keys){
        hit <- which(keys == ck)
        if(length(hit) > 0) return(df_names[hit[1]])
    }
    NA_character_
}

compute_pc1_generic <- function(df0, selected_cols, align_to = NULL, miss_keep = NULL){
    selected_cols <- unique(stats::na.omit(selected_cols))
    n <- nrow(df0)

    if(n < 3 || length(selected_cols) < 2){
        return(list(
            scores = rep(0, n),
            loadings = NULL,
            used_cols = selected_cols,
            var_prop = NA_real_
        ))
    }

    mat <- as.data.frame(lapply(df0[, selected_cols, drop = FALSE], function(v) suppressWarnings(as.numeric(v))))

    if(!is.null(miss_keep)){
        keep_miss <- sapply(mat, function(v) mean(is.finite(v)) >= miss_keep)
        mat <- mat[, keep_miss, drop = FALSE]
    }

    if(ncol(mat) < 2){
        return(list(
            scores = rep(0, n),
            loadings = NULL,
            used_cols = colnames(mat),
            var_prop = NA_real_
        ))
    }

    for(nm in names(mat)){
        v <- suppressWarnings(as.numeric(mat[[nm]]))
        med <- median(v[is.finite(v)], na.rm = TRUE)
        if(!is.finite(med)) med <- 0
        v[!is.finite(v)] <- med
        mat[[nm]] <- v
    }

    keep_sd <- sapply(mat, function(v){
        s <- sd(v, na.rm = TRUE)
        is.finite(s) && s > 0
    })
    mat <- mat[, keep_sd, drop = FALSE]

    if(ncol(mat) < 2){
        return(list(
            scores = rep(0, n),
            loadings = NULL,
            used_cols = colnames(mat),
            var_prop = NA_real_
        ))
    }

    pc <- prcomp(scale(as.matrix(mat)))
    pc1 <- as.numeric(pc$x[, 1])

    # Synchronization with the March 27 workflow: align PC1 to temperature first; if temperature is unavailable, align to rowSums(mat).
    if(is.null(align_to)){
        anchor <- rowSums(mat, na.rm = TRUE)
    } else {
        anchor <- suppressWarnings(as.numeric(align_to))
        med <- median(anchor[is.finite(anchor)], na.rm = TRUE)
        if(!is.finite(med)) med <- 0
        anchor[!is.finite(anchor)] <- med
    }

    flip <- 1
    cc <- suppressWarnings(cor(pc1, anchor, use = "pairwise.complete.obs"))
    if(is.finite(cc) && cc < 0) flip <- -1

    list(
        scores = pc1 * flip,
        loadings = pc$rotation[, 1] * flip,
        used_cols = colnames(mat),
        var_prop = (pc$sdev^2) / sum(pc$sdev^2)
    )
}




compute_site_shannon_from_spe <- function(path){
  spe <- read_tsv_safe(find_file(path))
  if(ncol(spe) < 3) stop("Species file has too few columns: ", path)
  
  site_hits_cols <- sum(names(spe)[-1] %in% env$Site)
  site_hits_rows <- sum(as.character(spe[[1]]) %in% env$Site)
  
  if(site_hits_cols >= site_hits_rows){
    rownames(spe) <- as.character(spe[[1]])
    spe <- spe[, -1, drop = FALSE]
    keep_sites <- intersect(colnames(spe), env$Site)
    if(length(keep_sites) < 2) stop("Too few overlapping sites in file: ", path)
    mat <- t(as.matrix(spe[, keep_sites, drop = FALSE]))
  } else {
    names(spe)[1] <- "Site"
    keep_sites <- intersect(as.character(spe$Site), env$Site)
    spe <- spe[spe$Site %in% keep_sites, , drop = FALSE]
    rownames(spe) <- as.character(spe$Site)
    spe <- spe[, setdiff(names(spe), "Site"), drop = FALSE]
    mat <- as.matrix(spe)
  }
  
  storage.mode(mat) <- "numeric"
  mat[!is.finite(mat)] <- 0
  sh <- vegan::diversity(mat, index = "shannon")
  tibble(Site = rownames(mat), Shannon = as.numeric(sh))
}

prepare_complete_case_data <- function(formulas_list, dat){
  vars_needed <- unique(unlist(lapply(formulas_list, function(f) {
    all.vars(as.formula(f))
  })))
  dat_use <- dat %>%
    filter(if_all(all_of(vars_needed), ~ !is.na(.x)))
  list(data = dat_use, vars_needed = vars_needed)
}

write_method_note <- function(out_dir){
    txt <- c(
        "IMPORTANT INTERPRETATION NOTE",
        "",
        "1) 3 mediators implemented: TSM, Nutrient PC1, PhysChem PC1.",
        "2) Added Shannon layer for bacteria, protist, metazoa and fish.",
        "3) Model controls for river_width (fixed effect) and river (random intercept).",
        "4) River Width is hidden from figures but remains in all fitted equations.",
        "5) Strict priors added: Intercept student_t(3,0,2.5), slopes normal(0,1), sigma/sd exponential(1).",
        "6) AB figures only show paths with |est| >= 0.20.",
        "7) CD figures show predictor-centered cascading effects (FCG, TSM, Nutrient PC1, PhysChem PC1, and four Shannon predictors).",
        "8) Stan control truly passed into brm(): adapt_delta and max_treedepth."
    )
    writeLines(txt, con = file.path(out_dir, "README_Method_and_Interpretation.txt"))
}

save_brms_diagnostics <- function(fit, fit_name, n_obs, out_dir){
  diag_dir <- file.path(out_dir, paste0("Diagnostics_", fit_name))
  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)
  
  writeLines(
    capture.output(summary(fit)),
    con = file.path(diag_dir, paste0(fit_name, "_summary.txt"))
  )
  
  draws_mat <- posterior::as_draws_matrix(fit)
  par_keep <- colnames(draws_mat)
  par_keep <- par_keep[grepl("^b_|^sigma|^sd_", par_keep)]
  
  if(length(par_keep) > 0){
    mat_sub <- draws_mat[, par_keep, drop = FALSE]
    
    diag_tab <- tibble(
      parameter = par_keep,
      median   = apply(mat_sub, 2, median, na.rm = TRUE),
      ci_lo    = apply(mat_sub, 2, quantile, probs = 0.025, na.rm = TRUE),
      ci_hi    = apply(mat_sub, 2, quantile, probs = 0.975, na.rm = TRUE),
      rhat     = as.numeric(posterior::rhat(mat_sub)),
      ess_bulk = as.numeric(posterior::ess_bulk(mat_sub)),
      ess_tail = as.numeric(posterior::ess_tail(mat_sub))
    )
    
    write.table(
      diag_tab,
      file = file.path(diag_dir, paste0(fit_name, "_parameter_diagnostics.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
  
  nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
  if(!is.null(nuts) && nrow(nuts) > 0){
    write.table(
      nuts,
      file = file.path(diag_dir, paste0(fit_name, "_nuts_params.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )
  }
}

fit_layer_model_separate <- function(formulas_list, dat, fit_name, out_dir){
  
  prep <- prepare_complete_case_data(formulas_list, dat)
  dat_use <- prep$data
  
  if(nrow(dat_use) < 20){
    stop(fit_name, ": complete-case sample size too small (n < 20).")
  }
  
  cat("\n  -> Fitting", length(formulas_list), "sub-models separately...\n")
  
  all_draws <- list()
  
  for(i in seq_along(formulas_list)){
    f_str <- formulas_list[[i]]
    resp <- trimws(strsplit(f_str, "~")[[1]][1])
    resp_clean <- gsub("_", "", resp)
    
    cat(sprintf("     [%d/%d] Fitting sub-model: %s\n", i, length(formulas_list), resp))
    
    fit <- brm(
      formula = bf(as.formula(f_str)),
      data    = dat_use,
      family  = gaussian(),
      prior   = PRIORS_STRICT,
      chains  = MCMC_CHAINS,
      cores   = max(1, parallel::detectCores() - 1),
      iter    = MCMC_ITER,
      warmup  = MCMC_WARMUP,
      backend = "rstan",
      silent  = 2,
      refresh = 0,
      seed    = 123 + i,
      control = list(
        adapt_delta   = ADAPT_DELTA,
        max_treedepth = MAX_TREEDEPTH
      )
    )
    
    save_brms_diagnostics(fit, paste0(fit_name, "_", resp), nrow(dat_use), out_dir)
    
    d <- posterior::as_draws_df(fit)
    names(d) <- gsub("^b_", paste0("b_", resp_clean, "_"), names(d))
    d <- d[, grepl(paste0("^b_", resp_clean, "_"), names(d)), drop = FALSE]
    
    all_draws[[resp]] <- d
  }
  
  draws_combined <- bind_cols(all_draws)
  
  list(
    draws_df = draws_combined,
    data = dat_use,
    vars_used = prep$vars_needed
  )
}

# -----------------------------
# 4) Load Data
# -----------------------------
write_method_note(OUT_DIR)

env_raw <- read_tsv_safe(find_file(ENV_FILE))
names(env_raw) <- trimws(names(env_raw))
if("PO₄" %in% names(env_raw)) env_raw <- env_raw %>% dplyr::rename(PO4 = `PO₄`)
env <- env_raw %>% janitor::clean_names()

if("site" %in% names(env)) names(env)[names(env) == "site"] <- "Site"

# Standardize the forest-cover gradient column name.
# The current manuscript uses FCG. A previous forest-gradient column name is also accepted
# to keep older input tables compatible with this script.
if("fcg" %in% names(env)){
  names(env)[names(env) == "fcg"] <- "FCG"
} else if("fi" %in% names(env)){
  names(env)[names(env) == "fi"] <- "FCG"
} else {
  stop("The environmental table must contain an 'FCG' column or the accepted previous forest-gradient column name.")
}

env$Site <- as.character(env$Site)
for(nm in names(env)){
  if(!nm %in% c("Site","river")) env[[nm]] <- suppressWarnings(as.numeric(env[[nm]]))
}

# ---- Compute PhysChem PC1 at the environmental-data stage using the updated March 27 workflow ----
t_col <- find_col_relaxed(names(env), c("T"))

physchem_cols <- c(
    find_col_relaxed(names(env), c("T")),
    find_col_relaxed(names(env), c("DO")),
    find_col_relaxed(names(env), c("pH")),
    find_col_relaxed(names(env), c("TDS")),
    find_col_relaxed(names(env), c("Cond", "c", "Conductivity"))
)

phys_res <- compute_pc1_generic(
    df0 = env,
    selected_cols = physchem_cols,
    align_to = if(!is.na(t_col)) env[[t_col]] else NULL,
    miss_keep = 0.8
)

env$physchem_pc1 <- phys_res$scores



gffg <- read_tsv_safe(find_file(GENUS_FFG)) %>%
  janitor::clean_names() %>%
  transmute(Genus = as.character(genus), FFG = as.character(ffg)) %>%
  tidyr::separate_rows(FFG, sep = ";") %>%
  mutate(FFG = stringr::str_squish(FFG)) %>%
  filter(!is.na(Genus), !is.na(FFG), Genus != "", FFG != "") %>%
  distinct()

ffg_edges_raw <- read_tsv_safe(find_file(FFG_EDGES))
names(ffg_edges_raw) <- janitor::make_clean_names(names(ffg_edges_raw))
ffg_edges <- tibble(
  FromFFG = as.character(ffg_edges_raw[[grep("^from(_)?ffg$", names(ffg_edges_raw), value = TRUE)[1]]]),
  ToFFG   = as.character(ffg_edges_raw[[grep("^to(_)?ffg$",   names(ffg_edges_raw), value = TRUE)[1]]])
) %>%
  mutate(FromFFG = trimws(FromFFG), ToFFG = trimws(ToFFG)) %>%
  filter(FromFFG != "", ToFFG != "", FromFFG != ToFFG) %>%
  distinct()

BASAL_ALWAYS <- character(0)
if(isTRUE(ALWAYS_INCLUDE_BASAL_RESOURCES)){
  BASAL_ALWAYS <- BASAL_ALWAYS_CAND[
    tolower(BASAL_ALWAYS_CAND) %in% tolower(unique(c(ffg_edges$FromFFG, ffg_edges$ToFFG)))
  ]
}

abund_raw <- read_tsv_safe(find_file(ABUND_FILE)) %>%
  dplyr::select(-any_of(c("tax", "Tax", "TAX", "taxonomy", "Taxonomy")))

if(sum(colnames(abund_raw) %in% env$Site) >= sum(as.character(abund_raw[[1]]) %in% env$Site)){
  names(abund_raw)[1] <- "Genus"
  abund_long <- tidyr::pivot_longer(
    abund_raw,
    cols = all_of(intersect(names(abund_raw)[-1], env$Site)),
    names_to = "Site", values_to = "Abund"
  )
} else {
  names(abund_raw)[1] <- "Site"
  abund_long <- tidyr::pivot_longer(
    abund_raw,
    cols = all_of(setdiff(names(abund_raw), "Site")),
    names_to = "Genus", values_to = "Abund"
  )
}

abund_long <- abund_long %>%
  mutate(
    Genus = as.character(Genus),
    Site  = as.character(Site),
    Abund = suppressWarnings(as.numeric(Abund))
  ) %>%
  filter(Site %in% env$Site, is.finite(Abund)) %>%
  group_by(Site, Genus) %>%
  summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop") %>%
  filter(Abund >= MIN_TOTAL_ABUND)

if(USE_PRESENCE_ONLY) abund_long$Abund <- 1

genus_richness_df <- abund_long %>%
  distinct(Site, Genus) %>%
  count(Site, name = "Genus_richness")

site_ffg <- abund_long %>%
  inner_join(gffg, by = "Genus") %>%
  distinct(Site, FFG)

ffg_richness_df <- site_ffg %>%
  distinct(Site, FFG) %>%
  count(Site, name = "Observed_FFG_richness")

# -----------------------------
# 4b) Multi-group Shannon layer
# -----------------------------
cat("\nCalculating Shannon indices for bacteria / protist / metazoa / fish...\n")

shannon_bacteria <- compute_site_shannon_from_spe(SPE_BACTERIA_FILE) %>% rename(bacteria_shannon = Shannon)
shannon_protist  <- compute_site_shannon_from_spe(SPE_PROTIST_FILE)  %>% rename(protist_shannon  = Shannon)
shannon_metazoa  <- compute_site_shannon_from_spe(SPE_METAZOA_FILE)  %>% rename(metazoa_shannon  = Shannon)
shannon_fish     <- compute_site_shannon_from_spe(SPE_FISH_FILE)     %>% rename(fish_shannon     = Shannon)

shannon_df <- list(shannon_bacteria, shannon_protist, shannon_metazoa, shannon_fish) %>%
  purrr::reduce(full_join, by = "Site")

# -----------------------------
# 5) Hierarchical Topology Calculator
# -----------------------------
calc_topology_hierarchical <- function(g_dir, mod_n){
  
  S <- vcount(g_dir)
  L <- ecount(g_dir)
  
  if(S < 2 || L < 1){
    return(tibble(
      FFG_richness = S, Connectance_dir = NA_real_, Degree_skewness = NA_real_,
      Mean_TL = NA_real_, Mean_Generality = NA_real_, Omnivory = NA_real_,
      Trophic_incoherence = NA_real_, Modularity = NA_real_,
      Nestedness_NODF = NA_real_, Niche_overlap = NA_real_
    ))
  }
  
  FFG_richness <- S
  Connectance_dir <- L / (S * (S - 1))
  
  deg_all <- degree(g_dir, mode = "all")
  Degree_skewness <- calc_skewness(deg_all)
  
  A <- as.matrix(igraph::as_adjacency_matrix(g_dir, sparse = FALSE))
  indeg_g <- colSums(A)
  
  W <- matrix(0, nrow = S, ncol = S)
  for(i in seq_len(S)){
    if(indeg_g[i] > 0) W[i, which(A[, i] > 0)] <- 1 / indeg_g[i]
  }
  
  TL <- as.numeric(MASS::ginv(diag(S) - W) %*% rep(1, S))
  TL[indeg_g == 0] <- 1
  Mean_TL <- mean(TL, na.rm = TRUE)
  
  consumers_idx <- which(indeg_g > 0)
  n_consumers <- length(consumers_idx)
  
  Mean_Generality <- ifelse(n_consumers > 0, mean(indeg_g[consumers_idx], na.rm = TRUE), NA_real_)
  
  Omnivory <- NA_real_
  if(n_consumers > 0){
    n_omni <- sum(sapply(consumers_idx, function(j){
      prey_tls <- TL[which(A[, j] > 0)]
      length(prey_tls) > 1 && (max(prey_tls, na.rm = TRUE) - min(prey_tls, na.rm = TRUE)) >= 2.0
    }))
    Omnivory <- n_omni / n_consumers
  }
  
  el <- as_edgelist(g_dir, names = FALSE)
  dx <- TL[el[,2]] - TL[el[,1]]
  dx <- dx[is.finite(dx)]
  Trophic_incoherence <- if(length(dx) > 0) sqrt(mean((dx - 1)^2)) else NA_real_
  
  Niche_overlap <- NA_real_
  if(n_consumers > 1){
    sim_mat <- igraph::similarity(g_dir, vids = consumers_idx, mode = "in", method = "jaccard")
    Niche_overlap <- mean(sim_mat[upper.tri(sim_mat)], na.rm = TRUE)
  }
  
  M <- (A[consumers_idx, , drop = FALSE] > 0) * 1
  Nestedness_NODF <- NA_real_
  if(nrow(M) > 1 && ncol(M) > 1){
    rs <- rowSums(M)
    cs <- colSums(M)
    M <- M[rs > 0, cs > 0, drop = FALSE]
    if(nrow(M) > 1 && ncol(M) > 1){
      rs <- rowSums(M); cs <- colSums(M)
      M  <- M[order(rs, decreasing = TRUE), order(cs, decreasing = TRUE), drop = FALSE]
      rs <- rowSums(M); cs <- colSums(M)
      
      acc <- 0; den <- 0
      for(i in 1:(nrow(M)-1)){
        for(j in (i+1):nrow(M)){
          if(rs[i] > rs[j] && rs[j] > 0){
            acc <- acc + sum(M[i,] * M[j,]) / rs[j]
            den <- den + 1
          }
        }
      }
      for(i in 1:(ncol(M)-1)){
        for(j in (i+1):ncol(M)){
          if(cs[i] > cs[j] && cs[j] > 0){
            acc <- acc + sum(M[,i] * M[,j]) / cs[j]
            den <- den + 1
          }
        }
      }
      Nestedness_NODF <- if(den == 0) 0 else acc / den
    }
  }
  
  g_und <- as_undirected(g_dir, mode = "collapse")
  Modularity <- if(ecount(g_und) > 0){
    mean(sapply(1:mod_n, function(x){
      set.seed(2025 + x)
      igraph::modularity(igraph::cluster_louvain(igraph::permute(g_und, sample.int(S))))
    }), na.rm = TRUE)
  } else {
    NA_real_
  }
  
  tibble(
    FFG_richness, Connectance_dir, Degree_skewness,
    Mean_TL, Mean_Generality, Omnivory, Trophic_incoherence,
    Modularity, Nestedness_NODF, Niche_overlap
  )
}

cat("\nCalculating FFG-level inferred metaweb metrics...\n")
sites_used <- sort(unique(site_ffg$Site))
pb <- txtProgressBar(min = 0, max = length(sites_used), style = 3)

metrics_df <- bind_rows(lapply(seq_along(sites_used), function(i){
  sid <- sites_used[i]
  nf  <- unique(c(site_ffg$FFG[site_ffg$Site == sid], BASAL_ALWAYS))
  ed  <- ffg_edges %>% filter(FromFFG %in% nf, ToFFG %in% nf)
  
  ed_dir <- if(EDGE_DIRECTION == "predator_to_prey"){
    ed %>% transmute(FromFFG = ToFFG, ToFFG = FromFFG)
  } else {
    ed
  }
  
  g <- graph_from_data_frame(
    ed_dir %>% rename(from = FromFFG, to = ToFFG),
    directed = TRUE,
    vertices = data.frame(name = nf)
  )
  
  setTxtProgressBar(pb, i)
  
  tibble(Site = sid) %>%
    bind_cols(calc_topology_hierarchical(g, MODULARITY_N))
}))
close(pb)

metrics_df <- metrics_df %>%
  left_join(genus_richness_df, by = "Site") %>%
  left_join(ffg_richness_df, by = "Site") %>%
  left_join(shannon_df, by = "Site") %>%
  left_join(env, by = "Site")

# -----------------------------
# 6) Data Scaling & Random Effect Prep
# -----------------------------
df <- metrics_df %>%
  clean_names() %>%
  rename(FCG = fcg) %>%
  filter(is.finite(FCG))

df <- do_pc1(df, NUTRIENT_VARS, "nutrient_pc1")
# PhysChem PC1 is not recalculated here; the synchronized env$physchem_pc1 computed above is used directly.

tsm_col <- janitor::make_clean_names(TSM_VAR)
df$tsm <- if(tsm_col %in% names(df)) as.numeric(df[[tsm_col]]) else NA_real_

width_col <- if("river_width" %in% names(df)) "river_width" else if("riverwidth" %in% names(df)) "riverwidth" else NA
if(!is.na(width_col)){
  df$river_width <- as.numeric(df[[width_col]])
  df$river_width_z <- z(df$river_width)
} else {
  stop("river_width variable not found in env.txt")
}

if("river" %in% names(df)){
  df$river <- as.factor(df$river)
} else {
  stop("river variable not found in env.txt (Required for random intercept)")
}

df$FCG_z           <- z(df$FCG)
df$tsm_z          <- z(df$tsm)
df$nutrient_pc1_z <- z(df$nutrient_pc1)
df$physchem_pc1_z <- z(df$physchem_pc1)
df$bacteria_shannon_z <- z(df$bacteria_shannon)
df$protist_shannon_z  <- z(df$protist_shannon)
df$metazoa_shannon_z  <- z(df$metazoa_shannon)
df$fish_shannon_z     <- z(df$fish_shannon)

vars_logit  <- c("connectance_dir", "omnivory", "niche_overlap", "nestedness_nodf", "modularity")
vars_log    <- c("genus_richness", "mean_generality", "trophic_incoherence")
vars_linear <- c("degree_skewness", "mean_tl")

for(v in c(vars_logit, vars_log, vars_linear)){
  if(v %in% names(df)){
    if(v %in% vars_logit){
      df[[paste0(v, "_z")]] <- z(logit01(df[[v]]))
    } else if(v %in% vars_log){
      df[[paste0(v, "_z")]] <- z(log(pmax(df[[v]], 1e-6)))
    } else {
      df[[paste0(v, "_z")]] <- z(df[[v]])
    }
  }
}

dat_brms <- df %>% filter(is.finite(FCG_z) & is.finite(river_width_z) & !is.na(river))

# -----------------------------
# 7) Bayesian helpers
# -----------------------------
extract_paths_with_draws <- function(draws, from_vars, to_vars){
  
  sum_list  <- list()
  draw_list <- list()
  
  for(to in to_vars){
    to_clean <- gsub("_", "", to)
    for(from in from_vars){
      col_name <- paste0("b_", to_clean, "_", from)
      if(col_name %in% names(draws)){
        x  <- as.numeric(draws[[col_name]])
        sm <- summarize_posterior_vec(x)
        
        sum_list[[length(sum_list) + 1]] <- tibble(
          from  = from,
          to    = to,
          est   = sm$est,
          ci_lo = sm$ci_lo,
          ci_hi = sm$ci_hi,
          pd    = sm$pd,
          sig   = sm$sig,
          pd_sig = sm$pd_sig,
          stars = sm$stars
        )
        
        draw_list[[paste0(gsub("_z$", "", from), "->", gsub("_z$", "", to))]] <- x
      }
    }
  }
  
  if(length(sum_list) == 0){
    sum_df <- tibble(
      from = character(), to = character(),
      est = numeric(), ci_lo = numeric(), ci_hi = numeric(),
      pd = numeric(), sig = logical(), pd_sig = logical(), stars = character()
    )
  } else {
    sum_df <- bind_rows(sum_list)
  }
  
  list(
    summary = sum_df,
    draws   = draw_list,
    n_draws = if(length(draw_list) > 0) length(draw_list[[1]]) else 0
  )
}

calc_full_cascading_effects_posterior <- function(edge_draws, preds, endpoint, label_map){
  
  if(length(edge_draws) == 0) stop("No posterior edge draws available.")
  
  edge_keys  <- names(edge_draws)
  edge_parts <- strsplit(edge_keys, "->", fixed = TRUE)
  edge_from  <- vapply(edge_parts, `[`, character(1), 1)
  edge_to    <- vapply(edge_parts, `[`, character(1), 2)
  
  preds_clean    <- gsub("_z$", "", preds)
  endpoint_clean <- gsub("_z$", "", endpoint)
  
  all_nodes <- unique(c(edge_from, edge_to, preds_clean, endpoint_clean))
  N <- length(all_nodes)
  
  n_draws <- length(edge_draws[[1]])
  
  direct_store   <- matrix(0, nrow = n_draws, ncol = length(preds_clean), dimnames = list(NULL, preds_clean))
  indirect_store <- matrix(0, nrow = n_draws, ncol = length(preds_clean), dimnames = list(NULL, preds_clean))
  total_store    <- matrix(0, nrow = n_draws, ncol = length(preds_clean), dimnames = list(NULL, preds_clean))
  
  for(d in seq_len(n_draws)){
    A <- matrix(0, nrow = N, ncol = N, dimnames = list(all_nodes, all_nodes))
    for(k in seq_along(edge_keys)){
      A[edge_from[k], edge_to[k]] <- edge_draws[[k]][d]
    }
    
    DIR <- A
    IND <- matrix(0, nrow = N, ncol = N, dimnames = list(all_nodes, all_nodes))
    
    if(N >= 3){
      Apow <- A %*% A
      IND  <- IND + Apow
      if(N >= 4){
        for(step in 3:(N - 1)){
          Apow <- Apow %*% A
          IND  <- IND + Apow
        }
      }
    }
    
    TOT <- DIR + IND
    
    for(i in seq_along(preds_clean)){
      pc <- preds_clean[i]
      if(pc %in% rownames(DIR) && endpoint_clean %in% colnames(DIR)){
        direct_store[d, i]   <- DIR[pc, endpoint_clean]
        indirect_store[d, i] <- IND[pc, endpoint_clean]
        total_store[d, i]    <- TOT[pc, endpoint_clean]
      }
    }
  }
  
  long_list <- list()
  for(i in seq_along(preds)){
    plab <- nice_label(preds[i], label_map)
    olab <- nice_label(endpoint, label_map)
    
    sm_dir <- summarize_posterior_vec(direct_store[, i])
    sm_ind <- summarize_posterior_vec(indirect_store[, i])
    sm_tot <- summarize_posterior_vec(total_store[, i])
    
    long_list[[length(long_list) + 1]] <- tibble(
      Predictor  = plab, Outcome = olab, EffectType = "Direct",
      est = sm_dir$est, ci_lo = sm_dir$ci_lo, ci_hi = sm_dir$ci_hi,
      pd = sm_dir$pd, sig = sm_dir$sig, pd_sig = sm_dir$pd_sig, stars = sm_dir$stars
    )
    long_list[[length(long_list) + 1]] <- tibble(
      Predictor  = plab, Outcome = olab, EffectType = "Indirect",
      est = sm_ind$est, ci_lo = sm_ind$ci_lo, ci_hi = sm_ind$ci_hi,
      pd = sm_ind$pd, sig = sm_ind$sig, pd_sig = sm_ind$pd_sig, stars = sm_ind$stars
    )
    long_list[[length(long_list) + 1]] <- tibble(
      Predictor  = plab, Outcome = olab, EffectType = "Total",
      est = sm_tot$est, ci_lo = sm_tot$ci_lo, ci_hi = sm_tot$ci_hi,
      pd = sm_tot$pd, sig = sm_tot$sig, pd_sig = sm_tot$pd_sig, stars = sm_tot$stars
    )
  }
  
  effects_long <- bind_rows(long_list)
  
  effects_draws <- bind_rows(lapply(seq_along(preds), function(i){
    tibble(
      draw      = seq_len(n_draws),
      Predictor = nice_label(preds[i], label_map),
      Outcome   = nice_label(endpoint, label_map),
      Direct    = direct_store[, i],
      Indirect  = indirect_store[, i],
      Total     = total_store[, i]
    )
  }))
  
  list(summary_long = effects_long, draws_long = effects_draws)
}


get_path_products <- function(A, start, end, visited = character()){
    if(!(start %in% rownames(A)) || !(end %in% colnames(A))) return(numeric(0))
    visited <- c(visited, start)
    nbrs <- colnames(A)[which(abs(A[start, ]) > 0)]
    if(length(nbrs) == 0) return(numeric(0))
    out <- numeric(0)
    for(nb in nbrs){
        w <- A[start, nb]
        if(nb == end){
            out <- c(out, w)
        } else if(!(nb %in% visited)){
            sub <- get_path_products(A, nb, end, visited)
            if(length(sub) > 0) out <- c(out, w * sub)
        }
    }
    out
}

calc_fi_channel_effects_posterior <- function(edge_draws, fi_node, channels, endpoint, label_map){
    if(length(edge_draws) == 0) stop("No posterior edge draws available.")
    edge_keys  <- names(edge_draws)
    edge_parts <- strsplit(edge_keys, "->", fixed = TRUE)
    edge_from  <- vapply(edge_parts, `[`, character(1), 1)
    edge_to    <- vapply(edge_parts, `[`, character(1), 2)
    all_nodes <- unique(c(edge_from, edge_to, fi_node, channels, endpoint))
    n_draws <- length(edge_draws[[1]])
    
    direct_fi <- rep(0, n_draws)
    channel_direct <- matrix(0, nrow = n_draws, ncol = length(channels), dimnames = list(NULL, channels))
    channel_indirect <- matrix(0, nrow = n_draws, ncol = length(channels), dimnames = list(NULL, channels))
    channel_total <- matrix(0, nrow = n_draws, ncol = length(channels), dimnames = list(NULL, channels))
    
    for(d in seq_len(n_draws)){
        A <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes), dimnames = list(all_nodes, all_nodes))
        for(k in seq_along(edge_keys)) A[edge_from[k], edge_to[k]] <- edge_draws[[k]][d]
        
        direct_fi[d] <- if(fi_node %in% rownames(A) && endpoint %in% colnames(A)) A[fi_node, endpoint] else 0
        
        for(i in seq_along(channels)){
            ch <- channels[i]
            w_fi_ch <- if(fi_node %in% rownames(A) && ch %in% colnames(A)) A[fi_node, ch] else 0
            if(!is.finite(w_fi_ch) || w_fi_ch == 0){
                channel_direct[d, i] <- 0
                channel_indirect[d, i] <- 0
                channel_total[d, i] <- 0
            } else {
                sub_paths <- get_path_products(A, ch, endpoint)
                tot <- if(length(sub_paths) > 0) w_fi_ch * sum(sub_paths) else 0
                dir2 <- if(ch %in% rownames(A) && endpoint %in% colnames(A)) w_fi_ch * A[ch, endpoint] else 0
                channel_direct[d, i] <- dir2
                channel_total[d, i] <- tot
                channel_indirect[d, i] <- tot - dir2
            }
        }
    }
    
    fi_indirect_draws <- rowSums(channel_total)
    fi_total_draws <- direct_fi + fi_indirect_draws
    
    ordered_labels <- c("FCG overall", paste0("FCG via ", vapply(channels, function(x) nice_label(paste0(x, "_z"), label_map), character(1))))
    
    long_list <- list()
    
    sm_fi_dir <- summarize_posterior_vec(direct_fi)
    sm_fi_ind <- summarize_posterior_vec(fi_indirect_draws)
    sm_fi_tot <- summarize_posterior_vec(fi_total_draws)
    
    long_list[[length(long_list) + 1]] <- tibble(
        Predictor = "FCG overall", Outcome = nice_label(paste0(endpoint, "_z"), label_map), EffectType = "Direct",
        est = sm_fi_dir$est, ci_lo = sm_fi_dir$ci_lo, ci_hi = sm_fi_dir$ci_hi,
        pd = sm_fi_dir$pd, sig = sm_fi_dir$sig, pd_sig = sm_fi_dir$pd_sig, stars = sm_fi_dir$stars
    )
    long_list[[length(long_list) + 1]] <- tibble(
        Predictor = "FCG overall", Outcome = nice_label(paste0(endpoint, "_z"), label_map), EffectType = "Indirect",
        est = sm_fi_ind$est, ci_lo = sm_fi_ind$ci_lo, ci_hi = sm_fi_ind$ci_hi,
        pd = sm_fi_ind$pd, sig = sm_fi_ind$sig, pd_sig = sm_fi_ind$pd_sig, stars = sm_fi_ind$stars
    )
    long_list[[length(long_list) + 1]] <- tibble(
        Predictor = "FCG overall", Outcome = nice_label(paste0(endpoint, "_z"), label_map), EffectType = "Total",
        est = sm_fi_tot$est, ci_lo = sm_fi_tot$ci_lo, ci_hi = sm_fi_tot$ci_hi,
        pd = sm_fi_tot$pd, sig = sm_fi_tot$sig, pd_sig = sm_fi_tot$pd_sig, stars = sm_fi_tot$stars
    )
    
    for(i in seq_along(channels)){
        lab <- paste0("FCG via ", nice_label(paste0(channels[i], "_z"), label_map))
        sm_dir <- summarize_posterior_vec(channel_direct[, i])
        sm_ind <- summarize_posterior_vec(channel_indirect[, i])
        sm_tot <- summarize_posterior_vec(channel_total[, i])
        long_list[[length(long_list) + 1]] <- tibble(
            Predictor = lab, Outcome = nice_label(paste0(endpoint, "_z"), label_map), EffectType = "Direct",
            est = sm_dir$est, ci_lo = sm_dir$ci_lo, ci_hi = sm_dir$ci_hi,
            pd = sm_dir$pd, sig = sm_dir$sig, pd_sig = sm_dir$pd_sig, stars = sm_dir$stars
        )
        long_list[[length(long_list) + 1]] <- tibble(
            Predictor = lab, Outcome = nice_label(paste0(endpoint, "_z"), label_map), EffectType = "Indirect",
            est = sm_ind$est, ci_lo = sm_ind$ci_lo, ci_hi = sm_ind$ci_hi,
            pd = sm_ind$pd, sig = sm_ind$sig, pd_sig = sm_ind$pd_sig, stars = sm_ind$stars
        )
        long_list[[length(long_list) + 1]] <- tibble(
            Predictor = lab, Outcome = nice_label(paste0(endpoint, "_z"), label_map), EffectType = "Total",
            est = sm_tot$est, ci_lo = sm_tot$ci_lo, ci_hi = sm_tot$ci_hi,
            pd = sm_tot$pd, sig = sm_tot$sig, pd_sig = sm_tot$pd_sig, stars = sm_tot$stars
        )
    }
    
    effects_long <- bind_rows(long_list) %>%
        mutate(Predictor = factor(Predictor, levels = ordered_labels),
               EffectType = factor(EffectType, levels = c("Direct", "Indirect", "Total"))) %>%
        arrange(Predictor, EffectType)
    
    draws_long <- bind_rows(
        tibble(draw = seq_len(n_draws), Predictor = "FCG overall", Outcome = nice_label(paste0(endpoint, "_z"), label_map),
               Direct = direct_fi, Indirect = fi_indirect_draws, Total = fi_total_draws),
        bind_rows(lapply(seq_along(channels), function(i){
            tibble(draw = seq_len(n_draws),
                   Predictor = paste0("FCG via ", nice_label(paste0(channels[i], "_z"), label_map)),
                   Outcome = nice_label(paste0(endpoint, "_z"), label_map),
                   Direct = channel_direct[, i],
                   Indirect = channel_indirect[, i],
                   Total = channel_total[, i])
        }))
    )
    
    fi_total_summary <- summarize_posterior_vec(fi_total_draws)
    
    list(summary_long = effects_long, draws_long = draws_long, fi_total_draws = fi_total_draws, fi_total_summary = fi_total_summary)
}

# -----------------------------
# 8) Plot helpers
# -----------------------------
plot_hierarchical_sem <- function(paths, nodes_df, title, subtitle = NULL, edge_cutoff = 0.20){
  
  paths_clean <- paths %>%
    mutate(
      from = gsub("_z$", "", from),
      to   = gsub("_z$", "", to)
    ) %>%
    filter(
      is.finite(est),
      abs(est) >= edge_cutoff
    )
  
  edges <- paths_clean %>%
    filter(from %in% nodes_df$name, to %in% nodes_df$name) %>%
    mutate(
      col       = ifelse(sig, ifelse(est >= 0, COL_POS, COL_NEG), "grey75"),
      lty       = ifelse(sig, "solid", "dashed"),
      alpha_v   = ifelse(sig, 1.0, 0.6),
      w         = ifelse(sig, 0.7 + abs(est) * 2.2, 0.45),
      elab      = sprintf("%.2f%s", est, stars),
      font_face = ifelse(sig, "bold", "plain"),
      text_col  = ifelse(sig, "black", "grey40")
    )
  
  edges_mid <- edges %>%
    left_join(nodes_df, by = c("from" = "name")) %>% rename(x1 = x, y1 = y) %>%
    left_join(nodes_df, by = c("to" = "name"))   %>% rename(x2 = x, y2 = y) %>%
    mutate(
      xm = 0.45 * x1 + 0.55 * x2,
      ym = 0.45 * y1 + 0.55 * y2
    ) %>%
    filter(is.finite(x1) & is.finite(y1) & is.finite(x2) & is.finite(y2))
  
  edge_df_for_graph <- edges %>%
    transmute(from, to, col, alpha_v, lty, w)
  
  if(nrow(edge_df_for_graph) == 0){
    edge_df_for_graph <- tibble(
      from = character(),
      to = character(),
      col = character(),
      alpha_v = numeric(),
      lty = character(),
      w = numeric()
    )
  }
  
  g <- igraph::graph_from_data_frame(
    edge_df_for_graph,
    directed = TRUE,
    vertices = nodes_df
  )
  g_tbl <- tidygraph::as_tbl_graph(g)
  
  p <- ggraph(g_tbl, layout = "manual", x = x, y = y) +
    geom_edge_link(
      aes(edge_colour = col, edge_alpha = alpha_v, edge_width = w, edge_linetype = lty),
      arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"),
      end_cap = ggraph::circle(8.5, "mm"),
      start_cap = ggraph::circle(8.5, "mm"),
      lineend = "round", show.legend = FALSE
    ) +
    geom_node_label(
      aes(label = label, fill = fill_col),
      size = 5.2, fontface = "bold", color = "grey20",
      label.padding = grid::unit(0.38, "lines"),
      label.r = grid::unit(0.25, "lines"), show.legend = FALSE
    ) +
    ggrepel::geom_label_repel(
      data = edges_mid,
      aes(x = xm, y = ym, label = elab, color = text_col, fontface = font_face),
      inherit.aes = FALSE,
      size = 4.6,
      label.size = 0, fill = scales::alpha("white", 0.82),
      label.padding = grid::unit(0.10, "lines"),
      box.padding = 0.025, point.padding = 0.00,
      force = 0.15, force_pull = 1.0,
      min.segment.length = 0, seed = 123, show.legend = FALSE
    ) +
    scale_fill_identity() + scale_color_identity() +
    scale_edge_colour_identity() + scale_edge_alpha_identity() +
    scale_edge_width_identity() + scale_edge_linetype_identity() +
    coord_cartesian(
      xlim = range(nodes_df$x) + c(-0.32, 0.32),
      ylim = range(nodes_df$y) + c(-0.20, 0.20), clip = "off"
    ) +
    theme_void() +
    labs(title = title, subtitle = subtitle) +
    theme(
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 13, hjust = 0.5, color = "grey25", margin = margin(b = 10)),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  p
}

plot_effect_sizes_fi_channels <- function(eff_long, title, ylim_range = NULL){
    eff_plot <- eff_long %>%
        mutate(
            Predictor = factor(Predictor, levels = unique(as.character(Predictor))),
            EffectType = factor(EffectType, levels = c("Direct", "Indirect", "Total"))
        )
    
    dodge <- position_dodge(width = 0.75)
    p_bar <- ggplot(eff_plot, aes(x = Predictor, y = est, fill = EffectType)) +
        geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
        geom_col(position = dodge, width = 0.70, color = "grey25", linewidth = 0.35) +
        scale_fill_manual(values = c("Direct" = "#9ECAE1", "Indirect" = "#A1D99B", "Total" = "#FDD0A2")) +
        theme_bw(base_size = 14) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 20, hjust = 1),
            axis.text.y = element_text(size = 12, color = "black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(face = "bold", size = 14),
            plot.title = element_text(face = "bold", size = 18, hjust = 0.5, margin = margin(b = 8)),
            plot.margin = margin(8, 8, 8, 8)
        ) +
        labs(y = "Posterior effect size\n(median)", fill = "Effect type", title = title)
    if(!is.null(ylim_range) && length(ylim_range) == 2 && all(is.finite(ylim_range))) p_bar <- p_bar + coord_cartesian(ylim = ylim_range)
    p_bar
}

# -----------------------------
# 9) Labels & Coordinates
# -----------------------------
LABEL_MAP <- c(
  FCG_z = "FCG",
  tsm_z = "TSM",
  nutrient_pc1_z = "Nutrient PC1",
  physchem_pc1_z = "PhysChem PC1",
  bacteria_shannon_z = "Bacteria H'",
  protist_shannon_z = "Protist H'",
  metazoa_shannon_z = "Metazoa H'",
  fish_shannon_z = "Fish H'",
  river_width_z = "River Width",
  degree_skewness_z = "Degree Skewness",
  mean_tl_z = "Mean TL",
  mean_generality_z = "Mean Generality",
  omnivory_z = "Omnivory",
  trophic_incoherence_z = "Trophic Incoherence",
  genus_richness_z = "True Genus Richness",
  connectance_dir_z = "Connectance",
  nestedness_nodf_z = "Nestedness",
  modularity_z = "Modularity",
  niche_overlap_z = "Niche Overlap"
)

nodes_comp <- tibble(
  name = c(
    "FCG", "tsm", "nutrient_pc1", "physchem_pc1",
    "bacteria_shannon", "protist_shannon", "metazoa_shannon", "fish_shannon",
    "degree_skewness", "mean_tl", "mean_generality", "omnivory", "trophic_incoherence"
  ),
  x    = c(0, -1.7, 0, 1.7, -2.5, -0.85, 0.85, 2.5, -0.95, 0.95, -0.95, 0.95, 0),
  y    = c(4.1, 3.1, 3.1, 3.1, 2.1, 2.1, 2.1, 2.1, 1.05, 1.05, 0.0, 0.0, -1.0),
  label = c(
    "FCG", "TSM", "Nutrient PC1", "PhysChem PC1",
    "Bacteria H'", "Protist H'", "Metazoa H'", "Fish H'",
    "Degree Skewness", "Mean TL", "Mean Generality", "Omnivory", "Trophic Incoherence"
  ),
  fill_col = c(
    "#E1F5FE", rep("#FFF8E1", 3), rep("#E8F5E9", 4),
    rep("#F3E5F5", 4), "#FFEBEE"
  )
)

nodes_struct <- tibble(
  name = c(
    "FCG", "tsm", "nutrient_pc1", "physchem_pc1",
    "bacteria_shannon", "protist_shannon", "metazoa_shannon", "fish_shannon",
    "genus_richness", "connectance_dir", "nestedness_nodf", "modularity", "niche_overlap"
  ),
  x    = c(0, -1.7, 0, 1.7, -2.5, -0.85, 0.85, 2.5, -0.95, 0.95, -0.95, 0.95, 0),
  y    = c(4.1, 3.1, 3.1, 3.1, 2.1, 2.1, 2.1, 2.1, 1.05, 1.05, 0.0, 0.0, -1.0),
  label = c(
    "FCG", "TSM", "Nutrient PC1", "PhysChem PC1",
    "Bacteria H'", "Protist H'", "Metazoa H'", "Fish H'",
    "True Genus Richness", "Connectance", "Nestedness (NODF)", "Modularity", "Niche Overlap"
  ),
  fill_col = c(
    "#E1F5FE", rep("#FFF8E1", 3), rep("#E8F5E9", 4),
    rep("#F3E5F5", 4), "#FFEBEE"
  )
)

# -----------------------------
# 10) Run Composition model
# -----------------------------
cat("\n=== Part A: Composition Model ===\n")

comp_eqs <- c(
  "tsm_z ~ FCG_z + river_width_z + (1 | river)",
  "nutrient_pc1_z ~ FCG_z + river_width_z + (1 | river)",
  "physchem_pc1_z ~ FCG_z + river_width_z + (1 | river)",
  
  "bacteria_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  "protist_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  "metazoa_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  "fish_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  
  "degree_skewness_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + river_width_z + (1 | river)",
  "mean_tl_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + degree_skewness_z + river_width_z + (1 | river)",
  "mean_generality_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + degree_skewness_z + mean_tl_z + river_width_z + (1 | river)",
  "omnivory_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + degree_skewness_z + mean_tl_z + river_width_z + (1 | river)",
  "trophic_incoherence_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + degree_skewness_z + mean_tl_z + mean_generality_z + omnivory_z + river_width_z + (1 | river)"
)

fit_comp_obj <- fit_layer_model_separate(comp_eqs, dat_brms, "Composition_Model", OUT_DIR)
draws_comp   <- fit_comp_obj$draws_df

comp_res <- extract_paths_with_draws(
  draws_comp,
  from_vars = c(
    "FCG_z", "river_width_z", "tsm_z", "nutrient_pc1_z", "physchem_pc1_z",
    "bacteria_shannon_z", "protist_shannon_z", "metazoa_shannon_z", "fish_shannon_z",
    "degree_skewness_z", "mean_tl_z", "mean_generality_z", "omnivory_z"
  ),
  to_vars   = c(
    "tsm_z", "nutrient_pc1_z", "physchem_pc1_z",
    "bacteria_shannon_z", "protist_shannon_z", "metazoa_shannon_z", "fish_shannon_z",
    "degree_skewness_z", "mean_tl_z", "mean_generality_z", "omnivory_z", "trophic_incoherence_z"
  )
)
paths_comp <- comp_res$summary

effect_predictors_z <- c(
    "FCG_z",
    "tsm_z",
    "nutrient_pc1_z",
    "physchem_pc1_z",
    "bacteria_shannon_z",
    "protist_shannon_z",
    "metazoa_shannon_z",
    "fish_shannon_z"
)

predictor_levels <- c(
    "FCG",
    "TSM",
    "Nutrient PC1",
    "PhysChem PC1",
    "Bacteria H'",
    "Protist H'",
    "Metazoa H'",
    "Fish H'"
)

effects_comp_mixed <- calc_full_cascading_effects_posterior(
    edge_draws = comp_res$draws,
    preds      = effect_predictors_z,
    endpoint   = "trophic_incoherence_z",
    label_map  = LABEL_MAP
)$summary_long %>%
    mutate(
        Predictor = factor(Predictor, levels = predictor_levels),
        EffectType = factor(EffectType, levels = c("Direct", "Indirect", "Total"))
    ) %>%
    arrange(Predictor, EffectType)

# -----------------------------
# 11) Run Structure model
# -----------------------------
cat("\n=== Part B: Structure Model ===\n")

struct_eqs <- c(
  "tsm_z ~ FCG_z + river_width_z + (1 | river)",
  "nutrient_pc1_z ~ FCG_z + river_width_z + (1 | river)",
  "physchem_pc1_z ~ FCG_z + river_width_z + (1 | river)",
  
  "bacteria_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  "protist_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  "metazoa_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  "fish_shannon_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + river_width_z + (1 | river)",
  
  "genus_richness_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + river_width_z + (1 | river)",
  "connectance_dir_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + genus_richness_z + river_width_z + (1 | river)",
  "nestedness_nodf_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + genus_richness_z + connectance_dir_z + river_width_z + (1 | river)",
  "modularity_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + genus_richness_z + connectance_dir_z + river_width_z + (1 | river)",
  "niche_overlap_z ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z + bacteria_shannon_z + protist_shannon_z + metazoa_shannon_z + fish_shannon_z + genus_richness_z + connectance_dir_z + nestedness_nodf_z + modularity_z + river_width_z + (1 | river)"
)

fit_struct_obj <- fit_layer_model_separate(struct_eqs, dat_brms, "Structure_Model", OUT_DIR)
draws_struct   <- fit_struct_obj$draws_df

struct_res <- extract_paths_with_draws(
  draws_struct,
  from_vars = c(
    "FCG_z", "river_width_z", "tsm_z", "nutrient_pc1_z", "physchem_pc1_z",
    "bacteria_shannon_z", "protist_shannon_z", "metazoa_shannon_z", "fish_shannon_z",
    "genus_richness_z", "connectance_dir_z", "nestedness_nodf_z", "modularity_z"
  ),
  to_vars   = c(
    "tsm_z", "nutrient_pc1_z", "physchem_pc1_z",
    "bacteria_shannon_z", "protist_shannon_z", "metazoa_shannon_z", "fish_shannon_z",
    "genus_richness_z", "connectance_dir_z", "nestedness_nodf_z", "modularity_z", "niche_overlap_z"
  )
)
paths_struct <- struct_res$summary

effects_struct_mixed <- calc_full_cascading_effects_posterior(
    edge_draws = struct_res$draws,
    preds      = effect_predictors_z,
    endpoint   = "niche_overlap_z",
    label_map  = LABEL_MAP
)$summary_long %>%
    mutate(
        Predictor = factor(Predictor, levels = predictor_levels),
        EffectType = factor(EffectType, levels = c("Direct", "Indirect", "Total"))
    ) %>%
    arrange(Predictor, EffectType)

# -----------------------------
# 12) Save summaries
# -----------------------------
write.table(
  paths_comp,
  file = file.path(OUT_DIR, "Composition_Path_Summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  effects_comp_mixed,
  file = file.path(OUT_DIR, "Composition_Predictor_Effects.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  paths_struct,
  file = file.path(OUT_DIR, "Structure_Path_Summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
write.table(
  effects_struct_mixed,
  file = file.path(OUT_DIR, "Structure_Predictor_Effects.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# 13) Pathway diagrams (Row 1)
# -----------------------------
p_comp_path <- plot_hierarchical_sem(
  paths_comp, nodes_comp,
  title = "A. Composition SEM",
  subtitle = "Controlled model | 3 Mediators + 4 Shannon layer",
  edge_cutoff = EDGE_DISPLAY_THRESHOLD_ROW1
)

p_struct_path <- plot_hierarchical_sem(
  paths_struct, nodes_struct,
  title = "B. Structure SEM",
  subtitle = "Controlled model | 3 Mediators + 4 Shannon layer",
  edge_cutoff = EDGE_DISPLAY_THRESHOLD_ROW1
)

# -----------------------------
# 14) Calculate Global Y-Limits & Plot predictor effects (Row 2)
# -----------------------------
combined_eff_df <- bind_rows(effects_comp_mixed, effects_struct_mixed)
global_y_min <- 0
global_y_max <- 0

if(nrow(combined_eff_df) > 0){
  global_y_min <- min(c(0, combined_eff_df$est), na.rm = TRUE)
  global_y_max <- max(c(0, combined_eff_df$est), na.rm = TRUE)

  pad <- (global_y_max - global_y_min) * 0.08
  if(pad == 0) pad <- 0.1

  global_y_min <- global_y_min - pad
  global_y_max <- global_y_max + pad
}
y_limits_global <- c(global_y_min, global_y_max)

p_comp_eff <- plot_effect_sizes_fi_channels(
  effects_comp_mixed,
  title = "C. Predictor effects on Trophic Incoherence",
  ylim_range = y_limits_global
)

p_struct_eff <- plot_effect_sizes_fi_channels(
  effects_struct_mixed,
  title = "D. Predictor effects on Niche Overlap",
  ylim_range = y_limits_global
)

# -----------------------------
# 15) Final integrated layout
# -----------------------------
cat("\n=== Generating Combined Integrated Output ===\n")

row1 <- p_comp_path | p_struct_path
row2 <- p_comp_eff  | p_struct_eff

final_layout <- row1 /
  ((row2 + plot_layout(guides = "collect")) & theme(legend.position = "bottom")) +
  plot_layout(heights = c(1.45, 1.00))

out_combined_pdf <- file.path(OUT_DIR, "Fig_Integrated_SEM_Controlled_HideWidth_StrictPrior_Cut020_PredictorEffects.pdf")
save_pdf_tmp(out_combined_pdf, width = 20.0, height = 18.8, function() print(final_layout))

cat("\n✅ SUCCESS! Final controlled pipeline with FCG + FCG-channel decomposition completed.\n")
cat("Output directory: ", OUT_DIR, "\n")
