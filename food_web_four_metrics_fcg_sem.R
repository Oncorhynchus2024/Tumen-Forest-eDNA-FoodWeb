# ============================================================
# Bayesian piecewise SEM for four food-web structural endpoints
# ============================================================
#
# Purpose:
#   This script evaluates how the Forest Cover Gradient (FCG) affects
#   four inferred aquatic food-web structural metrics through three
#   environmental mediators: TSM, Nutrient PC1, and PhysChem PC1.
#
# Endpoints:
#   1. Dist_between_nodes
#   2. LinkDensity
#   3. Modularity
#   4. Clustering_und
#
# Main path structure:
#   FCG -> TSM / Nutrient PC1 / PhysChem PC1 -> food-web endpoint
#
# Required input files in the working directory:
#   abund.txt       : genus-by-site abundance or presence table
#   genus_ffg.txt   : mapping between genus names and functional feeding groups
#   ffg_edges.txt   : metaweb edge list among functional feeding groups
#   env.txt         : site-level environmental variables
#
# Important variable-name note:
#   The formal forest-cover gradient variable is FCG.
#   For backward compatibility, the script can still read an old column
#   named FI from env.txt, but it is converted internally to FCG.
#
# Main outputs:
#   - Site-level network metrics
#   - PC1 scores and loadings for Nutrient PC1 and PhysChem PC1
#   - Bayesian SEM path summaries for each endpoint
#   - Direct, indirect, and total effect summaries
#   - A combined 2-row x 4-column SEM figure
#
# ============================================================

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
set.seed(123)

# -----------------------------
# 0) Run mode and MCMC settings
# -----------------------------
# Use 'fast' for testing and 'final' for publication-quality estimates.
# The final mode increases modularity iterations, MCMC iterations, and
# Stan/NUTS tuning parameters.
RUN_MODE <- "final"   # "fast" or "final"

MODULARITY_N  <- ifelse(RUN_MODE == "final", 1000, 50)
MCMC_ITER     <- ifelse(RUN_MODE == "final", 8000, 1500)
MCMC_WARMUP   <- ifelse(RUN_MODE == "final", 4000, 700)
MCMC_CHAINS   <- 4
ADAPT_DELTA   <- ifelse(RUN_MODE == "final", 0.99, 0.95)
MAX_TREEDEPTH <- ifelse(RUN_MODE == "final", 15, 12)

PD_AUX_THRESHOLD      <- 0.975
EDGE_DISPLAY_THRESHOLD <- 0.10

# -----------------------------
# 1) Package installation and loading
# -----------------------------
# Required packages are installed automatically if they are missing.
# For a frozen reproducible environment, record package versions separately.
pkgs <- c(
  "tidyverse","janitor","igraph","ggraph","tidygraph",
  "scales","ggrepel","patchwork","brms","posterior"
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
  library(brms)
  library(posterior)
})

# -----------------------------
# 2) User settings
# -----------------------------
# Edit these file names only if your input files use different names.
ABUND_FILE <- "abund.txt"
GENUS_FFG  <- "genus_ffg.txt"
FFG_EDGES  <- "ffg_edges.txt"
ENV_FILE   <- "env.txt"

OUT_DIR <- file.path(getwd(), "SEM_4Metrics_FCG_3Mediators")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_TOTAL_ABUND   <- 1
USE_PRESENCE_ONLY <- TRUE
EDGE_DIRECTION    <- "prey_to_predator"

ALWAYS_INCLUDE_BASAL_RESOURCES <- TRUE
BASAL_ALWAYS_CAND <- c(
  "Plants","Fungi","Plankton","AquaticDetritus",
  "TerrestrialDetritus","Aquatic_Detritus",
  "Terrestrial_Detritus","Detritus","Algae"
)

# Nutrient PC1 follows the original 4.4 workflow.
NUTRIENT_VARS <- c("NH4","NO2","NO3","PO4")

# PhysChem PC1 is synchronized with the 3.27 workflow and excludes salinity.
PHYSCHEM_VARS <- c("T","DO","pH","TDS","Cond","c","Conductivity")
TSM_VAR <- "TSM"

COL_POS <- "#2C7BB6"
COL_NEG <- "#D7191C"

# -----------------------------
# 3) Helper functions
# -----------------------------
# These utilities handle file detection, data type conversion, scaling,
# transformations, PC1 calculation, model diagnostics, and plotting.
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

numify <- function(v){
  if(is.factor(v)) v <- as.character(v)
  if(is.character(v)) v <- trimws(gsub(",", "", v))
  suppressWarnings(as.numeric(v))
}

z <- function(x){
  x <- numify(x)
  mu  <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if(!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

logit01 <- function(p){
  p <- numify(p)
  p <- pmin(pmax(p, 0), 1)
  qlogis(p * 0.998 + 0.001)
}

safe_log <- function(x){
  log(pmax(numify(x), 1e-6))
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

# ---- Original 4.4-style PC1 helper used for Nutrient PC1 ----
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

# ---- Robust column matching helper synchronized with the 3.27 workflow ----
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

# ---- PhysChem PC1 helper synchronized with the 3.27 workflow ----
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

  mat <- as.data.frame(lapply(df0[, selected_cols, drop = FALSE], numify))

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
    v <- numify(mat[[nm]])
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

  # For PhysChem PC1, align PC1 direction to temperature when available; otherwise align to row sums.
  if(is.null(align_to)){
    anchor <- rowSums(mat, na.rm = TRUE)
  } else {
    anchor <- numify(align_to)
    med <- median(anchor[is.finite(anchor)], na.rm = TRUE)
    if(!is.finite(med)) med <- 0
    anchor[!is.finite(anchor)] <- med
  }

  flip <- 1
  cc <- suppressWarnings(cor(pc1, anchor, use = "pairwise.complete.obs"))
  if(is.finite(cc) && cc < 0) flip <- -1

  list(
    scores   = pc1 * flip,
    loadings = pc$rotation[, 1] * flip,
    used_cols = colnames(mat),
    var_prop = (pc$sdev^2) / sum(pc$sdev^2)
  )
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

  tibble(
    est   = median(x),
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    pd    = pd,
    sig   = cr_interval_sig(ci_lo, ci_hi),
    pd_sig = is.finite(pd) && pd >= PD_AUX_THRESHOLD,
    stars = stars_from_cri(x)
  )
}

prepare_complete_case_data <- function(formulas_list, dat){
  vars_needed <- unique(unlist(lapply(formulas_list, function(f) all.vars(as.formula(f)))))
  dat_use <- dat %>%
    filter(if_all(all_of(vars_needed), ~ is.finite(.x)))
  list(data = dat_use, vars_needed = vars_needed)
}

save_brms_diagnostics <- function(fit, fit_name, out_dir){
  diag_dir <- file.path(out_dir, paste0("Diagnostics_", fit_name))
  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

  writeLines(
    capture.output(summary(fit)),
    con = file.path(diag_dir, paste0(fit_name, "_summary.txt"))
  )

  draws_mat <- posterior::as_draws_matrix(fit)
  par_keep <- colnames(draws_mat)
  par_keep <- par_keep[grepl("^b_|^sigma", par_keep)]

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
}

fit_layer_model <- function(formulas_list, dat, fit_name, out_dir){
  prep <- prepare_complete_case_data(formulas_list, dat)
  dat_use <- prep$data

  if(nrow(dat_use) < 20){
    stop(fit_name, ": complete-case sample size too small (n < 20).")
  }

  bform <- bf(as.formula(formulas_list[[1]]))
  for(i in 2:length(formulas_list)){
    bform <- bform + bf(as.formula(formulas_list[[i]]))
  }
  bform <- bform + set_rescor(FALSE)

  fit <- brm(
    formula = bform,
    data    = dat_use,
    family  = gaussian(),
    chains  = MCMC_CHAINS,
    cores   = max(1, parallel::detectCores() - 1),
    iter    = MCMC_ITER,
    warmup  = MCMC_WARMUP,
    backend = "rstan",
    save_pars = save_pars(all = TRUE),
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
    silent  = 2,
    refresh = 0,
    seed    = 123
  )

  save_brms_diagnostics(fit, fit_name, out_dir)

  list(
    fit = fit,
    data = dat_use,
    vars_used = prep$vars_needed
  )
}

get_b_draws <- function(draws, response, predictor){
  nm <- paste0("b_", gsub("_", "", response), "_", predictor)
  if(!(nm %in% names(draws))) stop("Missing posterior column: ", nm)
  as.numeric(draws[[nm]])
}

extract_paths_with_draws <- function(draws, from_vars, to_vars) {
  paths_list <- list()
  for(to in to_vars) {
    for(from in from_vars) {
      if(from != to) {
        nm <- paste0("b_", gsub("_", "", to), "_", from)
        if(nm %in% names(draws)) {
          vec <- as.numeric(draws[[nm]])
          sm <- summarize_posterior_vec(vec)
          paths_list[[length(paths_list) + 1]] <- tibble(
            from = from,
            to = to,
            est = sm$est,
            ci_lo = sm$ci_lo,
            ci_hi = sm$ci_hi,
            pd = sm$pd,
            sig = sm$sig,
            pd_sig = sm$pd_sig,
            stars = sm$stars
          )
        }
      }
    }
  }
  list(summary = bind_rows(paths_list))
}

calc_endpoint_effects_posterior <- function(draws, endpoint_z){

  b_fi_to_tsm <- get_b_draws(draws, "tsm_z", "FCG_z")
  b_fi_to_nut <- get_b_draws(draws, "nutrient_pc1_z", "FCG_z")
  b_fi_to_phy <- get_b_draws(draws, "physchem_pc1_z", "FCG_z")

  b_end_fi  <- get_b_draws(draws, endpoint_z, "FCG_z")
  b_end_tsm <- get_b_draws(draws, endpoint_z, "tsm_z")
  b_end_nut <- get_b_draws(draws, endpoint_z, "nutrient_pc1_z")
  b_end_phy <- get_b_draws(draws, endpoint_z, "physchem_pc1_z")

  n_draws <- length(b_end_fi)
  zero_vec <- rep(0, n_draws)

  eff_draws <- list(
    "FCG" = list(
      Direct   = b_end_fi,
      Indirect = b_fi_to_tsm * b_end_tsm + b_fi_to_nut * b_end_nut + b_fi_to_phy * b_end_phy,
      Total    = b_end_fi + b_fi_to_tsm * b_end_tsm + b_fi_to_nut * b_end_nut + b_fi_to_phy * b_end_phy
    ),
    "TSM" = list(
      Direct   = b_end_tsm,
      Indirect = zero_vec,
      Total    = b_end_tsm
    ),
    "Nutrient PC1" = list(
      Direct   = b_end_nut,
      Indirect = zero_vec,
      Total    = b_end_nut
    ),
    "PhysChem PC1" = list(
      Direct   = b_end_phy,
      Indirect = zero_vec,
      Total    = b_end_phy
    )
  )

  sum_list <- list()
  for(pred in names(eff_draws)){
    for(tp in names(eff_draws[[pred]])){
      sm <- summarize_posterior_vec(eff_draws[[pred]][[tp]])
      sum_list[[length(sum_list) + 1]] <- tibble(
        Predictor = pred,
        EffectType = tp,
        est = sm$est,
        ci_lo = sm$ci_lo,
        ci_hi = sm$ci_hi,
        pd = sm$pd,
        sig = sm$sig,
        pd_sig = sm$pd_sig,
        stars = sm$stars
      )
    }
  }

  list(
    summary_long = bind_rows(sum_list),
    draws = eff_draws
  )
}

plot_effect_sizes_posterior <- function(eff_long, title, ylim_range = NULL){

  eff_plot <- eff_long %>%
    mutate(
      EffectType = factor(EffectType, levels = c("Direct", "Indirect", "Total")),
      Predictor  = factor(Predictor, levels = c("FCG", "TSM", "Nutrient PC1", "PhysChem PC1"))
    )

  dodge <- position_dodge(width = 0.75)

  p_bar <- ggplot(eff_plot, aes(x = Predictor, y = est, fill = EffectType)) +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    geom_col(position = dodge, width = 0.68, color = "grey25", linewidth = 0.35) +
    scale_fill_manual(values = c(
      "Direct"   = "#9ECAE1",
      "Indirect" = "#A1D99B",
      "Total"    = "#FDD0A2"
    )) +
    theme_bw(base_size = 18) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      axis.text.x = element_text(size = 15, face = "bold", color = "black", angle = 15, hjust = 1),
      axis.text.y = element_text(size = 14, color = "black"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 16),
      plot.title = element_text(face = "bold", size = 22, hjust = 0.5, margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      y = "Posterior effect size\n(median)",
      fill = "Effect type",
      title = title
    )

  if(!is.null(ylim_range) && length(ylim_range) == 2 && all(is.finite(ylim_range))){
    p_bar <- p_bar + coord_cartesian(ylim = ylim_range)
  }

  p_bar
}

plot_single_endpoint_sem <- function(paths, nodes_df, title, subtitle = NULL, edge_cutoff = 0.10){

  paths_clean <- paths %>%
    mutate(
      from = gsub("_z$", "", from),
      to   = gsub("_z$", "", to)
    ) %>%
    filter(is.finite(est), abs(est) >= edge_cutoff)

  edges <- paths_clean %>%
    filter(from %in% nodes_df$name, to %in% nodes_df$name) %>%
    mutate(
      col       = ifelse(sig, ifelse(est >= 0, COL_POS, COL_NEG), "grey75"),
      lty       = ifelse(sig, "solid", "dashed"),
      alpha_v   = ifelse(sig, 1.0, 0.6),
      w         = ifelse(sig, 0.6 + abs(est) * 2.0, 0.45),
      elab      = sprintf("%.2f%s", est, stars),
      font_face = ifelse(sig, "bold", "plain"),
      text_col  = ifelse(sig, "black", "grey40")
    )

  edges_mid <- edges %>%
    left_join(nodes_df, by = c("from" = "name")) %>% rename(x1 = x, y1 = y) %>%
    left_join(nodes_df, by = c("to" = "name"))   %>% rename(x2 = x, y2 = y) %>%
    mutate(
      xm = 0.42 * x1 + 0.58 * x2,
      ym = 0.42 * y1 + 0.58 * y2
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
      arrow = grid::arrow(length = grid::unit(3.0, "mm"), type = "closed"),
      end_cap = ggraph::circle(10.0, "mm"),
      start_cap = ggraph::circle(10.0, "mm"),
      lineend = "round",
      show.legend = FALSE
    ) +
    geom_node_label(
      aes(label = label, fill = fill_col),
      size = 6.0,
      fontface = "bold",
      color = "grey20",
      label.padding = grid::unit(0.40, "lines"),
      label.r = grid::unit(0.22, "lines"),
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = edges_mid,
      aes(x = xm, y = ym, label = elab, color = text_col, fontface = font_face),
      inherit.aes = FALSE,
      size = ifelse(edges_mid$sig, 5.5, 4.5),
      box.padding = 0.10,
      point.padding = 0.03,
      min.segment.length = 0,
      seed = 123,
      show.legend = FALSE
    ) +
    scale_fill_identity() +
    scale_color_identity() +
    scale_edge_colour_identity() +
    scale_edge_alpha_identity() +
    scale_edge_width_identity() +
    scale_edge_linetype_identity() +
    coord_cartesian(
      xlim = range(nodes_df$x) + c(-0.65, 0.65),
      ylim = range(nodes_df$y) + c(-0.35, 0.35),
      clip = "off"
    ) +
    theme_void() +
    labs(title = title, subtitle = subtitle) +
    theme(
      plot.title = element_text(face = "bold", size = 24, hjust = 0.5, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 15, hjust = 0.5, color = "grey25", margin = margin(b = 10)),
      plot.margin = margin(10, 10, 10, 10)
    )

  p
}

run_one_endpoint_sem <- function(dat, endpoint_z, endpoint_label, out_dir){

  eqs <- c(
    "tsm_z ~ FCG_z",
    "nutrient_pc1_z ~ FCG_z",
    "physchem_pc1_z ~ FCG_z",
    paste0(endpoint_z, " ~ FCG_z + tsm_z + nutrient_pc1_z + physchem_pc1_z")
  )

  fit_name <- paste0("SEM_", endpoint_label)

  fit_obj <- fit_layer_model(eqs, dat, fit_name, out_dir)
  draws   <- as_draws_df(fit_obj$fit)

  res <- extract_paths_with_draws(
    draws,
    from_vars = c("FCG_z","tsm_z","nutrient_pc1_z","physchem_pc1_z"),
    to_vars   = c("tsm_z","nutrient_pc1_z","physchem_pc1_z", endpoint_z)
  )

  eff <- calc_endpoint_effects_posterior(draws, endpoint_z)

  write.table(
    res$summary,
    file = file.path(out_dir, paste0(fit_name, "_Path_Summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  write.table(
    eff$summary_long,
    file = file.path(out_dir, paste0(fit_name, "_Effects_Summary.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  nodes_df <- tibble(
    name = c("FCG", "tsm", "nutrient_pc1", "physchem_pc1", gsub("_z$", "", endpoint_z)),
    x    = c(0, -1.4, 0, 1.4, 0),
    y    = c(3.0, 2.0, 2.0, 2.0, 0.8),
    label = c("FCG", "TSM", "Nutrient PC1", "PhysChem PC1", endpoint_label),
    fill_col = c("#E1F5FE", "#FFF8E1", "#FFF8E1", "#FFF8E1", "#E8F5E9")
  )

  p_path <- plot_single_endpoint_sem(
    paths = res$summary,
    nodes_df = nodes_df,
    title = endpoint_label,
    subtitle = "FCG -> mediators -> endpoint",
    edge_cutoff = EDGE_DISPLAY_THRESHOLD
  )

  p_eff <- plot_effect_sizes_posterior(
    eff_long = eff$summary_long,
    title = paste0(endpoint_label, " effects")
  )

  list(
    plot = p_path,
    effect_plot = p_eff,
    fit = fit_obj$fit,
    paths = res$summary,
    effects = eff$summary_long
  )
}

# -----------------------------
# 4) Load ENV and synchronize PhysChem PC1 to 3.27
# -----------------------------
env_raw <- read_tsv_safe(find_file(ENV_FILE)) %>%
  janitor::clean_names()

if(!("site" %in% names(env_raw))) stop("env.txt must have a 'site' column.")
if(!("river" %in% names(env_raw))) env_raw$river <- "all_sites"

# Convert environmental variables to numeric values.
for(nm in names(env_raw)){
  if(!nm %in% c("site","river")) env_raw[[nm]] <- numify(env_raw[[nm]])
}

# Relaxed column matching synchronized with the 3.27 workflow.
fcg_col  <- find_col_relaxed(names(env_raw), c("FCG", "FI"))
tsm_col <- find_col_relaxed(names(env_raw), c("TSM"))
t_col   <- find_col_relaxed(names(env_raw), c("T"))

if(is.na(fcg_col))  stop("env.txt must contain either an FCG or FI column.")
if(is.na(tsm_col)) stop("env.txt must contain a TSM column.")

# Nutrient PC1 follows the original 4.4 workflow.
env_raw <- do_pc1(env_raw, NUTRIENT_VARS, "nutrient_pc1")

# PhysChem PC1 is fully synchronized with the 3.27 workflow.
physchem_cols <- c(
  find_col_relaxed(names(env_raw), c("T")),
  find_col_relaxed(names(env_raw), c("DO")),
  find_col_relaxed(names(env_raw), c("pH")),
  find_col_relaxed(names(env_raw), c("TDS")),
  find_col_relaxed(names(env_raw), c("Cond", "c", "Conductivity"))
)

phys_res <- compute_pc1_generic(
  df0 = env_raw,
  selected_cols = physchem_cols,
  align_to = if(!is.na(t_col)) env_raw[[t_col]] else NULL,
  miss_keep = 0.8
)

env_raw$physchem_pc1 <- phys_res$scores

# Export PC1 scores and loading information.
env_idx_export <- env_raw %>%
  dplyr::select(site, river, dplyr::all_of(fcg_col), dplyr::all_of(tsm_col), nutrient_pc1, physchem_pc1)

write.table(
  env_idx_export,
  file = file.path(OUT_DIR, "env_with_PC1_indices.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

if(!is.null(phys_res$loadings)){
  phys_load <- tibble::tibble(
    variable = names(phys_res$loadings),
    loading_PC1 = as.numeric(phys_res$loadings),
    contribution = (as.numeric(phys_res$loadings)^2) / sum((as.numeric(phys_res$loadings)^2))
  ) %>%
    dplyr::arrange(dplyr::desc(abs(loading_PC1)))

  write.table(
    phys_load,
    file = file.path(OUT_DIR, "PhysChem_PC1_Loadings_synced_to_3.27.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

write.table(
  data.frame(
    variable = c("FCG", "TSM", "nutrient_pc1", "physchem_pc1"),
    sd = c(
      sd(env_raw[[fcg_col]], na.rm = TRUE),
      sd(env_raw[[tsm_col]], na.rm = TRUE),
      sd(env_raw$nutrient_pc1, na.rm = TRUE),
      sd(env_raw$physchem_pc1, na.rm = TRUE)
    )
  ),
  file = file.path(OUT_DIR, "ENV_variation_check.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Convert the environmental table to the downstream structure used by this workflow.
env <- env_raw
if("site" %in% names(env)) names(env)[names(env) == "site"] <- "Site"
if("fcg" %in% names(env)) names(env)[names(env) == "fcg"] <- "FCG"
if("fi"  %in% names(env)) names(env)[names(env) == "fi"]  <- "FCG"
env$Site <- as.character(env$Site)

# -----------------------------
# 5) Load genus / FFG / abundance
# -----------------------------
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

site_ffg <- abund_long %>%
  inner_join(gffg, by = "Genus") %>%
  distinct(Site, FFG)

# -----------------------------
# 6) Calculate 4 target network metrics
# -----------------------------
calc_topology_target_metrics <- function(g_dir, mod_n){

  S <- vcount(g_dir)
  L <- ecount(g_dir)

  g_und <- as_undirected(g_dir, mode = "collapse")

  LinkDensity <- if(S > 0) L / S else NA_real_

  Dist_between_nodes <- if(S >= 2 && ecount(g_und) > 0){
    suppressWarnings(igraph::mean_distance(g_und, directed = FALSE, unconnected = TRUE))
  } else {
    NA_real_
  }

  Clustering_und <- if(S >= 3 && ecount(g_und) > 1){
    suppressWarnings(igraph::transitivity(g_und, type = "globalundirected"))
  } else {
    NA_real_
  }
  if(!is.finite(Clustering_und)) Clustering_und <- NA_real_

  Modularity <- if(S >= 3 && ecount(g_und) > 0){
    mean(sapply(1:mod_n, function(x){
      set.seed(2025 + x)
      igraph::modularity(igraph::cluster_louvain(igraph::permute(g_und, sample.int(S))))
    }), na.rm = TRUE)
  } else {
    NA_real_
  }

  tibble(
    Dist_between_nodes = Dist_between_nodes,
    LinkDensity = LinkDensity,
    Modularity = Modularity,
    Clustering_und = Clustering_und
  )
}

cat("\nCalculating target network metrics...\n")
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
    bind_cols(calc_topology_target_metrics(g, MODULARITY_N))
}))
close(pb)

metrics_df <- metrics_df %>%
  left_join(env, by = "Site")

write.table(
  metrics_df,
  file = file.path(OUT_DIR, "Network_4Metrics_By_Site.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# -----------------------------
# 7) Prepare SEM dataset
# -----------------------------
df <- metrics_df %>%
  filter(is.finite(FCG)) %>%
  clean_names() %>%
  rename(FCG = fcg)

# Do not recompute PhysChem PC1 here; use the synchronized PhysChem PC1 computed above.
if(!"nutrient_pc1" %in% names(df)) stop("Column 'nutrient_pc1' not found.")
if(!"physchem_pc1" %in% names(df)) stop("Column 'physchem_pc1' not found.")

tsm_col_clean <- janitor::make_clean_names(TSM_VAR)
df$tsm <- if(tsm_col_clean %in% names(df)) as.numeric(df[[tsm_col_clean]]) else NA_real_

df$FCG_z           <- z(df$FCG)
df$tsm_z          <- z(df$tsm)
df$nutrient_pc1_z <- z(df$nutrient_pc1)
df$physchem_pc1_z <- z(df$physchem_pc1)

# endpoint transforms
if(!"dist_between_nodes" %in% names(df)) stop("Column 'dist_between_nodes' not found after clean_names().")
if(!"link_density" %in% names(df)) stop("Column 'link_density' not found after clean_names().")
if(!"modularity" %in% names(df)) stop("Column 'modularity' not found after clean_names().")
if(!"clustering_und" %in% names(df)) stop("Column 'clustering_und' not found after clean_names().")

df$dist_between_nodes_z <- z(safe_log(df$dist_between_nodes))
df$linkdensity_z        <- z(safe_log(df$link_density))
df$modularity_z         <- z(logit01(df$modularity))
df$clustering_und_z     <- z(logit01(df$clustering_und))

dat_brms <- df %>%
  filter(is.finite(FCG_z))

# -----------------------------
# 8) Run four SEMs
# -----------------------------
cat("\n=== Running four SEMs ===\n")

res_dist <- run_one_endpoint_sem(
  dat = dat_brms,
  endpoint_z = "dist_between_nodes_z",
  endpoint_label = "Dist_between_nodes",
  out_dir = OUT_DIR
)

res_link <- run_one_endpoint_sem(
  dat = dat_brms,
  endpoint_z = "linkdensity_z",
  endpoint_label = "LinkDensity",
  out_dir = OUT_DIR
)

res_mod <- run_one_endpoint_sem(
  dat = dat_brms,
  endpoint_z = "modularity_z",
  endpoint_label = "Modularity",
  out_dir = OUT_DIR
)

res_clu <- run_one_endpoint_sem(
  dat = dat_brms,
  endpoint_z = "clustering_und_z",
  endpoint_label = "Clustering_und",
  out_dir = OUT_DIR
)

# -----------------------------
# 9) Combine plots: 2 rows x 4 columns
# -----------------------------
cat("\n=== Combining final figure ===\n")

combined_eff_df <- bind_rows(
  res_dist$effects,
  res_link$effects,
  res_mod$effects,
  res_clu$effects
)

global_y_min <- min(c(0, combined_eff_df$est), na.rm = TRUE)
global_y_max <- max(c(0, combined_eff_df$est), na.rm = TRUE)
pad <- (global_y_max - global_y_min) * 0.08
if(!is.finite(pad) || pad == 0) pad <- 0.1
y_limits_global <- c(global_y_min - pad, global_y_max + pad)

res_dist$effect_plot <- plot_effect_sizes_posterior(res_dist$effects, "Dist_between_nodes effects", y_limits_global)
res_link$effect_plot <- plot_effect_sizes_posterior(res_link$effects, "LinkDensity effects", y_limits_global)
res_mod$effect_plot  <- plot_effect_sizes_posterior(res_mod$effects,  "Modularity effects", y_limits_global)
res_clu$effect_plot  <- plot_effect_sizes_posterior(res_clu$effects,  "Clustering_und effects", y_limits_global)

row1 <- res_dist$plot | res_link$plot | res_mod$plot | res_clu$plot
row2 <- res_dist$effect_plot | res_link$effect_plot | res_mod$effect_plot | res_clu$effect_plot

final_plot <- row1 / ((row2 + plot_layout(guides = "collect")) & theme(legend.position = "bottom")) +
  plot_layout(heights = c(1.45, 1.00))

out_pdf <- file.path(OUT_DIR, "Fig_SEM_4Metrics_2rows4col_PhysChemSyncedTo3.27.pdf")
save_pdf_tmp(out_pdf, width = 26, height = 14.2, function() print(final_plot))

ggsave(
  filename = file.path(OUT_DIR, "Fig_SEM_4Metrics_2rows4col_PhysChemSyncedTo3.27.png"),
  plot = final_plot, width = 26, height = 14.2, dpi = 500, bg = "white"
)

cat("\n✅ SUCCESS! Finished 4-endpoint SEM pipeline.\n")
cat("Output directory: ", OUT_DIR, "\n")
cat("Main figure: ", out_pdf, "\n")
cat("\nNutrient PC1 (4.4 logic) variables:", paste(janitor::make_clean_names(NUTRIENT_VARS), collapse = ", "), "\n")
cat("PhysChem PC1 (synced to 3.27) used columns:", paste(phys_res$used_cols, collapse = ", "), "\n")
cat("PhysChem PC1 anchor:", if(!is.na(t_col)) t_col else "rowSums(mat)", "\n")
cat("PhysChem PC1 miss_keep:", 0.8, "\n")
cat("PhysChem PC1 includes Salinity?: NO\n")