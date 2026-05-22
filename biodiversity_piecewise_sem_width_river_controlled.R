# ============================================================
# Bayesian piecewise SEM for biodiversity responses across four groups
# ============================================================
#
# Purpose
# -------
# This script fits integrated Bayesian piecewise structural equation models
# for biodiversity responses across four taxonomic groups: bacteria, protists,
# metazoans, and fishes. The model evaluates how the forest gradient index
# (FI) affects Shannon diversity, evenness, and site-level beta diversity.
#
# Key statistical controls
# ------------------------
# 1. River identity is included as a random intercept when more than one river
#    is available after site matching.
# 2. River width and the FI-by-width interaction are always included as
#    covariates, so the reported FI effects are width-controlled.
# 3. Potential environmental mediators include Nutrient PC1, PhysChem PC1,
#    and TSM when they have sufficient variation.
#
# Beta-diversity definition
# -------------------------
# Within each river, all possible three-site combinations are generated. For
# each focal site, beta diversity is calculated as the average Bray-Curtis
# distance from that site to the other two sites within each triplet, then
# averaged across all triplets containing the focal site.
#
# Required input files
# --------------------
# - env.txt: site-level environmental data, including site, river, FI, TSM,
#   river_width, and the variables used to derive Nutrient PC1 and PhysChem PC1.
# - spe_bacteria.txt
# - spe_protist.txt
# - spe_Metazoa.txt
# - spe_fish.txt
#
# Main outputs
# ------------
# Outputs are written to SEM_Biodiv_Integrated_4Groups_Parallel_WidthControlled/.
# The folder contains model diagnostics, path tables, cascading effect tables,
# beta-triplet counts, and integrated SEM figures in PDF and PNG formats.
# ============================================================

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
set.seed(123)

# -----------------------------
# 0) RUN MODE & MCMC SETTINGS
# -----------------------------
RUN_MODE <- "final"   # "fast" or "final"

MCMC_ITER    <- ifelse(RUN_MODE == "final", 5000, 1000)
MCMC_WARMUP  <- ifelse(RUN_MODE == "final", 2000, 500)
MCMC_CHAINS  <- 4

ADAPT_DELTA   <- ifelse(RUN_MODE == "final", 0.99, 0.95)
MAX_TREEDEPTH <- ifelse(RUN_MODE == "final", 15, 12)

# Force a river-level random intercept to account for among-river differences.
USE_RIVER_RANDOM_EFFECT <- TRUE
CONTROL_READS_BIODIV    <- FALSE

# -----------------------------
# 1) Packages
# -----------------------------
pkgs <- c(
  "vegan", "tidyverse", "janitor", "scales", "ggrepel",
  "patchwork", "brms", "posterior", "bayesplot"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install) > 0){
  install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(vegan)
  library(tidyverse)
  library(janitor)
  library(scales)
  library(ggrepel)
  library(patchwork)
  library(brms)
  library(posterior)
  library(bayesplot)
})

# -----------------------------
# 2) User Settings
# -----------------------------
SPE_FILES <- list(
  "Bacteria" = "spe_bacteria.txt",
  "Protist"  = "spe_protist.txt",
  "Metazoa"  = "spe_Metazoa.txt",
  "Fish"     = "spe_fish.txt"
)

GROUP_ORDER <- c("Bacteria", "Protist", "Metazoa", "Fish")

GROUP_LABELS <- c(
  "Bacteria" = "Bacteria",
  "Protist"  = "Protist",
  "Metazoa"  = "Metazoa",
  "Fish"     = "Fish"
)

ENV_FILE <- "env.txt"

OUT_DIR <- file.path(getwd(), "SEM_Biodiv_Integrated_4Groups_Parallel_WidthControlled")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

MIN_TOTAL_READS_SITE <- 1
COL_POS <- "#2C7BB6"
COL_NEG <- "#D7191C"

# -----------------------------
# 3) Helpers
# -----------------------------
read_tsv_safe <- function(path){
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
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
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if(!is.finite(sdv) || sdv == 0) return(rep(0, length(x)))
  (x - mu) / sdv
}

logit01 <- function(p){
  p <- numify(p)
  qlogis(pmin(pmax(p, 0), 1) * 0.998 + 0.001)
}

save_pdf_tmp <- function(path, width, height, plot_fun){
  tmp <- paste0(path, ".tmp.pdf")
  if(isTRUE(capabilities("cairo"))){
    grDevices::cairo_pdf(tmp, width = width, height = height, onefile = TRUE)
  } else {
    grDevices::pdf(tmp, width = width, height = height, onefile = TRUE)
  }
  plot_fun()
  grDevices::dev.off()
  file.rename(tmp, path)
}

cr_interval_sig <- function(ci_lo, ci_hi){
  is.finite(ci_lo) && is.finite(ci_hi) && (ci_lo > 0 || ci_hi < 0)
}

stars_from_cri <- function(x){
  x <- x[is.finite(x)]
  if(length(x) < 10) return("")
  q <- quantile(
    x,
    probs = c(0.0005, 0.005, 0.025, 0.975, 0.995, 0.9995),
    na.rm = TRUE,
    names = FALSE
  )
  if(q[1] > 0 || q[6] < 0) return("***")
  if(q[2] > 0 || q[5] < 0) return("**")
  if(q[3] > 0 || q[4] < 0) return("*")
  ""
}

summarize_posterior_vec <- function(x){
  x <- as.numeric(x[is.finite(x)])
  if(length(x) < 2){
    return(tibble::tibble(
      est = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
      pd = NA_real_, sig = FALSE, stars = ""
    ))
  }
  ci_lo <- quantile(x, 0.025)
  ci_hi <- quantile(x, 0.975)
  pd <- max(mean(x > 0), mean(x < 0))
  tibble::tibble(
    est   = median(x),
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    pd    = pd,
    sig   = cr_interval_sig(ci_lo, ci_hi),
    stars = stars_from_cri(x)
  )
}

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

has_variation <- function(x){
  x <- x[is.finite(x)]
  if(length(x) < 2) return(FALSE)
  length(unique(x)) > 1 && isTRUE(sd(x) > 0)
}

pick_available_mediators <- function(dat){
  meds <- c("nutrient_pc1_z", "physchem_pc1_z", "tsm_z")
  meds[sapply(meds, function(v) v %in% names(dat) && has_variation(dat[[v]]))]
}

# Always include river width and the FI-by-width interaction in the base model.
build_rhs_base <- function(dat){
  "FI_z + river_width_z + int_FI_river_width_z"
}

paste_terms <- function(...){
  x <- c(...)
  x <- x[!is.na(x) & nzchar(x)]
  paste(x, collapse = " + ")
}

make_placeholder_plot <- function(title, subtitle = "No valid result"){
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0.5, y = 0.58, label = title, fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = 0.5, y = 0.42, label = subtitle, size = 4, color = "grey40") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

make_row_title <- function(txt){
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0.5, label = txt, hjust = 0, fontface = "bold", size = 6) +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
    ggplot2::theme_void()
}

# -----------------------------
# 4) Diagnostic Saver
# -----------------------------
save_brms_diagnostics <- function(fit, fit_name, out_dir){
  diag_dir <- file.path(out_dir, paste0("Diagnostics_", fit_name))
  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

  writeLines(
    capture.output(summary(fit)),
    con = file.path(diag_dir, paste0(fit_name, "_summary.txt"))
  )

  draws_mat <- posterior::as_draws_matrix(fit)
  par_keep <- colnames(draws_mat)[grepl("^b_|^sd_|^sigma", colnames(draws_mat))]

  if(length(par_keep) > 0){
    diag_tbl <- posterior::summarise_draws(
      draws_mat[, par_keep, drop = FALSE],
      "mean", "sd", "median", "mad", "rhat", "ess_bulk", "ess_tail"
    )

    write.table(
      as.data.frame(diag_tbl),
      file = file.path(diag_dir, paste0(fit_name, "_diagnostics.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    p_trace <- bayesplot::mcmc_trace(
      as.array(fit),
      pars = par_keep[1:min(12, length(par_keep))]
    ) + ggplot2::ggtitle(paste0("Trace plot: ", fit_name))

    save_pdf_tmp(
      file.path(diag_dir, paste0(fit_name, "_traceplot.pdf")),
      width = 12, height = 10,
      plot_fun = function() print(p_trace)
    )
  }
}

# -----------------------------
# 5) Biodiversity Metrics
# -----------------------------
# Core update:
# Within each river, beta diversity is computed from all possible three-site combinations.
# For each focal site, all triplets that contain the focal site are identified.
# The focal site's mean Bray-Curtis distance to the other two sites is averaged across all relevant triplets.

compute_beta_triplet_allcombn_within_river <- function(X_rel, site_river){
  site_names <- rownames(X_rel)

  beta_out <- setNames(rep(NA_real_, length(site_names)), site_names)
  triplet_n <- setNames(rep(0L, length(site_names)), site_names)

  rivers <- unique(site_river[!is.na(site_river)])

  for(rv in rivers){
    rv_sites <- site_names[site_river[site_names] == rv]

    # A river must contain at least three sites to support triplet-based beta diversity.
    if(length(rv_sites) < 3) next

    comb_list <- combn(rv_sites, 3, simplify = FALSE)

    site_vals <- setNames(vector("list", length(rv_sites)), rv_sites)
    for(i in seq_along(site_vals)) site_vals[[i]] <- numeric(0)

    for(cm in comb_list){
      sub_mat <- X_rel[cm, , drop = FALSE]
      d_sub <- as.matrix(vegan::vegdist(sub_mat, method = "bray"))
      diag(d_sub) <- NA_real_

      # For each site in the triplet, calculate its mean distance to the other two sites.
      local_beta <- rowMeans(d_sub, na.rm = TRUE)

      for(s in names(local_beta)){
        site_vals[[s]] <- c(site_vals[[s]], unname(local_beta[s]))
      }
    }

    for(s in rv_sites){
      vals <- site_vals[[s]]
      if(length(vals) > 0){
        beta_out[s] <- mean(vals, na.rm = TRUE)
        triplet_n[s] <- length(vals)
      }
    }
  }

  list(beta_mean = beta_out, beta_n_triplets = triplet_n)
}

compute_biodiv_metrics <- function(spe_file, env_info){
  spe <- read_tsv_safe(spe_file)

  rownames(spe) <- spe[[1]]
  spe <- spe[, -1, drop = FALSE]

  tax_idx <- which(tolower(colnames(spe)) %in% c("tax", "taxonomy"))
  if(length(tax_idx) > 0) spe <- spe[, -tax_idx[1], drop = FALSE]

  env_info <- env_info %>%
    dplyr::select(site, river) %>%
    dplyr::distinct(site, .keep_all = TRUE)

  site_cols <- intersect(colnames(spe), env_info$site)
  if(length(site_cols) < 2) stop(paste("Matched sites < 2 in", spe_file))

  X <- t(as.matrix(spe[, site_cols, drop = FALSE]))
  storage.mode(X) <- "numeric"
  rownames(X) <- site_cols

  X <- X[rowSums(X, na.rm = TRUE) >= MIN_TOTAL_READS_SITE, , drop = FALSE]
  X <- X[, colSums(X, na.rm = TRUE) > 0, drop = FALSE]

  if(nrow(X) < 3) stop(paste("Too few valid sites in", spe_file))

  env_sub <- env_info %>%
    dplyr::filter(site %in% rownames(X))

  site_river <- env_sub$river[match(rownames(X), env_sub$site)]
  names(site_river) <- rownames(X)

  if(any(is.na(site_river))){
    stop(paste("Some matched sites in", spe_file, "do not have river information in env.txt"))
  }

  shannon  <- vegan::diversity(X, index = "shannon")
  richness <- vegan::specnumber(X)
  evenness <- ifelse(richness >= 2, shannon / log(richness), NA_real_)

  X_rel <- vegan::decostand(X, method = "total", MARGIN = 1)

  # Updated beta-diversity algorithm based on all within-river site triplets.
  beta_res <- compute_beta_triplet_allcombn_within_river(
    X_rel = X_rel,
    site_river = site_river
  )

  tibble::tibble(
    site            = rownames(X),
    shannon         = as.numeric(shannon),
    evenness        = as.numeric(evenness),
    beta_mean       = as.numeric(beta_res$beta_mean[rownames(X)]),
    beta_n_triplets = as.integer(beta_res$beta_n_triplets[rownames(X)]),
    total_reads     = as.numeric(rowSums(X, na.rm = TRUE))
  )
}

# -----------------------------
# 6) Posterior Path Extraction + Cascading Effects
# -----------------------------
LABEL_MAP <- c(
  FI_z                 = "FI",
  river_width_z        = "Width",
  int_FI_river_width_z = "FI × Width",
  nutrient_pc1_z       = "Nutrient PC1",
  physchem_pc1_z       = "PhysChem PC1",
  tsm_z                = "TSM",
  shannon_z            = "Shannon",
  evenness_z           = "Evenness",
  beta_z               = "Beta"
)

nice_label <- function(x){
  y <- unname(LABEL_MAP[x])
  ifelse(is.na(y), x, y)
}

extract_paths_table <- function(fit, all_vars){
  draws <- posterior::as_draws_df(fit)
  out <- list()

  for(to in all_vars){
    for(from in all_vars){
      col_name <- paste0("b_", gsub("_", "", to), "_", from)
      if(col_name %in% names(draws)){
        x <- as.numeric(draws[[col_name]])
        sm <- summarize_posterior_vec(x)
        out[[length(out) + 1]] <- tibble::tibble(
          from  = from,
          to    = to,
          est   = sm$est,
          ci_lo = sm$ci_lo,
          ci_hi = sm$ci_hi,
          pd    = sm$pd,
          sig   = sm$sig,
          stars = sm$stars
        )
      }
    }
  }

  if(length(out) == 0) stop("No posterior paths found. Check brms parameter names.")
  dplyr::bind_rows(out)
}

extract_cascading_effects_multi <- function(fit, preds, endpoints, all_vars){
  draws <- posterior::as_draws_df(fit)

  edge_draws <- list()
  for(to in all_vars){
    for(from in all_vars){
      col_name <- paste0("b_", gsub("_", "", to), "_", from)
      if(col_name %in% names(draws)){
        edge_draws[[paste0(from, "->", to)]] <- as.numeric(draws[[col_name]])
      }
    }
  }

  if(length(edge_draws) == 0) stop("No posterior paths found. Check brms parameter names.")

  edge_keys <- names(edge_draws)
  edge_from <- sapply(strsplit(edge_keys, "->"), `[`, 1)
  edge_to   <- sapply(strsplit(edge_keys, "->"), `[`, 2)

  n_draws <- nrow(draws)
  N_vars  <- length(all_vars)

  out <- list()

  for(endpoint in endpoints){
    dir_store <- matrix(0, n_draws, length(preds), dimnames = list(NULL, preds))
    ind_store <- matrix(0, n_draws, length(preds), dimnames = list(NULL, preds))
    tot_store <- matrix(0, n_draws, length(preds), dimnames = list(NULL, preds))

    for(d in 1:n_draws){
      A <- matrix(0, N_vars, N_vars, dimnames = list(all_vars, all_vars))
      for(k in seq_along(edge_keys)){
        A[edge_from[k], edge_to[k]] <- edge_draws[[k]][d]
      }

      A2 <- A %*% A
      A3 <- A2 %*% A
      A4 <- A3 %*% A

      DIR <- A
      IND <- A2 + A3 + A4
      TOT <- DIR + IND

      for(i in seq_along(preds)){
        dir_store[d, i] <- DIR[preds[i], endpoint]
        ind_store[d, i] <- IND[preds[i], endpoint]
        tot_store[d, i] <- TOT[preds[i], endpoint]
      }
    }

    for(i in seq_along(preds)){
      s_dir <- summarize_posterior_vec(dir_store[, i])
      s_ind <- summarize_posterior_vec(ind_store[, i])
      s_tot <- summarize_posterior_vec(tot_store[, i])

      out[[length(out) + 1]] <- tibble::tibble(
        Predictor = nice_label(preds[i]),
        Outcome   = nice_label(endpoint),
        EffectType = "Direct",
        est = s_dir$est, ci_lo = s_dir$ci_lo, ci_hi = s_dir$ci_hi,
        pd = s_dir$pd, sig = s_dir$sig, stars = s_dir$stars
      )
      out[[length(out) + 1]] <- tibble::tibble(
        Predictor = nice_label(preds[i]),
        Outcome   = nice_label(endpoint),
        EffectType = "Indirect",
        est = s_ind$est, ci_lo = s_ind$ci_lo, ci_hi = s_ind$ci_hi,
        pd = s_ind$pd, sig = s_ind$sig, stars = s_ind$stars
      )
      out[[length(out) + 1]] <- tibble::tibble(
        Predictor = nice_label(preds[i]),
        Outcome   = nice_label(endpoint),
        EffectType = "Total",
        est = s_tot$est, ci_lo = s_tot$ci_lo, ci_hi = s_tot$ci_hi,
        pd = s_tot$pd, sig = s_tot$sig, stars = s_tot$stars
      )
    }
  }

  dplyr::bind_rows(out)
}

# -----------------------------
# 7) Plotting
# -----------------------------
plot_parallel_sem <- function(paths, group_title, available_mediators){
  med_labels <- c(
    nutrient_pc1_z = "Nutrient PC1",
    physchem_pc1_z = "PhysChem PC1",
    tsm_z = "TSM"
  )

  if(length(available_mediators) == 0){
    med_nodes <- tibble::tibble(
      name = character(0), x = numeric(0), y = numeric(0),
      label = character(0), fill_col = character(0)
    )
  } else {
    x_pos <- if(length(available_mediators) == 1) {
      0
    } else {
      seq(-2.3, 2.3, length.out = length(available_mediators))
    }

    med_nodes <- tibble::tibble(
      name = available_mediators,
      x = x_pos,
      y = 2.55,
      label = unname(med_labels[available_mediators]),
      fill_col = "#FFF8E1"
    )
  }

  nodes <- dplyr::bind_rows(
    tibble::tibble(
      name = "FI_z",
      x = 0,
      y = 4.0,
      label = "FI",
      fill_col = "#E1F5FE"
    ),
    med_nodes,
    tibble::tibble(
      name = c("shannon_z", "evenness_z", "beta_z"),
      x = c(-2.2, 0, 2.2),
      y = c(0.85, 0.85, 0.85),
      label = c("Shannon", "Evenness", "Beta"),
      fill_col = c("#E8F5E9", "#E8F5E9", "#F3E5F5")
    )
  )

  edges <- paths %>%
    dplyr::mutate(
      col  = ifelse(sig, ifelse(est >= 0, COL_POS, COL_NEG), "grey75"),
      lty  = ifelse(sig, "solid", "dashed"),
      w    = ifelse(sig, 0.7 + abs(est) * 2.0, 0.45),
      alpha = ifelse(sig, 1.0, 0.65),
      elab = sprintf("%.2f%s", est, stars),
      font_face = ifelse(sig, "bold", "plain"),
      text_col = ifelse(sig, "black", "grey45")
    ) %>%
    dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("from" = "name")) %>%
    dplyr::rename(x1 = x, y1 = y) %>%
    dplyr::left_join(nodes %>% dplyr::select(name, x, y), by = c("to" = "name")) %>%
    dplyr::rename(x2 = x, y2 = y) %>%
    dplyr::filter(is.finite(x1), is.finite(y1), is.finite(x2), is.finite(y2)) %>%
    dplyr::mutate(
      xm = 0.48 * x1 + 0.52 * x2,
      ym = 0.48 * y1 + 0.52 * y2
    )

  p <- ggplot2::ggplot() +
    ggplot2::xlim(-3.4, 3.4) + ggplot2::ylim(0.1, 4.45) +
    ggplot2::theme_void() +
    ggplot2::labs(title = group_title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
    )

  if(nrow(edges) > 0){
    p <- p +
      ggplot2::geom_segment(
        data = edges,
        ggplot2::aes(x = x1, y = y1, xend = x2, yend = y2,
                     colour = col, linewidth = w, linetype = lty, alpha = alpha),
        arrow = grid::arrow(length = grid::unit(2.3, "mm"), type = "closed"),
        lineend = "round",
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_linewidth_identity() +
      ggplot2::scale_linetype_identity() +
      ggplot2::scale_alpha_identity()
  }

  p <- p +
    ggplot2::geom_label(
      data = nodes,
      ggplot2::aes(x = x, y = y, label = label, fill = fill_col),
      fontface = "bold",
      size = 3.15,
      label.padding = grid::unit(0.20, "lines"),
      label.r = grid::unit(0.18, "lines"),
      color = "grey20",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_identity()

  if(nrow(edges) > 0){
    p <- p +
      ggrepel::geom_label_repel(
        data = edges,
        ggplot2::aes(x = xm, y = ym, label = elab, color = text_col, fontface = font_face),
        size = ifelse(edges$sig, 2.75, 2.25),
        fill = scales::alpha("white", 0.92),
        label.size = 0,
        box.padding = 0.12,
        point.padding = 0.05,
        min.segment.length = 0,
        seed = 123,
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity()
  }

  p
}

plot_fi_effects_group <- function(eff_df_group, group_title, y_limits = NULL){
  dat <- eff_df_group %>%
    dplyr::filter(Predictor == "FI") %>%
    dplyr::mutate(
      Outcome = factor(Outcome, levels = c("Shannon", "Evenness", "Beta")),
      EffectType = factor(EffectType, levels = c("Direct", "Indirect", "Total"))
    ) %>%
    dplyr::filter(!is.na(Outcome))

  if(nrow(dat) == 0){
    return(make_placeholder_plot(group_title, "No FI effect results"))
  }

  span <- max(dat$ci_hi, na.rm = TRUE) - min(dat$ci_lo, na.rm = TRUE)
  if(!is.finite(span) || span <= 0) span <- 1

  dat <- dat %>%
    dplyr::mutate(
      star_lab = ifelse(sig & nzchar(stars), stars, ""),
      lab_y = ifelse(est >= 0, ci_hi + 0.06 * span, ci_lo - 0.06 * span)
    )

  dodge <- ggplot2::position_dodge(width = 0.78)

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = Outcome, y = est, fill = EffectType)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey35", linewidth = 0.45) +
    ggplot2::geom_col(
      position = dodge,
      width = 0.68,
      color = "grey28",
      linewidth = 0.30
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ci_lo, ymax = ci_hi),
      position = dodge,
      width = 0.16,
      linewidth = 0.32,
      color = "grey20"
    ) +
    ggplot2::geom_text(
      data = dat %>% dplyr::filter(sig, nzchar(star_lab)),
      ggplot2::aes(y = lab_y, label = star_lab),
      position = dodge,
      size = 3.0,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = c(
      "Direct"   = "#9ECAE1",
      "Indirect" = "#A1D99B",
      "Total"    = "#FDD0A2"
    )) +
    ggplot2::labs(
      title = group_title,
      x = NULL,
      y = "Standardized effect of FI",
      fill = "Effect Type"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5),
      axis.text.x = ggplot2::element_text(face = "bold", size = 9.2, color = "black"),
      axis.text.y = ggplot2::element_text(size = 8.8, color = "black"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "grey35", linewidth = 0.55),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = 9.5, face = "bold"),
      legend.text = ggplot2::element_text(size = 8.8)
    )

  if(!is.null(y_limits)){
    p <- p + ggplot2::coord_cartesian(ylim = y_limits, clip = "off")
  }

  p
}

# -----------------------------
# 8) Read ENV + Build Indices
# -----------------------------
env_raw <- read_tsv_safe(ENV_FILE) %>% janitor::clean_names()

if(!("site" %in% names(env_raw)))  stop("env.txt must have a 'site' column.")
if(!("river" %in% names(env_raw))) env_raw$river <- "all_sites"

fi_col    <- find_col_relaxed(names(env_raw), c("FI"))
tsm_col   <- find_col_relaxed(names(env_raw), c("TSM"))
t_col     <- find_col_relaxed(names(env_raw), c("T"))
width_col <- find_col_relaxed(names(env_raw), c("river_width", "riverwidth", "width"))

if(is.na(fi_col))  stop("env.txt must contain an FI column.")
if(is.na(tsm_col)) stop("env.txt must contain a TSM column.")
if(is.na(width_col)) stop("env.txt must contain a river_width, riverwidth, or width column.")

nutrient_cols <- c(
  find_col_relaxed(names(env_raw), c("NH4")),
  find_col_relaxed(names(env_raw), c("NO2")),
  find_col_relaxed(names(env_raw), c("NO3")),
  find_col_relaxed(names(env_raw), c("PO4", "PO_4", "PO₄"))
)

physchem_cols <- c(
  find_col_relaxed(names(env_raw), c("T")),
  find_col_relaxed(names(env_raw), c("DO")),
  find_col_relaxed(names(env_raw), c("pH")),
  find_col_relaxed(names(env_raw), c("TDS")),
  find_col_relaxed(names(env_raw), c("Cond", "c", "Conductivity"))
)

nut_res <- compute_pc1_generic(
  df0 = env_raw,
  selected_cols = nutrient_cols,
  align_to = NULL,
  miss_keep = NULL
)

phys_res <- compute_pc1_generic(
  df0 = env_raw,
  selected_cols = physchem_cols,
  align_to = if(!is.na(t_col)) env_raw[[t_col]] else NULL,
  miss_keep = 0.8
)

env_raw$int_FI_river_width <- numify(env_raw[[fi_col]]) * numify(env_raw[[width_col]])

env_idx <- env_raw %>%
  dplyr::mutate(
    nutrient_pc1 = nut_res$scores,
    physchem_pc1 = phys_res$scores
  )

env_sem <- env_idx %>%
  dplyr::transmute(
    site = .data$site,
    river = as.character(.data$river),
    FI_z                 = z(.data[[fi_col]]),
    river_width_z        = z(.data[[width_col]]),
    int_FI_river_width_z = z(.data$int_FI_river_width),
    nutrient_pc1_z       = z(.data$nutrient_pc1),
    physchem_pc1_z       = z(.data$physchem_pc1),
    tsm_z                = z(.data[[tsm_col]])
  )

write.table(
  data.frame(
    variable = c("FI_z", "river_width_z", "int_FI_river_width_z", "nutrient_pc1_z", "physchem_pc1_z", "tsm_z"),
    sd = c(
      sd(env_sem$FI_z, na.rm = TRUE),
      sd(env_sem$river_width_z, na.rm = TRUE),
      sd(env_sem$int_FI_river_width_z, na.rm = TRUE),
      sd(env_sem$nutrient_pc1_z, na.rm = TRUE),
      sd(env_sem$physchem_pc1_z, na.rm = TRUE),
      sd(env_sem$tsm_z, na.rm = TRUE)
    )
  ),
  file = file.path(OUT_DIR, "ENV_variation_check.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

write.table(
  env_idx %>%
    dplyr::select(
      site, river,
      dplyr::all_of(fi_col),
      dplyr::all_of(tsm_col),
      dplyr::all_of(width_col),
      int_FI_river_width,
      nutrient_pc1, physchem_pc1
    ),
  file = file.path(OUT_DIR, "env_with_PC1_indices.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

if(!is.null(nut_res$loadings)){
  nut_load <- tibble::tibble(
    variable = names(nut_res$loadings),
    loading_PC1 = as.numeric(nut_res$loadings),
    contribution = (as.numeric(nut_res$loadings)^2) / sum(as.numeric(nut_res$loadings)^2)
  ) %>% dplyr::arrange(dplyr::desc(abs(loading_PC1)))

  write.table(
    nut_load,
    file = file.path(OUT_DIR, "Nutrient_PC1_Loadings.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

if(!is.null(phys_res$loadings)){
  phys_load <- tibble::tibble(
    variable = names(phys_res$loadings),
    loading_PC1 = as.numeric(phys_res$loadings),
    contribution = (as.numeric(phys_res$loadings)^2) / sum(as.numeric(phys_res$loadings)^2)
  ) %>% dplyr::arrange(dplyr::desc(abs(loading_PC1)))

  write.table(
    phys_load,
    file = file.path(OUT_DIR, "PhysChem_PC1_Loadings.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

write(
  paste0("SEM run started at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  file = file.path(OUT_DIR, "SEM_error_log.txt")
)

# -----------------------------
# 9) Main
# -----------------------------
all_paths_export   <- list()
all_effects_export <- list()
all_med_keep       <- list()

ENDPOINTS <- c("shannon_z", "evenness_z", "beta_z")

for(gname in GROUP_ORDER){
  cat(sprintf(
    "\n==============================\n>>> Strict SEM for %s <<<\n==============================\n",
    gname
  ))

  tryCatch({
    met <- compute_biodiv_metrics(
      spe_file = SPE_FILES[[gname]],
      env_info = env_sem %>% dplyr::select(site, river)
    ) %>%
      dplyr::left_join(env_sem, by = "site")

    df <- met %>%
      dplyr::mutate(
        log_total_reads_z = z(log(pmax(total_reads, 1))),
        shannon_z  = z(shannon),
        evenness_z = z(logit01(evenness)),
        beta_z     = z(logit01(beta_mean)),
        river      = factor(river)
      )

    dat_use <- df %>%
      dplyr::select(
        site, river,
        FI_z, river_width_z, int_FI_river_width_z, nutrient_pc1_z, physchem_pc1_z, tsm_z,
        shannon_z, evenness_z, beta_z, log_total_reads_z, beta_n_triplets
      ) %>%
      dplyr::filter(
        is.finite(FI_z),
        is.finite(river_width_z),
        is.finite(shannon_z),
        is.finite(evenness_z),
        is.finite(beta_z)
      )

    if(nrow(dat_use) < 15){
      message("Skipping ", gname, " due to insufficient data.")
      next
    }

    resp_check <- c("shannon_z", "evenness_z", "beta_z")
    bad_resp <- resp_check[!sapply(resp_check, function(v) has_variation(dat_use[[v]]))]
    if(length(bad_resp) > 0){
      stop(paste0("Zero variance in response(s): ", paste(bad_resp, collapse = ", ")))
    }

    if(!has_variation(dat_use$FI_z)){
      stop("FI_z has zero variance after site matching.")
    }

    med_keep <- pick_available_mediators(dat_use)
    all_med_keep[[gname]] <- med_keep

    rhs_base <- build_rhs_base(dat_use)

    rhs_re <- ""
    if("river" %in% names(dat_use)){
      rr <- droplevels(dat_use$river)
      if(nlevels(rr) > 1) rhs_re <- " + (1|river)"
    }

    med_term_str <- if(length(med_keep) > 0) paste(med_keep, collapse = " + ") else ""
    read_term <- if(CONTROL_READS_BIODIV && has_variation(dat_use$log_total_reads_z)) "log_total_reads_z" else ""

    rhs_shan <- paste_terms(rhs_base, med_term_str, read_term)
    rhs_even <- paste_terms(rhs_base, med_term_str, read_term)
    rhs_beta <- paste_terms(rhs_base, med_term_str, read_term)

    bforms <- list()

    if("nutrient_pc1_z" %in% med_keep){
      bforms[[length(bforms) + 1]] <- brms::bf(
        as.formula(paste0("nutrient_pc1_z ~ ", rhs_base, rhs_re))
      )
    }

    if("tsm_z" %in% med_keep){
      bforms[[length(bforms) + 1]] <- brms::bf(
        as.formula(paste0("tsm_z ~ ", rhs_base, rhs_re))
      )
    }

    if("physchem_pc1_z" %in% med_keep){
      bforms[[length(bforms) + 1]] <- brms::bf(
        as.formula(paste0("physchem_pc1_z ~ ", rhs_base, rhs_re))
      )
    }

    bforms[[length(bforms) + 1]] <- brms::bf(
      as.formula(paste0("shannon_z ~ ", rhs_shan, rhs_re))
    )
    bforms[[length(bforms) + 1]] <- brms::bf(
      as.formula(paste0("evenness_z ~ ", rhs_even, rhs_re))
    )
    bforms[[length(bforms) + 1]] <- brms::bf(
      as.formula(paste0("beta_z ~ ", rhs_beta, rhs_re))
    )

    bform <- bforms[[1]]
    if(length(bforms) > 1){
      for(ii in 2:length(bforms)){
        bform <- bform + bforms[[ii]]
      }
    }
    bform <- bform + brms::set_rescor(FALSE)

    write(
      paste0(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        " | ", gname,
        " | rhs_base = ", rhs_base,
        " | mediators_kept = ", if(length(med_keep) == 0) "none" else paste(med_keep, collapse = ", "),
        " | model = parallel responses"
      ),
      file = file.path(OUT_DIR, "SEM_error_log.txt"),
      append = TRUE
    )

    fit <- brms::brm(
      formula = bform,
      data = dat_use,
      chains = MCMC_CHAINS,
      cores = max(1, parallel::detectCores() - 1),
      iter = MCMC_ITER,
      warmup = MCMC_WARMUP,
      control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH),
      backend = "rstan",
      silent = 2,
      refresh = 0,
      seed = 123
    )

    fit_name <- paste0("BiodivSEM_Parallel_FIonly_", gname)
    save_brms_diagnostics(fit, fit_name, OUT_DIR)

    focal_preds <- "FI_z"
    sem_vars_for_cascade <- unique(c("FI_z", "river_width_z", "int_FI_river_width_z", med_keep, ENDPOINTS))

    paths_df <- extract_paths_table(fit, all_vars = sem_vars_for_cascade)
    eff_df <- extract_cascading_effects_multi(
      fit,
      preds = focal_preds,
      endpoints = ENDPOINTS,
      all_vars = sem_vars_for_cascade
    )

    paths_df <- paths_df %>% dplyr::mutate(Group = gname, .before = 1)
    eff_df   <- eff_df %>% dplyr::mutate(Group = gname, .before = 1)

    write.table(
      paths_df,
      file = file.path(OUT_DIR, paste0("SEM_", gname, "_Paths_Parallel.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    write.table(
      eff_df,
      file = file.path(OUT_DIR, paste0("SEM_", gname, "_FI_ParallelCascades.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    # Export the number of beta-diversity triplets per site for quality checking.
    write.table(
      dat_use %>% dplyr::select(site, river, beta_n_triplets),
      file = file.path(OUT_DIR, paste0("SEM_", gname, "_BetaTripletCounts.tsv")),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    all_paths_export[[gname]]   <- paths_df
    all_effects_export[[gname]] <- eff_df

    message("DONE: ", gname)

    rm(fit, paths_df, eff_df, dat_use, df, met)
    gc()

  }, error = function(e){
    msg <- paste0(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      " | FAILED: ", gname,
      " | ", conditionMessage(e)
    )

    message(msg)

    write(
      msg,
      file = file.path(OUT_DIR, "SEM_error_log.txt"),
      append = TRUE
    )
  })
}

# -----------------------------
# 9.5) Global Y-Limits for outcome plots
# -----------------------------
global_y_min <- 0
global_y_max <- 0

if(length(all_effects_export) > 0){
  combined_eff_df <- dplyr::bind_rows(all_effects_export) %>%
    dplyr::filter(Predictor == "FI")
  if(nrow(combined_eff_df) > 0){
    global_y_min <- min(c(0, combined_eff_df$ci_lo), na.rm = TRUE)
    global_y_max <- max(c(0, combined_eff_df$ci_hi), na.rm = TRUE)
    pad <- (global_y_max - global_y_min) * 0.08
    if(!is.finite(pad) || pad == 0) pad <- 0.1
    global_y_min <- global_y_min - pad
    global_y_max <- global_y_max + pad
  }
}

y_limits_global <- c(global_y_min, global_y_max)

# -----------------------------
# 9.6) Generate plots
# -----------------------------
net_plots <- vector("list", length(GROUP_ORDER))
names(net_plots) <- GROUP_ORDER

for(gname in GROUP_ORDER){
  if(gname %in% names(all_paths_export)){
    net_plots[[gname]] <- plot_parallel_sem(
      all_paths_export[[gname]],
      GROUP_LABELS[gname],
      available_mediators = all_med_keep[[gname]]
    )
  } else {
    net_plots[[gname]] <- make_placeholder_plot(GROUP_LABELS[gname], "Model failed/No data")
  }
}

all_effects_df <- NULL
if(length(all_effects_export) > 0){
  all_effects_df <- dplyr::bind_rows(all_effects_export)
}

fi_plots <- vector("list", length(GROUP_ORDER))
names(fi_plots) <- GROUP_ORDER

for(gname in GROUP_ORDER){
  if(!is.null(all_effects_df) && gname %in% names(all_effects_export)){
    fi_plots[[gname]] <- plot_fi_effects_group(
      eff_df_group = all_effects_export[[gname]],
      group_title = GROUP_LABELS[gname],
      y_limits = y_limits_global
    )
  } else {
    fi_plots[[gname]] <- make_placeholder_plot(GROUP_LABELS[gname], "No FI effect results")
  }
}

# -----------------------------
# 10) Export merged tables
# -----------------------------
if(length(all_paths_export) > 0){
  all_paths_df <- dplyr::bind_rows(all_paths_export)
  write.table(
    all_paths_df,
    file = file.path(OUT_DIR, "SEM_AllGroups_Paths_Parallel.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

if(length(all_effects_export) > 0){
  write.table(
    all_effects_df,
    file = file.path(OUT_DIR, "SEM_AllGroups_FI_ParallelCascades.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

# -----------------------------
# 11) Combined figure
# -----------------------------
top_row <- patchwork::wrap_plots(net_plots[GROUP_ORDER], ncol = 4)
bottom_row <- patchwork::wrap_plots(fi_plots[GROUP_ORDER], ncol = 4)

p_all <- make_row_title("Bayesian SEM (parallel biodiversity responses)") /
  top_row /
  make_row_title("Direct, indirect, and total effects of FI on Shannon, Evenness, and Beta") /
  bottom_row +
  patchwork::plot_layout(heights = c(0.05, 0.95, 0.05, 1.00), guides = "collect") +
  patchwork::plot_annotation(
    title = "Integrated Bayesian Piecewise SEM (Width & River Controlled)",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 18, hjust = 0.5)
    )
  ) &
  ggplot2::theme(legend.position = "bottom")

save_pdf_tmp(
  file.path(OUT_DIR, "Fig_StrictSEM_4Groups_ParallelResponses_WidthControlled.pdf"),
  width = 24, height = 12.5,
  plot_fun = function() print(p_all)
)

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "Fig_StrictSEM_4Groups_ParallelResponses_WidthControlled.png"),
  plot = p_all, width = 24, height = 12.5, dpi = 500, bg = "white"
)

cat("
All analyses are complete. Outputs saved to:
", OUT_DIR, "
")
cat("\nNutrient PC1 used columns:", paste(nut_res$used_cols, collapse = ", "), "\n")
cat("PhysChem PC1 used columns:", paste(phys_res$used_cols, collapse = ", "), "\n")
cat("USE_RIVER_RANDOM_EFFECT =", USE_RIVER_RANDOM_EFFECT, "\n")
cat("CONTROL_READS_BIODIV   =", CONTROL_READS_BIODIV, "\n")
cat("Beta diversity now uses ALL 3-site combinations within each river.\n")