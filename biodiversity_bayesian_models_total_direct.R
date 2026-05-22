# ============================================================
# Publication-ready Bayesian biodiversity models using JAGS
# Responses:
#   1. Shannon diversity
#   2. Pielou evenness
#   3. Site-level local beta diversity
#
# Model sets:
#   TOTAL  - excludes potential environmental mediators to estimate an overall association.
#   DIRECT - includes environmental variables to estimate conditional direct effects.
#
# Main design features:
#   - Controls for river identity using random intercepts.
#   - Controls for river width and mainstem status.
#   - Computes beta_local from all four-site combinations within the same river.
#   - Runs four organismal groups: bacteria, protists, metazoans, and fishes.
#   - Exports posterior summaries, diagnostics, trace plots, HOP plots, and PDP plots.
#
# Required input files in the working directory:
#   env.txt
#     Required columns: site, river, FI, river_width
#   spe_bacteria.txt
#   spe_protist.txt
#   spe_Metazoa.txt
#   spe_fish.txt
#
# Notes for public repositories:
#   - Exact sampling coordinates are not required by this script.
#   - Use anonymized site IDs if sampling locations are sensitive.
#   - Keep raw data and large sequencing files out of this repository unless intended.
# ============================================================

options(stringsAsFactors = FALSE)
set.seed(123)

# ------------------------------------------------------------
# 0) Speed and quality settings
# ------------------------------------------------------------
# Use "fast" for quick testing and "final" for publication-level MCMC settings.
run_mode <- "final"  # "fast" or "final"

mcmc_alpha_fast  <- list(n_chains=3, n_adapt=800,  n_update=1500, n_iter=5000,  thin=8)
mcmc_alpha_final <- list(n_chains=4, n_adapt=8000, n_update=15000, n_iter=50000, thin=8)

# PDP and HOP plotting settings
pdp_grid_fast    <- 60
pdp_grid_final   <- 100
n_lines_fast     <- 250
n_lines_final    <- 800

if(run_mode=="fast"){
  mcmc_alpha <- mcmc_alpha_fast
  pdp_grid <- pdp_grid_fast
  n_lines_show <- n_lines_fast
} else {
  mcmc_alpha <- mcmc_alpha_final
  pdp_grid <- pdp_grid_final
  n_lines_show <- n_lines_final
}

# ------------------------------------------------------------
# 1) Install and load required packages
# ------------------------------------------------------------
pkgs <- c("vegan","rjags","coda","gridExtra","grid","ggplot2","viridisLite")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install)>0) install.packages(to_install, repos="https://cloud.r-project.org")

suppressPackageStartupMessages({
  library(vegan)
  library(rjags)
  library(coda)
  library(gridExtra)
  library(grid)
  library(ggplot2)
  library(viridisLite)
})

# ------------------------------------------------------------
# 2) Global model settings
# ------------------------------------------------------------
# TOTAL estimates overall associations; DIRECT estimates conditional effects after controlling mediators.
MODEL_TYPES <- c("TOTAL","DIRECT")

# Main predictor of interest. Rename here if the forest gradient variable has a different column name.
focus_main <- "FI"
mainstem_var <- "mainstem"
width_var    <- "river_width"
int_terms <- c("int_FI_mainstem", "int_FI_river_width")
focus_report_terms <- c(focus_main, int_terms)
force_pred_terms  <- c(focus_main, mainstem_var, width_var, int_terms)

# Predictor screening and collinearity settings
cor_cut <- 0.85
max_alpha_predictors <- 12
min_N_complete_sites <- 20

# Missing-value handling
IMPUTE_ENV_DIRECT <- TRUE
IMPUTE_ENV_TOTAL  <- FALSE

# Technical covariate: sequencing depth
CONTROL_READS_ALPHA <- TRUE

# Beta-diversity settings based on four-site combinations
BETA_COMB_SIZE <- 4
BETA_USE_RELATIVE_ABUNDANCE <- TRUE
BETA_MIN_COMBS_PER_SITE <- 1

# Output size for overview PDFs
OV_PAGE_W <- 18.0
OV_PAGE_H <- 11.8

# Color-bar size
CBAR_W <- 0.22
CBAR_H <- 0.95

STRIP_OVERVIEW_TITLES <- TRUE

# MCMC diagnostic output
SAVE_TRACE_PDF <- TRUE
SAVE_DIAG_CSV  <- TRUE
Rhat_threshold <- 1.01
min_ESS_key    <- 800

# ------------------------------------------------------------
# 3) Check external JAGS installation
# ------------------------------------------------------------
check_jags_or_stop <- function(){
  ver <- tryCatch(rjags::jags.version(), error=function(e) e)
  if(inherits(ver, "error")){
    stop("JAGS is not available: ", conditionMessage(ver),
         "\nPlease install the JAGS software first, then install the rjags R package, and verify the installation with rjags::jags.version().")
  }
  message("✅ JAGS OK. version = ", as.character(ver))
}
check_jags_or_stop()

# ------------------------------------------------------------
# 4) Utility functions
# ------------------------------------------------------------
invlogit <- function(x) 1/(1+exp(-x))

safe_blank_plot <- function(title="FAILED"){
  ggplot() +
    theme_void() +
    annotate("text", x=0, y=0, label=title, size=4, fontface="bold") +
    xlim(-1,1) + ylim(-1,1)
}

strip_title <- function(p){
  if(inherits(p, "ggplot")) p + theme(plot.title=element_blank()) else p
}

# Standardize ggplot grob dimensions in overview panels.
# This keeps panel proportions consistent across response rows.
normalize_overview_plot_grobs <- function(grobs){
  grobs2 <- lapply(grobs, function(g){
    if(inherits(g, "ggplot")) ggplotGrob(g) else g
  })

  is_gtable <- vapply(grobs2, function(g) inherits(g, "gtable"), logical(1))
  if(any(is_gtable)){
    max_widths  <- Reduce(grid::unit.pmax, lapply(grobs2[is_gtable], function(g) g$widths))
    max_heights <- Reduce(grid::unit.pmax, lapply(grobs2[is_gtable], function(g) g$heights))

    grobs2[is_gtable] <- lapply(grobs2[is_gtable], function(g){
      g$widths  <- max_widths
      g$heights <- max_heights
      g
    })
  }
  grobs2
}

scale_fill_viridis_safe <- function(){
  ggplot2::scale_fill_gradientn(colours = viridisLite::viridis(256))
}

PLOT_FAMILY <- "sans"

open_pdf_tmp <- function(final_path, width, height){
  tmp_path <- paste0(final_path, ".tmp.pdf")
  if(file.exists(tmp_path)) suppressWarnings(file.remove(tmp_path))

  if(capabilities("cairo")){
    grDevices::cairo_pdf(tmp_path, width=width, height=height, onefile=TRUE, family=PLOT_FAMILY)
  } else {
    grDevices::pdf(tmp_path, width=width, height=height, onefile=TRUE, family=PLOT_FAMILY)
  }
  tmp_path
}

close_pdf_tmp <- function(tmp_path, final_path){
  grDevices::dev.off()
  if(file.exists(final_path)) suppressWarnings(file.remove(final_path))
  ok <- file.rename(tmp_path, final_path)
  if(!isTRUE(ok)){
    warning("Failed to overwrite the existing PDF file. The new temporary file was kept at: ", tmp_path,
            "\nPlease close the existing PDF and rename the temporary file to: ", final_path)
  }
}

# Variables listed here are treated as potential environmental mediators and are excluded from TOTAL models.
MEDIATOR_EXACT <- c(
  "t","temp","temperature","do","oxygen","dissolved_oxygen",
  "ph","tds","cond","conductivity",
  "tsm","tss","turbidity","ntu","secchi",
  "nh4","no2","no3","po4","si","din","tin","tn","tp",
  "chla","chlorophyll","chlorophyll_a",
  "sal","salinity","sss","sst"
)

is_mediator_var <- function(v){
  v2 <- tolower(v)
  v2 <- sub("^d_", "", v2)
  v2 <- gsub("[^a-z0-9_]+", "_", v2)
  v2 %in% MEDIATOR_EXACT
}

filter_predictors_by_model <- function(preds_all, force_in, model_type){
  preds_all <- unique(preds_all)
  force_in  <- unique(force_in)
  if(model_type == "TOTAL"){
    keep <- preds_all[!sapply(preds_all, is_mediator_var)]
    unique(c(force_in, keep))
  } else {
    unique(c(force_in, preds_all))
  }
}

join_by_site_keep_order <- function(df, env2){
  out <- df
  idx <- match(out$site, env2$site)
  for(nm in names(env2)){
    if(nm == "site") next
    out[[nm]] <- env2[[nm]][idx]
  }
  out
}

median_impute_numeric <- function(df, cols){
  impute_log <- data.frame(var=character(0), n_na_before=integer(0), n_imputed=integer(0), stringsAsFactors=FALSE)
  for(nm in cols){
    x <- df[[nm]]
    if(!is.numeric(x)) next
    na_idx <- which(!is.finite(x))
    if(length(na_idx) > 0){
      med <- suppressWarnings(median(x[is.finite(x)], na.rm=TRUE))
      if(!is.finite(med)) next
      df[[nm]][na_idx] <- med
      impute_log <- rbind(impute_log, data.frame(var=nm, n_na_before=length(na_idx), n_imputed=length(na_idx), stringsAsFactors=FALSE))
    }
  }
  list(df=df, log=impute_log)
}

collinearity_filter <- function(X0, cor_cut=0.85, focus_set=character()){
  C <- suppressWarnings(cor(X0))
  diag(C) <- 0
  Z_names <- colnames(X0)
  focus_set <- intersect(focus_set, Z_names)

  drop_vars <- c()
  if(ncol(C) > 2){
    repeat{
      absC <- abs(C)
      if(length(focus_set) >= 2){
        absC[focus_set, focus_set] <- 0
      }
      mx <- max(absC, na.rm=TRUE)
      if(!is.finite(mx) || mx <= cor_cut) break

      ij <- which(absC==mx, arr.ind=TRUE)[1,]
      v1 <- colnames(absC)[ij[1]]
      v2 <- colnames(absC)[ij[2]]

      v1_is_focus <- v1 %in% focus_set
      v2_is_focus <- v2 %in% focus_set

      if(v1_is_focus && !v2_is_focus){
        drop_one <- v2
      } else if(v2_is_focus && !v1_is_focus){
        drop_one <- v1
      } else {
        m1 <- mean(abs(C[v1,]), na.rm=TRUE)
        m2 <- mean(abs(C[v2,]), na.rm=TRUE)
        drop_one <- ifelse(m1 >= m2, v1, v2)
      }

      drop_vars <- c(drop_vars, drop_one)
      keep <- setdiff(colnames(C), drop_one)
      C <- C[keep, keep, drop=FALSE]
      focus_set <- intersect(focus_set, keep)
      if(ncol(C) < 2) break
    }
  }
  list(keep=setdiff(Z_names, unique(drop_vars)), drop=unique(drop_vars))
}

# ------------------------------------------------------------
# 5) Compute site-level local beta diversity from all four-site combinations within each river
# ------------------------------------------------------------
# This function converts the community matrix into a site-level beta-diversity response.
# For each river, it builds all possible four-site combinations and summarizes Bray-Curtis distances.
compute_site_beta_local <- function(X_asv, env2,
                                    comb_size = BETA_COMB_SIZE,
                                    use_relative_abundance = BETA_USE_RELATIVE_ABUNDANCE){
  stopifnot(all(rownames(X_asv) %in% env2$site))

  mat <- as.matrix(X_asv)
  mat <- mat[env2$site, , drop=FALSE]
  rownames(mat) <- env2$site

  rs <- rowSums(mat, na.rm=TRUE)
  keep <- rs > 0 & !is.na(env2$river)
  mat <- mat[keep, , drop=FALSE]
  meta <- env2[keep, c("site","river","FI","river_width","mainstem"), drop=FALSE]

  if(use_relative_abundance){
    mat <- sweep(mat, 1, rowSums(mat), "/")
    mat[!is.finite(mat)] <- 0
  }

  rec <- list()
  idx <- 1
  rivers <- unique(meta$river)

  for(rv in rivers){
    meta_r <- meta[meta$river == rv, , drop=FALSE]
    n_r <- nrow(meta_r)
    if(n_r < comb_size) next

    comb_list <- combn(meta_r$site, comb_size, simplify=FALSE)

    for(cm in comb_list){
      sub_meta <- meta_r[match(cm, meta_r$site), , drop=FALSE]
      sub_mat  <- mat[cm, , drop=FALSE]

      bc <- as.matrix(vegan::vegdist(sub_mat, method="bray"))
      diag(bc) <- NA_real_

      # For each site in a four-site combination, compute the mean Bray-Curtis distance to the other three sites.
      local_beta <- rowMeans(bc, na.rm=TRUE)

      fi_mean <- mean(sub_meta$FI, na.rm=TRUE)
      width_mean <- mean(sub_meta$river_width, na.rm=TRUE)
      mainstem_mean <- mean(sub_meta$mainstem, na.rm=TRUE)

      for(ss in names(local_beta)){
        rec[[idx]] <- data.frame(
          site = ss,
          river = rv,
          comb_beta_local = as.numeric(local_beta[ss]),
          comb_FI_mean = fi_mean,
          comb_river_width_mean = width_mean,
          comb_mainstem_mean = mainstem_mean,
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
  }

  if(length(rec) == 0){
    return(data.frame(site=env2$site, beta_local=NA_real_, n_beta_combs=0, stringsAsFactors=FALSE))
  }

  W <- do.call(rbind, rec)

  agg_mean <- aggregate(comb_beta_local ~ site, data=W, FUN=mean)
  agg_n    <- aggregate(comb_beta_local ~ site, data=W, FUN=length)

  names(agg_mean)[2] <- "beta_local"
  names(agg_n)[2]    <- "n_beta_combs"

  out <- merge(data.frame(site=env2$site, stringsAsFactors=FALSE), agg_mean, by="site", all.x=TRUE)
  out <- merge(out, agg_n, by="site", all.x=TRUE)
  out$n_beta_combs[is.na(out$n_beta_combs)] <- 0
  out
}

# ------------------------------------------------------------
# 6) Plotting helper functions
# ------------------------------------------------------------
posterior_density_gg <- function(draws){
  if(is.null(draws)) return(safe_blank_plot("NO DRAWS"))
  draws <- as.numeric(draws)
  draws <- draws[is.finite(draws)]
  if(length(draws) < 50) return(safe_blank_plot("NO DRAWS"))

  d <- density(draws)
  df <- data.frame(x=d$x, y=d$y)

  q <- quantile(draws, c(0.05,0.95), na.rm=TRUE)
  mu <- mean(draws, na.rm=TRUE)
  col_line <- ifelse(mu >= 0, "dodgerblue3", "firebrick3")

  ggplot(df, aes(x,y)) +
    geom_line(linewidth=1.1, color=col_line) +
    geom_vline(xintercept=0, linetype=2, color="grey55") +
    geom_segment(aes(x=q[1], xend=q[2], y=0, yend=0), linewidth=2.5, color="grey35") +
    geom_point(aes(x=mu, y=0), size=2.5) +
    theme_bw(base_size=12, base_family=PLOT_FAMILY) +
    theme(
      axis.title=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank(),
      plot.margin=margin(10, 10, 10, 10)
    )
}

make_overview_grid <- function(grobs, row_labels, col_labels, title, file,
                               n_rows, n_cols, page_w=OV_PAGE_W, page_h=OV_PAGE_H,
                               left_w=0.9, top_h=0.55, title_h=0.70){

  stopifnot(length(grobs) == n_rows*n_cols)
  grobs <- normalize_overview_plot_grobs(grobs)
  blank <- grid::nullGrob()

  cell_w <- (page_w - left_w) / n_cols
  cell_h <- (page_h - title_h - top_h) / n_rows

  col_grobs <- lapply(col_labels, function(x)
    grid::textGrob(x, gp=grid::gpar(fontsize=11, fontface="bold", fontfamily=PLOT_FAMILY)))
  row_grobs <- lapply(row_labels, function(x)
    grid::textGrob(x, rot=90, gp=grid::gpar(fontsize=11, fontface="bold", fontfamily=PLOT_FAMILY)))

  items <- list(blank)
  items <- c(items, col_grobs)

  k <- 1
  for(r in 1:n_rows){
    items <- c(items, list(row_grobs[[r]]))
    for(c in 1:n_cols){
      items <- c(items, list(grobs[[k]]))
      k <- k + 1
    }
  }

  g_main <- gridExtra::arrangeGrob(
    grobs=items, ncol=(n_cols+1),
    widths = grid::unit.c(grid::unit(left_w,"in"), rep(grid::unit(cell_w,"in"), n_cols)),
    heights= grid::unit.c(grid::unit(top_h,"in"), rep(grid::unit(cell_h,"in"), n_rows))
  )

  g_all <- gridExtra::arrangeGrob(
    grid::textGrob(title, gp=grid::gpar(fontsize=13, fontface="bold", fontfamily=PLOT_FAMILY)),
    g_main,
    ncol=1,
    heights = grid::unit.c(grid::unit(title_h,"in"), grid::unit(1,"null"))
  )

  tmp <- open_pdf_tmp(file, width=page_w, height=page_h)
  grid::grid.newpage()
  grid::grid.draw(g_all)
  close_pdf_tmp(tmp, file)
}

format_prob3 <- function(x, digits=3){
  if(length(x)==0 || !is.finite(x)) return("NA")
  sprintf(paste0("%.", digits, "f"), x)
}

make_pd_hop_text <- function(draws, var_label, digits=3){
  if(is.null(draws)) return(paste0("pd(", var_label, ">0)=NA | HOP=NA"))
  draws <- as.numeric(draws)
  draws <- draws[is.finite(draws)]
  if(length(draws) < 10) return(paste0("pd(", var_label, ">0)=NA | HOP=NA"))

  pd_pos <- mean(draws > 0)
  hop    <- max(pd_pos, 1 - pd_pos)

  paste0(
    "pd(", var_label, ">0)=", format_prob3(pd_pos, digits),
    " | HOP=", format_prob3(hop, digits)
  )
}

hop_multi <- function(pred_raw, metric=c("shannon","evenness","beta_local"),
                      dir_mat, Z_names, Z_center, Z_scale, dat_z, fitted_rows,
                      n_lines=300, stat_text=NULL){
  metric <- match.arg(metric)
  pred_z <- paste0(pred_raw, "_z")
  if(!pred_z %in% Z_names) return(safe_blank_plot(paste0(pred_raw, " missing")))
  j <- which(Z_names == pred_z)

  x_raw <- dat_z[[pred_raw]][fitted_rows]
  x_raw <- x_raw[is.finite(x_raw)]
  if(length(x_raw) < 10) return(safe_blank_plot(paste0(pred_raw, " no data")))

  xgrid_raw <- seq(min(x_raw), max(x_raw), length.out=140)
  xgrid_z <- (xgrid_raw - Z_center[pred_raw]) / Z_scale[pred_raw]

  if(metric=="shannon"){
    a <- dir_mat[, "alpha[1]"]
    b <- dir_mat[, paste0("beta[", j, ",1]")]
    pred_mat <- sweep(outer(b, xgrid_z, "*"), 1, a, "+")
    ylab <- "Predicted Shannon (H')"
    title <- paste0("HOP: ", pred_raw, " -> Shannon")
  } else if(metric=="evenness"){
    a <- dir_mat[, "alpha[2]"]
    b <- dir_mat[, paste0("beta[", j, ",2]")]
    pred_logit <- sweep(outer(b, xgrid_z, "*"), 1, a, "+")
    pred_mat <- invlogit(pred_logit)
    ylab <- "Predicted Evenness"
    title <- paste0("HOP: ", pred_raw, " -> Evenness")
  } else {
    a <- dir_mat[, "alpha[3]"]
    b <- dir_mat[, paste0("beta[", j, ",3]")]
    pred_logit <- sweep(outer(b, xgrid_z, "*"), 1, a, "+")
    pred_mat <- invlogit(pred_logit)
    ylab <- "Predicted Beta_local"
    title <- paste0("HOP: ", pred_raw, " -> Beta_local")
  }

  keep_n <- min(n_lines, nrow(pred_mat))
  idx_show <- sample(seq_len(nrow(pred_mat)), keep_n)

  sum_df <- data.frame(
    x_raw=xgrid_raw,
    mid=apply(pred_mat,2,median),
    lo =apply(pred_mat,2,quantile,0.05),
    hi =apply(pred_mat,2,quantile,0.95)
  )
  line_df <- do.call(rbind, lapply(idx_show, function(ii){
    data.frame(x_raw=xgrid_raw, y=pred_mat[ii,], draw=ii)
  }))

  p <- ggplot() +
    geom_line(data=line_df, aes(x=x_raw, y=y, group=draw), color="dodgerblue", alpha=0.05) +
    geom_ribbon(data=sum_df, aes(x=x_raw, ymin=lo, ymax=hi), fill="dodgerblue", alpha=0.30) +
    geom_line(data=sum_df, aes(x=x_raw, y=mid), color="black", linewidth=1.15) +
    theme_bw(base_size=11, base_family=PLOT_FAMILY) +
    theme(plot.title=element_text(size=10, face="bold", hjust=0.5, lineheight=0.95),
          plot.margin=margin(5,5,5,5)) +
    labs(x=paste0(pred_raw, " (raw)"), y=ylab, title=title)

  if(!is.null(stat_text) && nzchar(stat_text)){
    p <- p + labs(subtitle = stat_text) +
      theme(plot.subtitle = element_text(size=9, face="bold", color="firebrick3", hjust=0.5))
  }

  return(p)
}

pdp_multi_simple <- function(pred1_raw, pred2_raw, metric=c("shannon","evenness","beta_local"),
                             dir_mat, Z_names, Z_center, Z_scale, dat_z, fitted_rows, n_grid=80){

  metric <- match.arg(metric)
  p1z <- paste0(pred1_raw, "_z")
  p2z <- paste0(pred2_raw, "_z")
  if(!(p1z %in% Z_names) || !(p2z %in% Z_names)) return(safe_blank_plot("Missing pred"))
  j1 <- which(Z_names==p1z); j2 <- which(Z_names==p2z)

  x1 <- dat_z[[pred1_raw]][fitted_rows]; x1 <- x1[is.finite(x1)]
  x2 <- dat_z[[pred2_raw]][fitted_rows]; x2 <- x2[is.finite(x2)]
  if(length(x1)<10 || length(x2)<10) return(safe_blank_plot("No data"))

  g1 <- seq(min(x1), max(x1), length.out=n_grid)
  g2 <- seq(min(x2), max(x2), length.out=n_grid)

  df <- expand.grid(x1=g1, x2=g2)
  df$z1 <- (df$x1 - Z_center[pred1_raw]) / Z_scale[pred1_raw]
  df$z2 <- (df$x2 - Z_center[pred2_raw]) / Z_scale[pred2_raw]

  if(metric=="shannon"){
    a0 <- mean(dir_mat[, "alpha[1]"])
    b1 <- mean(dir_mat[, paste0("beta[", j1, ",1]")])
    b2 <- mean(dir_mat[, paste0("beta[", j2, ",1]")])
    df$y <- a0 + b1*df$z1 + b2*df$z2
    fill_lab <- "Shannon"
    ttl <- paste0("PDP: ", pred1_raw, " × ", pred2_raw, " -> Shannon")
  } else if(metric=="evenness") {
    a0 <- mean(dir_mat[, "alpha[2]"])
    b1 <- mean(dir_mat[, paste0("beta[", j1, ",2]")])
    b2 <- mean(dir_mat[, paste0("beta[", j2, ",2]")])
    df$logit <- a0 + b1*df$z1 + b2*df$z2
    df$y <- invlogit(df$logit)
    fill_lab <- "Evenness"
    ttl <- paste0("PDP: ", pred1_raw, " × ", pred2_raw, " -> Evenness")
  } else {
    a0 <- mean(dir_mat[, "alpha[3]"])
    b1 <- mean(dir_mat[, paste0("beta[", j1, ",3]")])
    b2 <- mean(dir_mat[, paste0("beta[", j2, ",3]")])
    df$logit <- a0 + b1*df$z1 + b2*df$z2
    df$y <- invlogit(df$logit)
    fill_lab <- "Beta_local"
    ttl <- paste0("PDP: ", pred1_raw, " × ", pred2_raw, " -> Beta_local")
  }

  ggplot(df, aes(x1, x2, fill=y)) +
    geom_raster() +
    scale_fill_viridis_safe() +
    guides(fill=guide_colorbar(barwidth=unit(CBAR_W,"in"), barheight=unit(CBAR_H,"in"))) +
    theme_bw(base_size=10, base_family=PLOT_FAMILY) +
    theme(
      plot.title=element_text(size=10, face="bold", hjust=0.5),
      axis.title=element_text(size=9),
      legend.title=element_text(size=9),
      legend.text=element_text(size=8),
      plot.margin=margin(4,4,4,4),
      aspect.ratio=1
    ) +
    labs(x=pred1_raw, y=pred2_raw, fill=fill_lab, title=ttl)
}

get_beta_summary_multi <- function(samp, Z_names){
  M <- as.matrix(samp)
  out <- list()
  resp_names <- c("shannon","evenness","beta_local")
  for(k in 1:3){
    for(j in seq_along(Z_names)){
      nm <- paste0("beta[", j, ",", k, "]")
      if(!nm %in% colnames(M)) next
      v <- M[, nm]
      hpd <- coda::HPDinterval(as.mcmc(v), prob=0.90)
      out[[length(out)+1]] <- data.frame(
        response=resp_names[k],
        predictor=Z_names[j],
        post_mean=mean(v),
        CI90_low=as.numeric(hpd[1,1]),
        CI90_high=as.numeric(hpd[1,2]),
        P_gt_0=mean(v>0),
        stringsAsFactors=FALSE
      )
    }
  }
  if(length(out)==0) return(data.frame())
  do.call(rbind, out)
}

# ------------------------------------------------------------
# 7) Load environmental data
# ------------------------------------------------------------
if(!file.exists("env.txt")) stop("env.txt not found!")
env <- read.delim("env.txt", sep="\t", check.names=FALSE, stringsAsFactors=FALSE)

names(env) <- trimws(names(env))
if(!"site" %in% names(env)) stop("env.txt must contain 'site' column.")
if(!"river" %in% names(env)) stop("env.txt must contain 'river' column.")
if(!"FI" %in% names(env) && "fi" %in% names(env)) env$FI <- env$fi
if(!"FI" %in% names(env)) stop("env.txt must contain 'FI' column.")
if(!"river_width" %in% names(env) && "riverwidth" %in% names(env)) env$river_width <- env$riverwidth
if(!"river_width" %in% names(env)) stop("env.txt must contain 'river_width' (or 'riverwidth') column.")

group_candidates <- c("season","campaign","time","year","basin","group")
group_col <- group_candidates[group_candidates %in% names(env)]
group_col <- if(length(group_col)>0) group_col[1] else NULL

env_num <- env
to_num_cols <- setdiff(names(env_num), c("site","river", group_col))
for(nm in to_num_cols){
  env_num[[nm]] <- suppressWarnings(as.numeric(trimws(env_num[[nm]])))
}

env_num$mainstem <- ifelse(env_num$river == "tmj" | grepl("^tmj", env_num$site), 1L, 0L)
env_num$int_FI_mainstem    <- env_num$FI * env_num$mainstem
env_num$int_FI_river_width <- env_num$FI * env_num$river_width

# ------------------------------------------------------------
# 8) MCMC diagnostic helper functions
# ------------------------------------------------------------
mcmc_diag_key <- function(samp, key_params){
  out <- data.frame(param=key_params, mean=NA_real_, sd=NA_real_,
                    q05=NA_real_, q50=NA_real_, q95=NA_real_,
                    ESS=NA_real_, Rhat=NA_real_, stringsAsFactors=FALSE)

  M <- as.matrix(samp)
  cn <- colnames(M)
  keep <- key_params[key_params %in% cn]
  if(length(keep) == 0) return(out)

  ess_all <- tryCatch(coda::effectiveSize(samp), error=function(e) NULL)
  rhat_all <- tryCatch({
    gd <- coda::gelman.diag(samp[, keep, drop=FALSE], multivariate=FALSE)
    as.numeric(gd$psrf[, "Point est."])
  }, error=function(e) NULL)
  rhat_names <- tryCatch({
    gd <- coda::gelman.diag(samp[, keep, drop=FALSE], multivariate=FALSE)
    rownames(gd$psrf)
  }, error=function(e) NULL)

  for(i in seq_along(keep)){
    p <- keep[i]
    v <- M[, p]
    out[out$param==p, "mean"] <- mean(v)
    out[out$param==p, "sd"]   <- sd(v)
    qs <- quantile(v, c(0.05,0.5,0.95))
    out[out$param==p, c("q05","q50","q95")] <- as.numeric(qs)

    if(!is.null(ess_all) && p %in% names(ess_all)){
      out[out$param==p, "ESS"] <- as.numeric(ess_all[p])
    }
    if(!is.null(rhat_all) && !is.null(rhat_names) && p %in% rhat_names){
      out[out$param==p, "Rhat"] <- rhat_all[match(p, rhat_names)]
    }
  }

  out
}

save_trace_pdf <- function(samp, params, file, w=10.5, h=7.5){
  params <- params[params %in% colnames(as.matrix(samp))]
  if(length(params)==0) return(invisible(NULL))
  tmp <- open_pdf_tmp(file, width=w, height=h)
  try({
    par(mfrow=c(ceiling(length(params)/2), 2), mar=c(3,3,2,1))
    for(p in params){
      try(coda::traceplot(samp[, p, drop=FALSE], main=p), silent=TRUE)
    }
  }, silent=TRUE)
  close_pdf_tmp(tmp, file)
}

# ------------------------------------------------------------
# 9) Core analysis function for one organism group and one model type
# ------------------------------------------------------------
# Main wrapper for one organismal group. It reads the species table, calculates response variables,
# constructs predictors, fits the JAGS model, and exports all outputs.
run_one_group <- function(prefix, spe_file, env_num, model_type){

  stopifnot(model_type %in% MODEL_TYPES)

  out_dir <- paste0("OUT_", prefix, "_", model_type)
  if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

  out_files <- list(
    csv_data     = file.path(out_dir, paste0(prefix, "_", model_type, "_alpha_beta_env_all.csv")),
    pdf_grid     = file.path(out_dir, paste0(prefix, "_", model_type, "_1_Posterior_Grid_3responses.pdf")),
    pdf_hop      = file.path(out_dir, paste0(prefix, "_", model_type, "_2_HOP_", focus_main, "_3responses.pdf")),
    pdf_pdp      = file.path(out_dir, paste0(prefix, "_", model_type, "_3_PDP_", focus_main, "_x_Other_3responses.pdf")),
    pdf_pdp_int  = file.path(out_dir, paste0(prefix, "_", model_type, "_3b_PDP_", focus_main, "_x_Other_withInt_3responses.pdf")),
    csv_sum_dir  = file.path(out_dir, paste0(prefix, "_", model_type, "_posterior_summary_directional_3responses.csv")),
    csv_runlog   = file.path(out_dir, paste0(prefix, "_", model_type, "_RUNLOG_3responses.csv")),
    csv_diag_a   = file.path(out_dir, paste0(prefix, "_", model_type, "_DIAG_3responses.csv")),
    pdf_trace_a  = file.path(out_dir, paste0(prefix, "_", model_type, "_TRACE_3responses.pdf"))
  )

  tryCatch({

    message("\n==============================")
    message("RUN GROUP: ", prefix, " | MODEL=", model_type, " | file=", spe_file)
    message("==============================")

    if(!file.exists(spe_file)) stop("spe file not found: ", spe_file)

    spe <- read.delim(spe_file, sep="\t", check.names=FALSE, stringsAsFactors=FALSE)
    if(ncol(spe) < 3) stop("spe file too few columns.")
    rownames(spe) <- spe[[1]]
    spe <- spe[, -1, drop=FALSE]

    sites <- intersect(colnames(spe), env_num$site)
    if(length(sites) < 10) stop("overlapping sites too few (<10).")

    spe2 <- spe[, sites, drop=FALSE]
    env2 <- env_num[match(sites, env_num$site), , drop=FALSE]

    # =========================================================
    # Part A: Calculate alpha-diversity metrics and site-level beta_local.
    # =========================================================
    X_asv <- t(spe2)
    richness <- vegan::specnumber(X_asv)
    shannon  <- vegan::diversity(X_asv, index="shannon")
    evenness <- ifelse(richness >= 2, shannon / log(richness), NA_real_)

    beta_site <- compute_site_beta_local(
      X_asv = X_asv,
      env2 = env2,
      comb_size = BETA_COMB_SIZE,
      use_relative_abundance = BETA_USE_RELATIVE_ABUNDANCE
    )

    alpha_df <- data.frame(
      site = rownames(X_asv),
      richness = as.numeric(richness),
      shannon  = as.numeric(shannon),
      evenness = as.numeric(evenness),
      total_reads = as.numeric(rowSums(X_asv)),
      stringsAsFactors = FALSE
    )

    dat <- join_by_site_keep_order(alpha_df, env2)
    dat$beta_local <- beta_site$beta_local[match(dat$site, beta_site$site)]
    dat$n_beta_combs <- beta_site$n_beta_combs[match(dat$site, beta_site$site)]

    if(CONTROL_READS_ALPHA){
      dat$log_total_reads <- log(pmax(dat$total_reads, 1))
    }

    write.csv(dat, out_files$csv_data, row.names=FALSE)

    exclude_cols <- c("site","richness","shannon","evenness","beta_local","n_beta_combs","total_reads")
    if(CONTROL_READS_ALPHA) exclude_cols <- c(exclude_cols, "log_total_reads")
    if(!is.null(group_col)) exclude_cols <- c(exclude_cols, group_col)
    exclude_cols <- unique(c(exclude_cols, "river"))

    preds0 <- setdiff(names(dat), exclude_cols)
    preds0 <- preds0[sapply(dat[preds0], is.numeric)]
    preds0 <- preds0[sapply(preds0, function(nm){
      x <- dat[[nm]]
      x <- x[is.finite(x)]
      length(x) >= 10 && sd(x, na.rm=TRUE) > 1e-8
    })]

    force_in <- intersect(force_pred_terms, names(dat))
    preds_all <- filter_predictors_by_model(preds0, force_in, model_type)
    if(CONTROL_READS_ALPHA) preds_all <- unique(c(preds_all, "log_total_reads"))

    preds_all <- preds_all[sapply(preds_all, function(nm){
      x <- dat[[nm]]
      x <- x[is.finite(x)]
      length(x) >= 10 && sd(x, na.rm=TRUE) > 1e-8
    })]

    if(!("FI" %in% preds_all)) stop("FI not in predictors after filtering (unexpected).")
    if(length(preds_all) < 1) stop("Need >=1 valid numeric predictors.")

    even_clamp <- pmin(pmax(dat$evenness, 0.001), 0.999)
    beta_clamp <- pmin(pmax(dat$beta_local, 0.001), 0.999)
    logit_even <- qlogis(even_clamp)
    logit_beta <- qlogis(beta_clamp)

    corr_score <- rep(0, length(preds_all)); names(corr_score) <- preds_all
    for(v in preds_all){
      x <- dat[[v]]
      ok <- is.finite(x) & is.finite(dat$shannon) & is.finite(logit_even) & is.finite(logit_beta)
      if(sum(ok) >= 12){
        c1 <- suppressWarnings(cor(x[ok], dat$shannon[ok]))
        c2 <- suppressWarnings(cor(x[ok], logit_even[ok]))
        c3 <- suppressWarnings(cor(x[ok], logit_beta[ok]))
        c1[!is.finite(c1)] <- 0; c2[!is.finite(c2)] <- 0; c3[!is.finite(c3)] <- 0
        corr_score[v] <- mean(abs(c(c1,c2,c3)))
      }
    }

    must_keep <- unique(c(force_in, if(CONTROL_READS_ALPHA) "log_total_reads" else character(0)))
    must_keep <- intersect(must_keep, preds_all)

    others <- setdiff(preds_all, must_keep)
    others <- others[order(corr_score[others], decreasing=TRUE)]

    n_left <- max(0, max_alpha_predictors - length(must_keep))
    preds_all_final <- unique(c(must_keep, head(others, n_left)))

    message(prefix, " [", model_type, "] Predictors used (cap=", max_alpha_predictors, "): ",
            paste(preds_all_final, collapse=", "))

    impute_log_a <- data.frame()
    if(model_type=="DIRECT" && IMPUTE_ENV_DIRECT){
      imp <- median_impute_numeric(dat, preds_all_final)
      dat <- imp$df
      impute_log_a <- imp$log
    }
    if(model_type=="TOTAL" && IMPUTE_ENV_TOTAL){
      imp <- median_impute_numeric(dat, preds_all_final)
      dat <- imp$df
      impute_log_a <- imp$log
    }

    Z_mat_a <- scale(dat[, preds_all_final, drop=FALSE])
    Z_center_a <- attr(Z_mat_a,"scaled:center")
    Z_scale_a  <- attr(Z_mat_a,"scaled:scale")
    Z_a <- as.data.frame(Z_mat_a)
    colnames(Z_a) <- paste0(preds_all_final, "_z")

    dat_z <- cbind(dat, Z_a)
    dat_z$row_id <- seq_len(nrow(dat_z))
    dat_z$shannon_y <- dat_z$shannon
    dat_z$evenness_clamp <- pmin(pmax(dat_z$evenness, 0.001), 0.999)
    dat_z$beta_clamp <- pmin(pmax(dat_z$beta_local, 0.001), 0.999)
    dat_z$logit_even <- qlogis(dat_z$evenness_clamp)
    dat_z$logit_beta <- qlogis(dat_z$beta_clamp)

    Z_names_a <- colnames(Z_a)

    need_cols_a <- c("row_id","river","n_beta_combs","shannon_y","logit_even","logit_beta", Z_names_a)
    if(!is.null(group_col)) need_cols_a <- c(need_cols_a, group_col)

    missA <- setdiff(need_cols_a, names(dat_z))
    if(length(missA)>0) stop("Missing columns: ", paste(missA, collapse=", "))

    df_fit_a <- dat_z[, unique(need_cols_a), drop=FALSE]
    if("n_beta_combs" %in% names(df_fit_a)) df_fit_a <- df_fit_a[df_fit_a$n_beta_combs >= BETA_MIN_COMBS_PER_SITE, , drop=FALSE]
    df_fit_a <- df_fit_a[stats::complete.cases(df_fit_a), , drop=FALSE]
    N_a <- nrow(df_fit_a)
    if(N_a < min_N_complete_sites) stop("Complete-case N < ", min_N_complete_sites)

    X0_a <- as.matrix(df_fit_a[, Z_names_a, drop=FALSE])
    focus_set_a <- paste0(intersect(force_pred_terms, preds_all_final), "_z")
    flt_a <- collinearity_filter(X0_a, cor_cut=cor_cut, focus_set=focus_set_a)
    Z_names_a <- flt_a$keep
    drop_a <- flt_a$drop

    if(!(paste0(focus_main, "_z") %in% Z_names_a)) stop("FI removed by collinearity filtering (unexpected).")

    Y_a <- as.matrix(df_fit_a[, c("shannon_y","logit_even","logit_beta"), drop=FALSE])
    Xmat_a <- as.matrix(df_fit_a[, Z_names_a, drop=FALSE])
    P_a <- ncol(Xmat_a)

    river_factor <- as.factor(df_fit_a$river)
    G_river <- nlevels(river_factor)
    river_id <- as.integer(river_factor)

    w_dir_a <- matrix(rep(c(0.25,0.50,0.25), P_a), nrow=P_a, byrow=TRUE)
    colnames(w_dir_a) <- c("neg","neu","pos")
    rownames(w_dir_a) <- Z_names_a

    idx_fi <- rownames(w_dir_a) == paste0(focus_main, "_z")
    if(any(idx_fi)) w_dir_a[idx_fi,] <- matrix(c(0.20,0.30,0.50), nrow=sum(idx_fi), ncol=3, byrow=TRUE)

    # JAGS model: three responses share the same predictor matrix but have response-specific coefficients.
# River-level random intercepts control non-independence among sites within the same river.
model_alpha <- "
model{
  mu_dir[1] <- -m0
  mu_dir[2] <-  0
  mu_dir[3] <-  m0
  tau_dir[1] <- 1/(sd_neg*sd_neg)
  tau_dir[2] <- 1/(sd_neu*sd_neu)
  tau_dir[3] <- 1/(sd_pos*sd_pos)

  for(j in 1:P){
    dir[j] ~ dcat(w_dir[j,1:3])
    mu_beta[j] ~ dnorm(mu_dir[dir[j]], tau_dir[dir[j]])
  }

  for(k in 1:3){
    alpha[k] ~ dnorm(0, 1.0E-4)
    tau_y[k] ~ dgamma(0.001, 0.001)
    sigma_y[k] <- 1/sqrt(tau_y[k])
    tau_river[k] ~ dgamma(0.5, 0.5)

    for(g in 1:G_river){
      a_river[g,k] ~ dnorm(0, tau_river[k])
    }

    for(j in 1:P){
      beta[j,k] ~ dnorm(mu_beta[j], tau_beta)
    }

    for(i in 1:N){
      mu[i,k] <- alpha[k] + inprod(X[i,1:P], beta[1:P,k]) + a_river[river_id[i],k]
      Y[i,k] ~ dnorm(mu[i,k], tau_y[k])
    }
  }
}
"

    jags_data_a <- list(
      N=N_a, P=P_a, G_river=G_river,
      Y=Y_a, X=Xmat_a, w_dir=w_dir_a, river_id=river_id,
      m0=0.25, sd_neg=0.25, sd_neu=0.70, sd_pos=0.25,
      tau_beta=1/(0.35^2)
    )

    message(prefix, " Running 3-response model (", model_type, ") ...")
    m_a <- jags.model(textConnection(model_alpha), data=jags_data_a,
                      n.chains=mcmc_alpha$n_chains, n.adapt=mcmc_alpha$n_adapt, quiet=TRUE)
    update(m_a, mcmc_alpha$n_update)
    samp_dir_a <- coda.samples(m_a, variable.names=c("alpha","beta","sigma_y"),
                               n.iter=mcmc_alpha$n_iter, thin=mcmc_alpha$thin)

    dir_mat_a <- as.matrix(samp_dir_a)
    sum_dir_a_all <- get_beta_summary_multi(samp_dir_a, Z_names_a)

    if(nrow(sum_dir_a_all)>0){
      sum_dir_a_focus <- sum_dir_a_all[sum_dir_a_all$predictor %in% paste0(focus_report_terms, "_z"), , drop=FALSE]
      write.csv(sum_dir_a_focus, out_files$csv_sum_dir, row.names=FALSE)
    } else {
      sum_dir_a_focus <- data.frame()
      write.csv(sum_dir_a_focus, out_files$csv_sum_dir, row.names=FALSE)
    }

    keyA <- c("alpha[1]","alpha[2]","alpha[3]","sigma_y[1]","sigma_y[2]","sigma_y[3]")
    for(tt in focus_report_terms){
      nmz <- paste0(tt, "_z")
      if(nmz %in% Z_names_a){
        j <- which(Z_names_a==nmz)
        keyA <- c(keyA,
                  paste0("beta[", j, ",1]"),
                  paste0("beta[", j, ",2]"),
                  paste0("beta[", j, ",3]"))
      }
    }
    keyA <- unique(keyA)

    diagA <- mcmc_diag_key(samp_dir_a, keyA)
    if(SAVE_DIAG_CSV) write.csv(diagA, out_files$csv_diag_a, row.names=FALSE)
    if(SAVE_TRACE_PDF) save_trace_pdf(samp_dir_a, keyA, out_files$pdf_trace_a)

    diag_flag <- function(diag_df){
      if(is.null(diag_df) || nrow(diag_df)==0) return("NO_DIAG")
      bad_rhat <- any(is.finite(diag_df$Rhat) & diag_df$Rhat > Rhat_threshold, na.rm=TRUE)
      bad_ess  <- any(is.finite(diag_df$ESS)  & diag_df$ESS < min_ESS_key, na.rm=TRUE)
      if(bad_rhat && bad_ess) return("BAD_RHAT_BAD_ESS")
      if(bad_rhat) return("BAD_RHAT")
      if(bad_ess)  return("BAD_ESS")
      "OK"
    }

    runlog <- data.frame(
      group=prefix,
      model_type=model_type,
      N_sites_model=N_a,
      beta_comb_size=BETA_COMB_SIZE,
      predictors=paste(preds_all_final, collapse=";"),
      Z_keep=paste(Z_names_a, collapse=";"),
      Z_drop=paste(drop_a, collapse=";"),
      impute_n = ifelse(is.null(impute_log_a) || nrow(impute_log_a)==0, 0, sum(impute_log_a$n_imputed)),
      diag_flag = diag_flag(diagA),
      stringsAsFactors=FALSE
    )
    write.csv(runlog, out_files$csv_runlog, row.names=FALSE)

    # =========================================================
    # Part B: Export per-group PDFs for Shannon, evenness, and beta_local.
    # =========================================================
    get_draws_coef <- function(pred_raw, resp_idx){
      pred_z <- paste0(pred_raw, "_z")
      if(!pred_z %in% Z_names_a) return(NULL)
      j <- which(Z_names_a==pred_z)
      nm <- paste0("beta[", j, ",", resp_idx, "]")
      if(!nm %in% colnames(dir_mat_a)) return(NULL)
      dir_mat_a[, nm]
    }

    grobs <- list()
    for(tt in focus_report_terms){
      grobs[[length(grobs)+1]] <- posterior_density_gg(get_draws_coef(tt, 1))
      grobs[[length(grobs)+1]] <- posterior_density_gg(get_draws_coef(tt, 2))
      grobs[[length(grobs)+1]] <- posterior_density_gg(get_draws_coef(tt, 3))
    }

    row_labels <- focus_report_terms
    col_labels_grid <- c("Shannon","Evenness","Beta_local")

    make_overview_grid(
      grobs=grobs,
      row_labels=row_labels,
      col_labels=col_labels_grid,
      title=paste0(prefix, " Posterior Grid (", model_type, ")"),
      file=out_files$pdf_grid,
      n_rows=length(row_labels), n_cols=3,
      page_w=13.6, page_h=3.1 + 2.55*length(row_labels)
    )

    fitted_rows <- df_fit_a$row_id
    draw_hop_shan <- get_draws_coef(focus_main, 1)
    draw_hop_even <- get_draws_coef(focus_main, 2)
    draw_hop_beta <- get_draws_coef(focus_main, 3)

    stat_hop_shan <- make_pd_hop_text(draw_hop_shan, focus_main)
    stat_hop_even <- make_pd_hop_text(draw_hop_even, focus_main)
    stat_hop_beta <- make_pd_hop_text(draw_hop_beta, focus_main)

    p_hop_shan <- hop_multi(focus_main, "shannon", dir_mat_a, Z_names_a, Z_center_a, Z_scale_a,
                            dat_z, fitted_rows, n_lines=n_lines_show, stat_text=stat_hop_shan)
    p_hop_even <- hop_multi(focus_main, "evenness", dir_mat_a, Z_names_a, Z_center_a, Z_scale_a,
                            dat_z, fitted_rows, n_lines=n_lines_show, stat_text=stat_hop_even)
    p_hop_beta <- hop_multi(focus_main, "beta_local", dir_mat_a, Z_names_a, Z_center_a, Z_scale_a,
                            dat_z, fitted_rows, n_lines=n_lines_show, stat_text=stat_hop_beta)

    hop_subtitle <- paste0(
      "Shannon: ", stat_hop_shan,
      "    | Evenness: ", stat_hop_even,
      "    | Beta_local: ", stat_hop_beta
    )

    tmp <- open_pdf_tmp(out_files$pdf_hop, width=18.0, height=7.8)
    grid::grid.newpage()
    grid::grid.draw(gridExtra::arrangeGrob(
      gridExtra::arrangeGrob(
        grid::textGrob(
          paste0(prefix, " HOP (", model_type, "): FI -> Shannon / Evenness / Beta_local"),
          gp=grid::gpar(fontsize=14, fontface="bold", fontfamily=PLOT_FAMILY)
        ),
        grid::textGrob(
          hop_subtitle,
          gp=grid::gpar(fontsize=10, col="grey25", fontfamily=PLOT_FAMILY)
        ),
        ncol=1, heights=c(0.60, 0.40)
      ),
      gridExtra::arrangeGrob(grobs=list(p_hop_shan, p_hop_even, p_hop_beta), ncol=3),
      ncol=1, heights=c(0.18, 0.82)
    ))
    close_pdf_tmp(tmp, out_files$pdf_hop)

    sumA <- sum_dir_a_all
    sumA$pred_raw <- sub("_z$", "", sumA$predictor)

    exclude_other <- unique(c(focus_main, int_terms, force_pred_terms, if(CONTROL_READS_ALPHA) "log_total_reads" else character(0)))
    candA <- unique(sumA$pred_raw[!(sumA$pred_raw %in% exclude_other)])

    other_alpha <- NA_character_
    if(length(candA)>0){
      scoreA <- sapply(candA, function(v){
        vv <- sumA$post_mean[sumA$pred_raw==v]
        mean(abs(vv), na.rm=TRUE)
      })
      candA <- candA[order(scoreA, decreasing=TRUE)]
      other_alpha <- candA[1]
    }

    p_pdp_shan <- if(is.na(other_alpha)) safe_blank_plot("NO OTHER") else
      pdp_multi_simple(focus_main, other_alpha, "shannon", dir_mat_a, Z_names_a, Z_center_a, Z_scale_a, dat_z, fitted_rows, n_grid=pdp_grid)
    p_pdp_even <- if(is.na(other_alpha)) safe_blank_plot("NO OTHER") else
      pdp_multi_simple(focus_main, other_alpha, "evenness", dir_mat_a, Z_names_a, Z_center_a, Z_scale_a, dat_z, fitted_rows, n_grid=pdp_grid)
    p_pdp_beta <- if(is.na(other_alpha)) safe_blank_plot("NO OTHER") else
      pdp_multi_simple(focus_main, other_alpha, "beta_local", dir_mat_a, Z_names_a, Z_center_a, Z_scale_a, dat_z, fitted_rows, n_grid=pdp_grid)

    tmp <- open_pdf_tmp(out_files$pdf_pdp, width=18.0, height=7.5)
    grid::grid.newpage()
    grid::grid.draw(gridExtra::arrangeGrob(
      grid::textGrob(
        paste0(prefix, " PDP (", model_type, "): FI × other = ", other_alpha),
        gp=grid::gpar(fontsize=14, fontface="bold", fontfamily=PLOT_FAMILY)
      ),
      gridExtra::arrangeGrob(grobs=list(p_pdp_shan, p_pdp_even, p_pdp_beta), ncol=3),
      ncol=1, heights=c(0.15, 0.85)
    ))
    close_pdf_tmp(tmp, out_files$pdf_pdp)

    tmp <- open_pdf_tmp(out_files$pdf_pdp_int, width=18.0, height=7.5)
    grid::grid.newpage()
    grid::grid.draw(gridExtra::arrangeGrob(
      grid::textGrob(
        paste0(prefix, " PDP(withInt) (", model_type, "): FI × other = ", other_alpha),
        gp=grid::gpar(fontsize=14, fontface="bold", fontfamily=PLOT_FAMILY)
      ),
      gridExtra::arrangeGrob(grobs=list(p_pdp_shan, p_pdp_even, p_pdp_beta), ncol=3),
      ncol=1, heights=c(0.15, 0.85)
    ))
    close_pdf_tmp(tmp, out_files$pdf_pdp_int)

    list(
      ok=TRUE,
      group=prefix,
      model_type=model_type,
      out_files=out_files,
      plots=list(
        hop=list(shannon=p_hop_shan, evenness=p_hop_even, beta_local=p_hop_beta),
        pdp=list(shannon=p_pdp_shan, evenness=p_pdp_even, beta_local=p_pdp_beta),
        other=list(alpha=other_alpha)
      ),
      sum_dir_a_focus=sum_dir_a_focus
    )

  }, error=function(e){
    list(ok=FALSE, group=prefix, model_type=model_type, error=conditionMessage(e), out_files=out_files)
  })
}

# ------------------------------------------------------------
# 10) Run all organism groups under TOTAL and DIRECT model sets
# ------------------------------------------------------------
groups <- list(
  list(prefix="bacteria",  file="spe_bacteria.txt"),
  list(prefix="protists",  file="spe_protist.txt"),
  list(prefix="metazoans", file="spe_Metazoa.txt"),
  list(prefix="fishes",    file="spe_fish.txt")
)

results <- list()
for(mt in MODEL_TYPES){
  results[[mt]] <- list()
  for(g in groups){
    results[[mt]][[g$prefix]] <- run_one_group(g$prefix, g$file, env_num, mt)
  }

  fail_df <- do.call(rbind, lapply(names(results[[mt]]), function(nm){
    obj <- results[[mt]][[nm]]
    if(is.null(obj) || isFALSE(obj$ok)){
      data.frame(model_type=mt, group=nm, error=ifelse(is.null(obj$error), "unknown", obj$error),
                 stringsAsFactors=FALSE)
    } else NULL
  }))
  if(!is.null(fail_df) && nrow(fail_df)>0){
    write.csv(fail_df, paste0("failed_groups_log_", mt, ".csv"), row.names=FALSE)
  }
}

# ------------------------------------------------------------
# 11) Export overview PDFs for each model type
# ------------------------------------------------------------
group_names <- sapply(groups, `[[`, "prefix")
row_metrics <- c("Shannon","Evenness","Beta_local")

for(mt in MODEL_TYPES){

  hop_grobs <- list()
  for(rm in row_metrics){
    key <- tolower(rm)
    for(gname in group_names){
      obj <- results[[mt]][[gname]]
      if(isTRUE(obj$ok)){
        p <- obj$plots$hop[[key]]
        if(is.null(p)) p <- safe_blank_plot(paste0(gname, "\nHOP missing"))
      } else {
        p <- safe_blank_plot(paste0(gname, "\nFAILED"))
      }
      if(STRIP_OVERVIEW_TITLES) p <- strip_title(p)
      hop_grobs[[length(hop_grobs)+1]] <- p
    }
  }
  make_overview_grid(
    grobs=hop_grobs,
    row_labels=row_metrics,
    col_labels=group_names,
    title=paste0("HOP (", mt, "): ", focus_main, " (rows=Shannon/Evenness/Beta_local ; cols=4 groups)"),
    file=paste0("overview_HOP_", focus_main, "_3responses_", mt, ".pdf"),
    n_rows=3, n_cols=length(group_names),
    page_w=OV_PAGE_W, page_h=OV_PAGE_H
  )

  pdp_grobs <- list()
  for(rm in row_metrics){
    key <- tolower(rm)
    for(gname in group_names){
      obj <- results[[mt]][[gname]]
      if(isTRUE(obj$ok)){
        p <- obj$plots$pdp[[key]]
        if(is.null(p)) p <- safe_blank_plot(paste0(gname, "\nPDP missing"))
      } else {
        p <- safe_blank_plot(paste0(gname, "\nFAILED"))
      }
      if(STRIP_OVERVIEW_TITLES) p <- strip_title(p)
      pdp_grobs[[length(pdp_grobs)+1]] <- p
    }
  }
  make_overview_grid(
    grobs=pdp_grobs,
    row_labels=row_metrics,
    col_labels=group_names,
    title=paste0("PDP (", mt, "): ", focus_main, " × other (rows=Shannon/Evenness/Beta_local ; cols=4 groups)"),
    file=paste0("overview_PDP_", focus_main, "_x_other_3responses_", mt, ".pdf"),
    n_rows=3, n_cols=length(group_names),
    page_w=OV_PAGE_W, page_h=OV_PAGE_H
  )
}

# ------------------------------------------------------------
# 12) Export cross-group posterior summaries and compare TOTAL vs DIRECT models
# ------------------------------------------------------------
collect_focus_summary <- function(mt){
  out <- list()
  for(gname in group_names){
    obj <- results[[mt]][[gname]]
    if(isTRUE(obj$ok)){
      if(!is.null(obj$sum_dir_a_focus) && nrow(obj$sum_dir_a_focus)>0){
        tmp <- obj$sum_dir_a_focus
        tmp$group <- gname
        tmp$model_type <- mt
        out[[length(out)+1]] <- tmp
      }
    }
  }
  if(length(out)==0) return(data.frame())
  do.call(rbind, out)
}

sum_TOTAL  <- collect_focus_summary("TOTAL")
sum_DIRECT <- collect_focus_summary("DIRECT")
write.csv(sum_TOTAL,  "overview_focus_posterior_summary_all_TOTAL_3responses.csv",  row.names=FALSE)
write.csv(sum_DIRECT, "overview_focus_posterior_summary_all_DIRECT_3responses.csv", row.names=FALSE)

if(nrow(sum_TOTAL)>0 && nrow(sum_DIRECT)>0){
  key <- c("group","response","predictor")
  A <- sum_TOTAL
  B <- sum_DIRECT
  A$key <- do.call(paste, c(A[key], sep="|"))
  B$key <- do.call(paste, c(B[key], sep="|"))

  common <- intersect(A$key, B$key)
  A2 <- A[match(common, A$key), , drop=FALSE]
  B2 <- B[match(common, B$key), , drop=FALSE]

  comp <- data.frame(
    group=A2$group,
    response=A2$response,
    predictor=A2$predictor,
    post_mean_TOTAL=A2$post_mean,
    CI90_low_TOTAL=A2$CI90_low,
    CI90_high_TOTAL=A2$CI90_high,
    P_gt_0_TOTAL=A2$P_gt_0,
    post_mean_DIRECT=B2$post_mean,
    CI90_low_DIRECT=B2$CI90_low,
    CI90_high_DIRECT=B2$CI90_high,
    P_gt_0_DIRECT=B2$P_gt_0,
    delta_mean_DIRECT_minus_TOTAL = B2$post_mean - A2$post_mean,
    delta_Pgt0_DIRECT_minus_TOTAL = B2$P_gt_0   - A2$P_gt_0,
    stringsAsFactors=FALSE
  )
  write.csv(comp, "COMPARE_TOTAL_vs_DIRECT_focus_terms_3responses.csv", row.names=FALSE)
}

# ------------------------------------------------------------
# 13) Export master run metadata
# ------------------------------------------------------------
meta <- data.frame(
  run_mode=run_mode,
  mcmc_alpha_chains=mcmc_alpha$n_chains,
  mcmc_alpha_adapt=mcmc_alpha$n_adapt,
  mcmc_alpha_update=mcmc_alpha$n_update,
  mcmc_alpha_iter=mcmc_alpha$n_iter,
  mcmc_alpha_thin=mcmc_alpha$thin,
  max_predictors=max_alpha_predictors,
  cor_cut=cor_cut,
  IMPUTE_ENV_DIRECT=IMPUTE_ENV_DIRECT,
  CONTROL_READS_ALPHA=CONTROL_READS_ALPHA,
  BETA_COMB_SIZE=BETA_COMB_SIZE,
  BETA_USE_RELATIVE_ABUNDANCE=BETA_USE_RELATIVE_ABUNDANCE,
  stringsAsFactors=FALSE
)
write.csv(meta, "RUN_METADATA_SCI_3responses.csv", row.names=FALSE)

message("\nALL DONE ✅  (TOTAL + DIRECT, Shannon + Evenness + Beta_local).")
message("Note: beta_local is computed from all four-site combinations within each river.")
message("Note: diagnostics for each group/model are saved in OUT_* directories: DIAG_3responses.csv, TRACE_3responses.pdf, and RUNLOG_3responses.csv.")
