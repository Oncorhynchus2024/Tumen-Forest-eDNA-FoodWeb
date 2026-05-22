# ============================================================
# One-click pipeline for food-web network metrics and FI models
# TOTAL vs DIRECT Bayesian modelling with brms and optional glmmTMB fallback
#
# Purpose
#   This script converts genus-level abundance data into functional feeding
#   group (FFG) food webs, calculates site-level network metrics, and fits
#   Bayesian models to evaluate the association between the forest cover
#   gradient (FI) and food-web structure.
#
# Model types
#   TOTAL  : estimates FI effects while controlling for mainstem status and
#            river width; potential environmental mediators are excluded.
#   DIRECT : estimates conditional FI effects after including available
#            environmental covariates and, when possible, river random effects.
#
# Main inputs
#   abund.txt      : genus-by-site or site-by-genus abundance table.
#   genus_ffg.txt  : mapping between Genus and functional feeding group (FFG).
#   ffg_edges.txt  : directed FFG trophic links with FromFFG and ToFFG columns.
#   env.txt        : site-level environmental data, including FI and river.
#
# Main outputs
#   FI_FFG_TOTAL/ and FI_FFG_DIRECT/
#     - Network_metrics_by_site_FFG.tsv
#     - FI_model_HOP_table_raw.tsv
#     - diagnostic files, LOO files, fitted model RDS files
#     - FI partial-effect summary PDFs
#
# Notes
#   - The script uses presence/absence data by default.
#   - The edge direction is interpreted as prey -> predator.
#   - Genus_richness can be restricted to genera mapped to an FFG.
#   - Null-model metrics are optionally calculated for network comparison.
# ============================================================

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
set.seed(123)

# -------------------------
# 0) Packages
# -------------------------
pkgs <- c(
    "tidyverse", "janitor", "igraph", "patchwork",
    "brms", "posterior", "glmmTMB", "MASS",
    "rstan", "loo"
)

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(to_install) > 0) {
    install.packages(to_install, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(igraph)
    library(patchwork)
    library(brms)
    library(posterior)
    library(glmmTMB)
    library(MASS)
    library(rstan)
    library(loo)
})

# -------------------------
# 1) User settings & RUN MODE
# -------------------------

# Available run modes:
# "fast"         : quick testing and trend checking
# "final"        : publication-oriented run
# "final_strict" : stricter submission-level run
RUN_MODE <- "final_strict"

ABUND_FILE <- "abund.txt"
GENUS_FFG  <- "genus_ffg.txt"
FFG_EDGES  <- "ffg_edges.txt"
ENV_FILE   <- "env.txt"

MODEL_TYPES <- c("TOTAL", "DIRECT")

MIN_TOTAL_ABUND <- 1
USE_PRESENCE_ONLY <- TRUE

# River random effects are used only in DIRECT models.
USE_RIVER_RE_DIRECT <- TRUE

USE_LINK_DENSITY <- TRUE

CORES <- max(1, parallel::detectCores() - 1)
MAX_TDEPTH <- 15
EDGE_DIRECTION <- "prey_to_predator"

# Predator classification threshold based on trophic level.
PRED_TL_THRESHOLD <- 2.5
RUN_PRED_TL_SENSITIVITY <- TRUE
PRED_TL_SENS_THRESHOLDS <- c(2.5, 3.0, 3.5)

ALWAYS_INCLUDE_BASAL_RESOURCES <- TRUE
BASAL_ALWAYS_CAND <- c(
    "Plants","Fungi","Plankton",
    "AquaticDetritus","TerrestrialDetritus",
    "Aquatic_Detritus","Terrestrial_Detritus",
    "Detritus","Algae"
)

ADD_DISTANCE_METRICS <- TRUE
ADD_DIAMETER <- FALSE

RUN_NULL_MODEL <- TRUE
NULL_KEEP_BASAL_MODE <- "fixed_list"
BASAL_FIXED <- c(
    "Detritus","Algae","Plants","Fungi","Plankton",
    "AquaticDetritus","TerrestrialDetritus",
    "Terrestrial_Detritus","Aquatic_Detritus"
)
NULL_MAX_TRIES <- 10000

RESPONSE_MODE <- "raw"   # "raw" / "devnull" / "ses"
IMPUTE_MISSING_ENV <- TRUE

SAVE_DIAGNOSTICS <- TRUE
SAVE_PPC_PLOTS   <- FALSE
RUN_LOO_COMPARE  <- TRUE

# Whether genus richness should count only genera mapped to an FFG.
GENUS_RICHNESS_MAPPED_ONLY <- TRUE

# Fixed Stan seed for reproducibility.
STAN_SEED <- 20260307

# Strict diagnostic thresholds for Bayesian model quality.
STRICT_RHAT_MAX <- 1.01
STRICT_ESS_MIN  <- 400
STRICT_DIV_MAX  <- 0
STRICT_BFMI_MIN <- 0.30

if (RUN_MODE == "final_strict") {
    message(">>> Current run mode: FINAL_STRICT (submission-level settings)")
    FI_GRID_N       <- 300
    CHAINS          <- 4
    ITER            <- 20000
    WARMUP          <- 10000  # Warmup is set to half of the total iterations.
    ADAPT_DELTA     <- 0.9995
    MAX_TDEPTH      <- 15
    NDRAW_SPAGHETTI <- 800
    MODULARITY_N    <- 2000
    NULL_N          <- 9999
    ALLOW_FALLBACK_GLMMTMB <- FALSE
} else if (RUN_MODE == "final") {
    message(">>> Current run mode: FINAL (publication-oriented settings)")
    FI_GRID_N       <- 300
    CHAINS          <- 4
    ITER            <- 10000
    WARMUP          <- 5000
    ADAPT_DELTA     <- 0.999
    MAX_TDEPTH      <- 15
    NDRAW_SPAGHETTI <- 600
    MODULARITY_N    <- 1000
    NULL_N          <- 4999
    ALLOW_FALLBACK_GLMMTMB <- FALSE
} else {
    message(">>> Current run mode: FAST (quick testing)")
    FI_GRID_N       <- 100
    CHAINS          <- 4
    ITER            <- 4000
    WARMUP          <- 2000
    ADAPT_DELTA     <- 0.99
    MAX_TDEPTH      <- 12
    NDRAW_SPAGHETTI <- 300
    MODULARITY_N    <- 100
    NULL_N          <- 499
    ALLOW_FALLBACK_GLMMTMB <- TRUE
}

# -------------------------
# 2) Utilities
# -------------------------
scale_z <- function(x){
    x <- suppressWarnings(as.numeric(x))
    mu <- mean(x, na.rm = TRUE)
    sdv <- sd(x, na.rm = TRUE)
    if(!is.finite(sdv) || sdv == 0) return(list(z = x * 0, mu = mu, sd = 1))
    list(z = (x - mu) / sdv, mu = mu, sd = sdv)
}

squeeze01 <- function(y, n){
    y <- pmin(pmax(as.numeric(y), 0), 1)
    (y * (n - 1) + 0.5) / n
}

thr_tag <- function(thr){
    s <- as.character(thr)
    s <- sub("\\.0+$", "", s)
    s <- gsub("\\.", "p", s)
    s
}

skewness_fp <- function(x){
    x <- suppressWarnings(as.numeric(x))
    x <- x[is.finite(x)]
    n <- length(x)
    if(n < 3) return(NA_real_)
    m <- mean(x)
    s <- stats::sd(x)
    if(!is.finite(s) || s == 0) return(0)
    sum((x - m)^3) / n / (s^3)
}

match_ignore_case <- function(cand, pool){
    if(length(cand) == 0 || length(pool) == 0) return(character(0))
    pool[tolower(pool) %in% tolower(cand)]
}

build_diet_matrix <- function(nodes, ed_dir){
    nodes <- unique(as.character(nodes))
    if(nrow(ed_dir) < 1) return(list(M = NULL, consumers = character(0), resources = character(0)))
    
    consumers <- unique(as.character(ed_dir$ToFFG))
    resources <- unique(as.character(ed_dir$FromFFG))
    consumers <- consumers[consumers %in% nodes]
    resources <- resources[resources %in% nodes]
    
    if(length(consumers) < 1 || length(resources) < 1){
        return(list(M = NULL, consumers = consumers, resources = resources))
    }
    
    M <- matrix(0, nrow = length(consumers), ncol = length(resources),
                dimnames = list(consumers, resources))
    idx_c <- match(as.character(ed_dir$ToFFG), consumers)
    idx_r <- match(as.character(ed_dir$FromFFG), resources)
    ok <- is.finite(idx_c) & is.finite(idx_r)
    
    if(any(ok)){
        for(k in which(ok)){
            M[idx_c[k], idx_r[k]] <- 1
        }
    }
    list(M = M, consumers = consumers, resources = resources)
}

nodf01 <- function(M){
    if(is.null(M)) return(NA_real_)
    M <- (M > 0) * 1
    if(nrow(M) < 2 || ncol(M) < 2) return(NA_real_)
    
    rs <- rowSums(M); cs <- colSums(M)
    keep_r <- rs > 0; keep_c <- cs > 0
    M <- M[keep_r, keep_c, drop = FALSE]
    if(nrow(M) < 2 || ncol(M) < 2) return(NA_real_)
    
    rs <- rowSums(M); cs <- colSums(M)
    or <- order(rs, decreasing = TRUE)
    oc <- order(cs, decreasing = TRUE)
    M <- M[or, oc, drop = FALSE]
    rs <- rs[or]; cs <- cs[oc]
    
    acc_r <- 0; cnt_r <- 0
    for(i in 1:(nrow(M)-1)){
        for(j in (i+1):nrow(M)){
            if(rs[i] > rs[j] && rs[j] > 0){
                ov <- sum(M[i,] * M[j,])
                acc_r <- acc_r + (ov / rs[j])
                cnt_r <- cnt_r + 1
            }
        }
    }
    
    acc_c <- 0; cnt_c <- 0
    for(i in 1:(ncol(M)-1)){
        for(j in (i+1):ncol(M)){
            if(cs[i] > cs[j] && cs[j] > 0){
                ov <- sum(M[,i] * M[,j])
                acc_c <- acc_c + (ov / cs[j])
                cnt_c <- cnt_c + 1
            }
        }
    }
    
    denom <- cnt_r + cnt_c
    if(denom == 0) return(0)
    (acc_r + acc_c) / denom
}

mean_jaccard_diet <- function(M){
    if(is.null(M)) return(NA_real_)
    M <- (M > 0) * 1
    rs <- rowSums(M)
    M <- M[rs > 0, , drop = FALSE]
    if(nrow(M) < 2) return(NA_real_)
    
    acc <- 0; cnt <- 0
    for(i in 1:(nrow(M)-1)){
        for(j in (i+1):nrow(M)){
            inter <- sum(M[i,] * M[j,])
            uni <- sum((M[i,] + M[j,]) > 0)
            if(uni > 0){
                acc <- acc + inter / uni
                cnt <- cnt + 1
            }
        }
    }
    
    if(cnt == 0) return(NA_real_)
    acc / cnt
}

modularity_louvain_mean <- function(g_und, n = 100, seed_base = 2025){
    if(vcount(g_und) < 2 || ecount(g_und) < 1){
        return(list(mod_mean = NA_real_, nmod_mean = NA_real_))
    }
    mods <- numeric(n)
    nmods <- numeric(n)
    for(i in seq_len(n)){
        set.seed(seed_base + i)
        perm_vec <- sample.int(vcount(g_und))
        g2 <- igraph::permute(g_und, perm_vec)
        cl <- igraph::cluster_louvain(g2)
        mods[i] <- igraph::modularity(cl)
        nmods[i] <- length(unique(igraph::membership(cl)))
    }
    list(mod_mean = mean(mods, na.rm = TRUE), nmod_mean = mean(nmods, na.rm = TRUE))
}

louvain_partition_once <- function(g_und, seed = 2025){
    if(vcount(g_und) < 2 || ecount(g_und) < 1){
        mem <- rep(NA_integer_, vcount(g_und))
        names(mem) <- V(g_und)$name
        return(list(
            membership = mem,
            modularity = NA_real_,
            n_modules = NA_real_
        ))
    }
    set.seed(seed)
    cl <- igraph::cluster_louvain(g_und)
    mem <- igraph::membership(cl)
    mem <- mem[V(g_und)$name]
    list(
        membership = mem,
        modularity = igraph::modularity(cl),
        n_modules = length(unique(mem))
    )
}

participation_coefficient_und <- function(g_und, membership_vec){
    vnm <- V(g_und)$name
    out <- rep(NA_real_, length(vnm))
    names(out) <- vnm
    
    if(vcount(g_und) < 2 || ecount(g_und) < 1) return(out)
    if(length(membership_vec) < 1) return(out)
    
    membership_vec <- membership_vec[vnm]
    if(all(is.na(membership_vec))) return(out)
    
    A <- as.matrix(igraph::as_adjacency_matrix(g_und, sparse = FALSE))
    diag(A) <- 0
    mods <- sort(unique(as.numeric(membership_vec[is.finite(membership_vec)])))
    
    for(i in seq_along(vnm)){
        ki <- sum(A[i, ], na.rm = TRUE)
        if(!is.finite(ki) || ki <= 0){
            out[i] <- NA_real_
        } else {
            acc <- 0
            for(m in mods){
                idx_m <- which(membership_vec == m)
                kim <- sum(A[i, idx_m], na.rm = TRUE)
                acc <- acc + (kim / ki)^2
            }
            out[i] <- 1 - acc
        }
    }
    out
}

is_prop_metric <- function(m){
    if(grepl("^Prop_Predators_TL", m)) return(TRUE)
    if(grepl("^Prop_BasalConsumers_TL", m)) return(TRUE)
    m %in% c(
        "Connectance_dir","Prop_Basal",
        "Prop_Predators","Prop_BasalConsumers",
        "Prop_Omnivores",
        "Nestedness_NODF","Niche_overlap",
        "Inter_module_link_ratio",
        "Mean_PC_Consumers",
        "Mean_PC_TopPredators",
        "Prop_CrossModulePredators_k2",
        "Prop_CrossModulePredators_k3"
    )
}

is_count_metric <- function(m){
    m %in% c("S_nodes","L_edges","Genus_richness")
}

inv_link_mu <- function(eta, fam){
    if(fam == "beta")   return(plogis(eta))
    if(fam == "negbin") return(exp(eta))
    eta
}


has_valid_term <- function(x){
    if(is.null(x)) return(FALSE)
    if(is.factor(x) || is.character(x)){
        ux <- unique(as.character(x[!is.na(x)]))
        ux <- ux[ux != ""]
        return(length(ux) >= 2)
    }
    x <- suppressWarnings(as.numeric(x))
    return(sum(is.finite(x)) > 1 && length(unique(x[is.finite(x)])) >= 2)
}

ref_value_for_prediction <- function(x){
    if(is.factor(x)){
        tab <- table(x, useNA = "no")
        if(length(tab) < 1) return(factor(NA, levels = levels(x)))
        lev <- names(which.max(tab))[1]
        return(factor(lev, levels = levels(x)))
    }
    if(is.character(x)){
        tab <- sort(table(x, useNA = "no"), decreasing = TRUE)
        if(length(tab) < 1) return(NA_character_)
        return(names(tab)[1])
    }
    x_num <- suppressWarnings(as.numeric(x))
    x_num <- x_num[is.finite(x_num)]
    if(length(x_num) < 1) return(NA_real_)
    stats::median(x_num, na.rm = TRUE)
}

HAS_MAKE <- nzchar(Sys.which("make"))
message("Toolchain check: make = ", ifelse(HAS_MAKE, Sys.which("make"), "<NOT FOUND>"))
CAN_TRY_BRMS <- HAS_MAKE

open_pdf_tmp <- function(final_path, width, height){
    tmp_path <- paste0(final_path, ".tmp.pdf")
    if(file.exists(tmp_path)) suppressWarnings(file.remove(tmp_path))
    pdf(tmp_path, width = width, height = height, onefile = TRUE)
    return(tmp_path)
}

close_pdf_tmp <- function(tmp_path, final_path){
    dev.off()
    if(file.exists(final_path)) suppressWarnings(file.remove(final_path))
    ok <- file.rename(tmp_path, final_path)
    if(!isTRUE(ok)){
        warning("PDF overwrite failed, possibly because the target file is open. New file retained: ", tmp_path)
    }
}

safe_write_tsv <- function(x, file){
    write.table(x, file, sep = "\t", row.names = FALSE, quote = FALSE)
}

ebfmi_one <- function(energy){
    energy <- as.numeric(energy)
    energy <- energy[is.finite(energy)]
    if(length(energy) < 2) return(NA_real_)
    num <- sum(diff(energy)^2)
    den <- sum((energy - mean(energy))^2)
    if(!is.finite(den) || den <= 0) return(NA_real_)
    num / den
}

extract_draws_diag <- function(fit){
    s <- posterior::summarise_draws(as_draws_array(fit))
    s <- as.data.frame(s)
    s <- tibble::rownames_to_column(s, "param")
    s
}

extract_sampler_diag <- function(fit, max_treedepth){
    sp_list <- rstan::get_sampler_params(fit$fit, inc_warmup = FALSE)
    if(length(sp_list) < 1){
        return(data.frame(
            chain = integer(0), n_divergent = integer(0),
            n_max_treedepth = integer(0), ebfmi = numeric(0), stringsAsFactors = FALSE
        ))
    }
    out <- lapply(seq_along(sp_list), function(i){
        sp <- as.data.frame(sp_list[[i]])
        data.frame(
            chain = i,
            n_divergent = sum(sp$divergent__ > 0, na.rm = TRUE),
            n_max_treedepth = sum(sp$treedepth__ >= max_treedepth, na.rm = TRUE),
            ebfmi = ebfmi_one(sp$energy__),
            stringsAsFactors = FALSE
        )
    })
    bind_rows(out)
}

diagnose_brms_fit <- function(fit, metric, fam_tag, max_treedepth){
    ddraw <- extract_draws_diag(fit)
    dsamp <- extract_sampler_diag(fit, max_treedepth)
    
    max_rhat <- suppressWarnings(max(ddraw$rhat, na.rm = TRUE))
    min_ess_bulk <- suppressWarnings(min(ddraw$ess_bulk, na.rm = TRUE))
    min_ess_tail <- suppressWarnings(min(ddraw$ess_tail, na.rm = TRUE))
    
    if(!is.finite(max_rhat)) max_rhat <- NA_real_
    if(!is.finite(min_ess_bulk)) min_ess_bulk <- NA_real_
    if(!is.finite(min_ess_tail)) min_ess_tail <- NA_real_
    
    n_divergent_total <- if(nrow(dsamp) > 0) sum(dsamp$n_divergent, na.rm = TRUE) else NA_integer_
    n_max_treedepth_total <- if(nrow(dsamp) > 0) sum(dsamp$n_max_treedepth, na.rm = TRUE) else NA_integer_
    min_ebfmi <- if(nrow(dsamp) > 0) min(dsamp$ebfmi, na.rm = TRUE) else NA_real_
    
    pass_strict <- isTRUE(
        is.finite(max_rhat) && max_rhat < STRICT_RHAT_MAX &&
            is.finite(min_ess_bulk) && min_ess_bulk > STRICT_ESS_MIN &&
            is.finite(min_ess_tail) && min_ess_tail > STRICT_ESS_MIN &&
            is.finite(n_divergent_total) && n_divergent_total <= STRICT_DIV_MAX &&
            (is.na(min_ebfmi) || min_ebfmi >= STRICT_BFMI_MIN)
    )
    
    list(
        detail = ddraw,
        chain_diag = dsamp,
        summary = data.frame(
            Metric = metric, family = fam_tag, max_rhat = max_rhat,
            min_ess_bulk = min_ess_bulk, min_ess_tail = min_ess_tail,
            n_divergent = n_divergent_total, n_max_treedepth = n_max_treedepth_total,
            min_ebfmi = min_ebfmi, pass_strict = pass_strict, stringsAsFactors = FALSE
        )
    )
}

save_loo_outputs <- function(fit, metric, out_prefix){
    loo_res <- tryCatch({
        if(is.null(fit$criteria$loo)){
            fit2 <- add_criterion(fit, criterion = "loo")
            list(fit = fit2, loo = fit2$criteria$loo)
        } else {
            list(fit = fit, loo = fit$criteria$loo)
        }
    }, error = function(e) NULL)
    
    if(is.null(loo_res)){
        return(list(
            fit = fit,
            ok = FALSE,
            summary = data.frame(Metric = metric, loo_ok = FALSE, stringsAsFactors = FALSE)
        ))
    }
    
    fit <- loo_res$fit
    loo_obj <- loo_res$loo
    est <- as.data.frame(loo_obj$estimates)
    est$measure <- rownames(est)
    rownames(est) <- NULL
    est <- est[, c("measure","Estimate","SE"), drop = FALSE]
    names(est) <- c("measure","estimate","se")
    est$Metric <- metric
    
    pk <- tryCatch(loo::pareto_k_values(loo_obj), error = function(e) NULL)
    if(!is.null(pk)){
        pareto_df <- data.frame(
            Metric = metric,
            n_pareto_k_gt_0.5 = sum(pk > 0.5, na.rm = TRUE),
            n_pareto_k_gt_0.7 = sum(pk > 0.7, na.rm = TRUE),
            n_pareto_k_gt_1.0 = sum(pk > 1.0, na.rm = TRUE),
            max_pareto_k = suppressWarnings(max(pk, na.rm = TRUE)),
            stringsAsFactors = FALSE
        )
    } else {
        pareto_df <- data.frame(
            Metric = metric,
            n_pareto_k_gt_0.5 = NA_integer_,
            n_pareto_k_gt_0.7 = NA_integer_,
            n_pareto_k_gt_1.0 = NA_integer_,
            max_pareto_k = NA_real_,
            stringsAsFactors = FALSE
        )
    }
    
    safe_write_tsv(est, paste0(out_prefix, "_loo_summary.tsv"))
    safe_write_tsv(pareto_df, paste0(out_prefix, "_loo_pareto_k.tsv"))
    
    sum1 <- est %>% pivot_wider(names_from = measure, values_from = c(estimate, se), names_sep = "_")
    out_sum <- bind_cols(
        data.frame(Metric = metric, loo_ok = TRUE, stringsAsFactors = FALSE),
        sum1,
        pareto_df[, setdiff(names(pareto_df), "Metric"), drop = FALSE]
    )
    list(fit = fit, ok = TRUE, summary = out_sum)
}

make_priors <- function(fam_tag, has_river){
    p <- c(
        set_prior("student_t(3, 0, 2.5)", class = "Intercept"),
        set_prior("normal(0, 1)", class = "b")
    )
    if(has_river) p <- c(p, set_prior("exponential(1)", class = "sd"))
    if(fam_tag == "gaussian") p <- c(p, set_prior("exponential(1)", class = "sigma"))
    if(fam_tag == "negbin")   p <- c(p, set_prior("exponential(1)", class = "shape"))
    if(fam_tag == "beta")     p <- c(p, set_prior("exponential(1)", class = "phi"))
    p
}

# -------------------------
# 3) Load env.txt
# -------------------------
message("Loading env.txt ...")
env_raw <- read.delim(ENV_FILE, sep = "\t", header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(env_raw) <- trimws(names(env_raw))

if("PO₄" %in% names(env_raw)) names(env_raw)[names(env_raw) == "PO₄"] <- "PO4"
if("c" %in% names(env_raw)) names(env_raw)[names(env_raw) == "c"] <- "Cond"
if("pH" %in% names(env_raw)) names(env_raw)[names(env_raw) == "pH"] <- "ph"

env <- env_raw %>% clean_names()

if("p_h" %in% names(env)) names(env)[names(env) == "p_h"] <- "ph"
if("site" %in% names(env)) names(env)[names(env) == "site"] <- "Site"
if("fi" %in% names(env))   names(env)[names(env) == "fi"] <- "FI"

env$Site <- as.character(env$Site)
env$river <- if("river" %in% names(env)) as.character(env$river) else NA_character_

for(nm in names(env)){
    if(nm %in% c("Site","river")) next
    env[[nm]] <- suppressWarnings(as.numeric(env[[nm]]))
}

# -------------------------
# 4) Load genus_ffg.txt
# -------------------------
message("Loading genus_ffg.txt ...")
gffg_raw <- read.delim(GENUS_FFG, sep = "\t", header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(gffg_raw) <- trimws(names(gffg_raw))
gffg <- gffg_raw %>% clean_names()

gffg <- gffg %>%
    transmute(Genus = as.character(genus), FFG = as.character(ffg))

if("separate_rows" %in% getNamespaceExports("tidyr")){
    gffg <- gffg %>%
        tidyr::separate_rows(FFG, sep = ";") %>%
        mutate(FFG = stringr::str_squish(FFG)) %>%
        filter(!is.na(Genus), !is.na(FFG), Genus != "", FFG != "") %>%
        distinct()
} else {
    spl <- strsplit(gffg$FFG, ";", fixed = TRUE)
    out <- lapply(seq_along(spl), function(i){
        data.frame(Genus = gffg$Genus[i], FFG = trimws(spl[[i]]), stringsAsFactors = FALSE)
    })
    gffg <- bind_rows(out) %>%
        mutate(FFG = stringr::str_squish(FFG)) %>%
        filter(!is.na(Genus), !is.na(FFG), Genus != "", FFG != "") %>%
        distinct()
}

# -------------------------
# 5) Load ffg_edges.txt
# -------------------------
message("Loading ffg_edges.txt ...")
ffg_edges_raw <- read.delim(FFG_EDGES, sep = "\t", header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(ffg_edges_raw) <- trimws(names(ffg_edges_raw))
names(ffg_edges_raw) <- janitor::make_clean_names(names(ffg_edges_raw))

from_cand <- grep("^from(_)?ffg$", names(ffg_edges_raw), value = TRUE)
to_cand   <- grep("^to(_)?ffg$", names(ffg_edges_raw), value = TRUE)
if(length(from_cand) < 1 || length(to_cand) < 1){
    stop("Could not identify FromFFG/ToFFG columns in ffg_edges.txt.")
}

ffg_edges <- data.frame(
    FromFFG = as.character(ffg_edges_raw[[from_cand[1]]]),
    ToFFG   = as.character(ffg_edges_raw[[to_cand[1]]]),
    stringsAsFactors = FALSE
) %>%
    mutate(FromFFG = trimws(FromFFG), ToFFG = trimws(ToFFG)) %>%
    filter(!is.na(FromFFG), !is.na(ToFFG), FromFFG != "", ToFFG != "") %>%
    distinct() %>% filter(FromFFG != ToFFG)

META_NODES <- unique(c(ffg_edges$FromFFG, ffg_edges$ToFFG))
BASAL_ALWAYS <- character(0)
if(isTRUE(ALWAYS_INCLUDE_BASAL_RESOURCES)){
    BASAL_ALWAYS <- match_ignore_case(BASAL_ALWAYS_CAND, META_NODES)
}

# -------------------------
# 6) Load abund.txt and reshape
# -------------------------
message("Loading abund.txt ...")
abund_raw <- read.delim(ABUND_FILE, sep = "\t", header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(abund_raw) <- trimws(names(abund_raw))

abund_raw <- abund_raw %>%
    dplyr::select(-any_of(c("tax", "Tax", "TAX", "taxonomy", "Taxonomy")))

env_sites <- env$Site

if(sum(colnames(abund_raw) %in% env_sites) >= sum(as.character(abund_raw[[1]]) %in% env_sites)){
    names(abund_raw)[1] <- "Genus"
    cols_to_pivot <- intersect(names(abund_raw)[-1], env_sites)
    abund_long <- tidyr::pivot_longer(
        abund_raw, cols = all_of(cols_to_pivot), names_to = "Site", values_to = "Abund"
    )
} else {
    names(abund_raw)[1] <- "Site"
    cols_to_pivot <- setdiff(names(abund_raw), "Site")
    abund_long <- tidyr::pivot_longer(
        abund_raw, cols = all_of(cols_to_pivot), names_to = "Genus", values_to = "Abund"
    )
}

abund_long <- abund_long %>%
    mutate(
        Genus = as.character(Genus),
        Site  = as.character(Site),
        Abund = suppressWarnings(as.numeric(as.character(Abund)))
    ) %>%
    filter(Site %in% env_sites, is.finite(Abund)) %>%
    group_by(Site, Genus) %>%
    summarise(Abund = sum(Abund, na.rm = TRUE), .groups = "drop") %>%
    filter(Abund >= MIN_TOTAL_ABUND)

if(USE_PRESENCE_ONLY) abund_long$Abund <- 1

# Genus richness calculation.
if(GENUS_RICHNESS_MAPPED_ONLY){
    genus_pool_for_richness <- abund_long %>%
        inner_join(gffg %>% distinct(Genus), by = "Genus") %>%
        distinct(Site, Genus)
} else {
    genus_pool_for_richness <- abund_long %>%
        distinct(Site, Genus)
}

site_genus_richness <- genus_pool_for_richness %>%
    count(Site, name = "Genus_richness")

site_ffg <- abund_long %>%
    inner_join(gffg, by = "Genus") %>%
    distinct(Site, FFG)

# -------------------------
# 7) Build per-site FFG networks & metrics
# -------------------------
message("Building per-site FFG networks & metrics ...")

calc_trophic_level <- function(g_dir){
    v <- V(g_dir)$name
    S <- length(v)
    if(S == 0) return(setNames(numeric(0), character(0)))
    
    A <- as.matrix(igraph::as_adjacency_matrix(g_dir, sparse = FALSE))
    indeg <- colSums(A)
    basal <- indeg == 0
    
    W <- matrix(0, nrow = S, ncol = S, dimnames = list(v, v))
    for(i in seq_len(S)){
        ki <- indeg[i]
        if(ki > 0){
            prey_idx <- which(A[, i] > 0)
            W[i, prey_idx] <- 1 / ki
        }
    }
    
    TL <- as.numeric(MASS::ginv(diag(S) - W) %*% rep(1, S))
    names(TL) <- v
    TL[basal] <- 1
    TL
}

randomize_edges_simple <- function(nodes, L, basal_no_pred = character(0), max_tries = 3000){
    nodes <- unique(as.character(nodes))
    S <- length(nodes)
    if(S < 2 || L <= 0) return(data.frame(FromFFG = character(0), ToFFG = character(0)))
    
    predator_pool <- setdiff(nodes, basal_no_pred)
    if(length(predator_pool) < 1) predator_pool <- nodes
    
    possible <- expand.grid(FromFFG = nodes, ToFFG = predator_pool, stringsAsFactors = FALSE)
    possible <- possible[possible$FromFFG != possible$ToFFG, , drop = FALSE]
    if(nrow(possible) < 1) return(data.frame(FromFFG = character(0), ToFFG = character(0)))
    
    sel <- possible[0, , drop = FALSE]
    tries <- 0
    while(nrow(sel) < L && tries < max_tries){
        need <- L - nrow(sel)
        add_n <- min(max(50, need * 2), 5000)
        add <- possible[sample.int(nrow(possible), size = add_n, replace = TRUE), , drop = FALSE]
        sel <- dplyr::bind_rows(sel, add) %>% dplyr::distinct(FromFFG, ToFFG)
        tries <- tries + 1
    }
    if(nrow(sel) > L) sel <- sel[sample.int(nrow(sel), L, replace = FALSE), , drop = FALSE]
    sel
}

calc_metrics_from_edges <- function(nodes_ffg, ed_dir, pred_thr_set, modularity_iter = 1){
    nodes_ffg <- unique(as.character(nodes_ffg))
    
    g_dir <- graph_from_data_frame(
        ed_dir %>% rename(from = FromFFG, to = ToFFG),
        directed = TRUE,
        vertices = data.frame(name = nodes_ffg, stringsAsFactors = FALSE)
    )
    
    S <- vcount(g_dir)
    L <- ecount(g_dir)
    
    LinkDensity <- if(S > 0) L / S else NA_real_
    Cdir <- if(S > 1) L / (S * (S - 1)) else NA_real_
    
    g_und <- as_undirected(g_dir, mode = "collapse")
    
    mod <- NA_real_
    nmods_mean <- NA_real_
    if(S > 1 && ecount(g_und) > 0){
        if(modularity_iter <= 1){
            cl <- igraph::cluster_louvain(g_und)
            mod <- igraph::modularity(cl)
            nmods_mean <- length(unique(igraph::membership(cl)))
        } else {
            ms <- modularity_louvain_mean(g_und, n = modularity_iter, seed_base = 2025)
            mod <- ms$mod_mean
            nmods_mean <- ms$nmod_mean
        }
    }
    
    part_once <- louvain_partition_once(g_und, seed = 2025)
    membership_vec <- part_once$membership
    
    clust <- if(S > 2 && ecount(g_und) > 0) transitivity(g_und, type = "global") else NA_real_
    
    indeg_g <- degree(g_dir, mode = "in")
    outdeg_g <- degree(g_dir, mode = "out")
    
    basal_g <- indeg_g == 0
    Prop_Basal <- if(S > 0) sum(basal_g) / S else NA_real_
    
    deg_total <- indeg_g + outdeg_g
    Degree_skewness <- skewness_fp(deg_total)
    
    TL <- calc_trophic_level(g_dir)
    consumers_idx <- which(indeg_g > 0)
    n_consumers <- length(consumers_idx)
    
    Mean_TL <- if(n_consumers > 0) mean(TL[consumers_idx], na.rm = TRUE) else NA_real_
    TL_max  <- if(n_consumers > 0) max(TL[consumers_idx], na.rm = TRUE) else NA_real_
    Mean_Generality <- if(n_consumers > 0) mean(indeg_g[consumers_idx], na.rm = TRUE) else NA_real_
    
    Omnivory <- NA_real_
    Prop_Omnivores <- NA_real_
    
    if(S > 0 && L > 0 && n_consumers > 0){
        A <- as.matrix(igraph::as_adjacency_matrix(g_dir, sparse = FALSE))
        sd_list <- c()
        n_omni <- 0
        
        for(j in consumers_idx){
            prey_idx <- which(A[, j] > 0)
            if(length(prey_idx) >= 1){
                prey_tl <- TL[prey_idx]
                
                if(length(prey_tl) == 1){
                    sd_list <- c(sd_list, 0)
                } else {
                    sd_list <- c(sd_list, sd(prey_tl, na.rm = TRUE))
                }
                
                if(length(unique(round(prey_tl, 3))) >= 2){
                    n_omni <- n_omni + 1
                }
            }
        }
        
        sd_list <- sd_list[is.finite(sd_list)]
        Omnivory <- if(length(sd_list) > 0) mean(sd_list) else NA_real_
        Prop_Omnivores <- n_omni / S
    }
    
    Trophic_incoherence <- NA_real_
    if(L > 0 && length(TL) > 0){
        dx <- TL[as.character(ed_dir$ToFFG)] - TL[as.character(ed_dir$FromFFG)]
        dx <- dx[is.finite(dx)]
        if(length(dx) > 0) Trophic_incoherence <- sqrt(mean((dx - 1)^2))
    }
    
    diet <- build_diet_matrix(nodes_ffg, ed_dir)
    Nestedness_NODF <- nodf01(diet$M)
    Niche_overlap <- mean_jaccard_diet(diet$M)
    
    Dist_between_nodes <- NA_real_
    if(ADD_DISTANCE_METRICS && S > 1 && ecount(g_und) > 0){
        Dist_between_nodes <- igraph::mean_distance(g_und, directed = FALSE, unconnected = TRUE)
    }
    
    Inter_module_link_ratio <- NA_real_
    Mean_PC_Consumers <- NA_real_
    Mean_PC_TopPredators <- NA_real_
    
    # Cross-module feeding metrics for predators.
    Prop_CrossModulePredators_k2 <- NA_real_
    Prop_CrossModulePredators_k3 <- NA_real_
    Mean_PreyModuleRichness_Predators <- NA_real_
    
    if(S > 1 && ecount(g_und) > 0 && nrow(ed_dir) > 0 && !all(is.na(membership_vec))){
        edge_from_mod <- membership_vec[as.character(ed_dir$FromFFG)]
        edge_to_mod   <- membership_vec[as.character(ed_dir$ToFFG)]
        ok_edge_mod <- is.finite(edge_from_mod) & is.finite(edge_to_mod)
        
        if(any(ok_edge_mod)){
            Inter_module_link_ratio <- mean(edge_from_mod[ok_edge_mod] != edge_to_mod[ok_edge_mod], na.rm = TRUE)
        }
        
        pc_vec <- participation_coefficient_und(g_und, membership_vec)
        
        if(n_consumers > 0){
            consumer_names <- names(indeg_g)[consumers_idx]
            consumer_pc <- pc_vec[consumer_names]
            consumer_pc <- consumer_pc[is.finite(consumer_pc)]
            Mean_PC_Consumers <- if(length(consumer_pc) > 0) mean(consumer_pc, na.rm = TRUE) else NA_real_
            
            top_pred_names <- consumer_names[TL[consumer_names] > PRED_TL_THRESHOLD]
            top_pred_pc <- pc_vec[top_pred_names]
            top_pred_pc <- top_pred_pc[is.finite(top_pred_pc)]
            Mean_PC_TopPredators <- if(length(top_pred_pc) > 0) mean(top_pred_pc, na.rm = TRUE) else NA_real_
        }
        
        predator_names <- names(indeg_g)[consumers_idx][TL[consumers_idx] > PRED_TL_THRESHOLD]
        
        if(length(predator_names) > 0){
            prey_module_counts <- c()
            
            for(pred in predator_names){
                prey_nodes <- ed_dir$FromFFG[ed_dir$ToFFG == pred]
                prey_nodes <- unique(prey_nodes)
                
                prey_mods <- membership_vec[prey_nodes]
                prey_mods <- unique(prey_mods[is.finite(prey_mods)])
                
                prey_module_counts <- c(prey_module_counts, length(prey_mods))
            }
            
            prey_module_counts <- prey_module_counts[is.finite(prey_module_counts)]
            
            if(length(prey_module_counts) > 0){
                Mean_PreyModuleRichness_Predators <- mean(prey_module_counts, na.rm = TRUE)
                Prop_CrossModulePredators_k2 <- mean(prey_module_counts >= 2, na.rm = TRUE)
                Prop_CrossModulePredators_k3 <- mean(prey_module_counts >= 3, na.rm = TRUE)
            }
        }
    }
    
    extra_cols_all <- list()
    for(thr in pred_thr_set){
        tag <- thr_tag(thr)
        if(S > 0 && n_consumers > 0){
            n_pred <- sum(TL[consumers_idx] > thr, na.rm = TRUE)
            n_basal_cons <- sum(TL[consumers_idx] <= thr, na.rm = TRUE)
            extra_cols_all[[paste0("Prop_Predators_TL", tag)]] <- n_pred / S
            extra_cols_all[[paste0("Prop_BasalConsumers_TL", tag)]] <- n_basal_cons / S
        } else {
            extra_cols_all[[paste0("Prop_Predators_TL", tag)]] <- NA_real_
            extra_cols_all[[paste0("Prop_BasalConsumers_TL", tag)]] <- NA_real_
        }
    }
    
    tag0 <- thr_tag(PRED_TL_THRESHOLD)
    Prop_Predators <- extra_cols_all[[paste0("Prop_Predators_TL", tag0)]]
    Prop_BasalConsumers <- extra_cols_all[[paste0("Prop_BasalConsumers_TL", tag0)]]
    
    out <- tibble(
        S_nodes = S,
        L_edges = L,
        LinkDensity = LinkDensity,
        Connectance_dir = Cdir,
        Modularity = mod,
        N_modules_mean = nmods_mean,
        Clustering_und = clust,
        Dist_between_nodes = Dist_between_nodes,
        Degree_skewness = Degree_skewness,
        Mean_TL = Mean_TL,
        TL_max = TL_max,
        Mean_Generality = Mean_Generality,
        Omnivory = Omnivory,
        Prop_Omnivores = Prop_Omnivores,
        Nestedness_NODF = Nestedness_NODF,
        Niche_overlap = Niche_overlap,
        Trophic_incoherence = Trophic_incoherence,
        Prop_Basal = Prop_Basal,
        Prop_Predators = Prop_Predators,
        Prop_BasalConsumers = Prop_BasalConsumers,
        Inter_module_link_ratio = Inter_module_link_ratio,
        Mean_PC_Consumers = Mean_PC_Consumers,
        Mean_PC_TopPredators = Mean_PC_TopPredators,
        Prop_CrossModulePredators_k2 = Prop_CrossModulePredators_k2,
        Prop_CrossModulePredators_k3 = Prop_CrossModulePredators_k3,
        Mean_PreyModuleRichness_Predators = Mean_PreyModuleRichness_Predators
    ) %>% bind_cols(as_tibble(extra_cols_all))
    
    out
}

match_basal_fixed <- function(nodes, basal_fixed){
    if(length(nodes) == 0 || length(basal_fixed) == 0) return(character(0))
    nodes[tolower(nodes) %in% tolower(basal_fixed)]
}

calc_metrics_one <- function(nodes_ffg, edges_template){
    nodes_ffg <- unique(c(nodes_ffg, BASAL_ALWAYS))
    pred_thr_set <- unique(sort(c(PRED_TL_THRESHOLD, PRED_TL_SENS_THRESHOLDS)))
    
    if(length(nodes_ffg) == 0){
        extra_cols_all <- list()
        for(thr in pred_thr_set){
            tag <- thr_tag(thr)
            extra_cols_all[[paste0("Prop_Predators_TL", tag)]] <- NA_real_
            extra_cols_all[[paste0("Prop_BasalConsumers_TL", tag)]] <- NA_real_
        }
        return(
            tibble(
                S_nodes = 0,
                L_edges = 0,
                LinkDensity = NA_real_,
                Connectance_dir = NA_real_,
                Modularity = NA_real_,
                N_modules_mean = NA_real_,
                Clustering_und = NA_real_,
                Dist_between_nodes = NA_real_,
                Degree_skewness = NA_real_,
                Mean_TL = NA_real_,
                TL_max = NA_real_,
                Mean_Generality = NA_real_,
                Omnivory = NA_real_,
                Prop_Omnivores = NA_real_,
                Nestedness_NODF = NA_real_,
                Niche_overlap = NA_real_,
                Trophic_incoherence = NA_real_,
                Prop_Basal = NA_real_,
                Prop_Predators = NA_real_,
                Prop_BasalConsumers = NA_real_,
                Inter_module_link_ratio = NA_real_,
                Mean_PC_Consumers = NA_real_,
                Mean_PC_TopPredators = NA_real_,
                Prop_CrossModulePredators_k2 = NA_real_,
                Prop_CrossModulePredators_k3 = NA_real_,
                Mean_PreyModuleRichness_Predators = NA_real_
            ) %>% bind_cols(as_tibble(extra_cols_all))
        )
    }
    
    ed0 <- edges_template %>%
        filter(FromFFG %in% nodes_ffg, ToFFG %in% nodes_ffg) %>%
        distinct(FromFFG, ToFFG) %>%
        filter(FromFFG != ToFFG)
    
    ed_dir <- ed0
    if(EDGE_DIRECTION == "predator_to_prey"){
        ed_dir <- ed0 %>% transmute(FromFFG = ToFFG, ToFFG = FromFFG)
    }
    
    obs <- calc_metrics_from_edges(nodes_ffg, ed_dir, pred_thr_set, MODULARITY_N)
    
    if(isTRUE(RUN_NULL_MODEL) && obs$S_nodes >= 2 && obs$L_edges >= 1){
        basal_keep <- if(NULL_KEEP_BASAL_MODE == "fixed_list") match_basal_fixed(nodes_ffg, BASAL_FIXED) else character(0)
        L <- obs$L_edges
        null_list <- vector("list", NULL_N)
        for(k in seq_len(NULL_N)){
            set.seed(100000 + k)
            ed_null <- randomize_edges_simple(nodes_ffg, L = L, basal_no_pred = basal_keep, max_tries = NULL_MAX_TRIES)
            null_list[[k]] <- calc_metrics_from_edges(nodes_ffg, ed_null, pred_thr_set, 1)
        }
        null_df <- bind_rows(null_list)
        
        for(mn in names(obs)){
            vnull <- null_df[[mn]]
            mu <- suppressWarnings(mean(vnull, na.rm = TRUE))
            sdv <- suppressWarnings(sd(vnull, na.rm = TRUE))
            obs[[paste0(mn, "_null")]] <- mu
            obs[[paste0(mn, "_null_sd")]] <- sdv
            obs[[paste0(mn, "_devnull")]] <- as.numeric(obs[[mn]] - mu)
            obs[[paste0(mn, "_SES")]] <- ifelse(is.finite(sdv) && sdv > 0, as.numeric((obs[[mn]] - mu) / sdv), NA_real_)
        }
        obs$null_n <- NULL_N
    } else {
        obs$null_n <- ifelse(isTRUE(RUN_NULL_MODEL), 0, NA_integer_)
    }
    
    obs
}

sites <- sort(unique(site_ffg$Site))
metrics_list <- vector("list", length(sites))
pb <- txtProgressBar(min = 0, max = length(sites), style = 3)

for(i in seq_along(sites)){
    sid <- sites[i]
    nodes_ffg <- site_ffg %>% filter(Site == sid) %>% pull(FFG)
    m <- calc_metrics_one(nodes_ffg, ffg_edges)
    metrics_list[[i]] <- tibble(Site = sid) %>% bind_cols(m)
    setTxtProgressBar(pb, i)
}
close(pb)

metrics_df <- bind_rows(metrics_list) %>%
    left_join(site_genus_richness, by = "Site") %>%
    mutate(Genus_richness = ifelse(is.na(Genus_richness), 0L, as.integer(Genus_richness))) %>%
    left_join(env, by = "Site")

# -------------------------
# 8) Prepare Covariates
# -------------------------
message("Preparing covariates ...")
df <- metrics_df[is.finite(metrics_df$FI), , drop = FALSE]

FI_sc <- scale_z(df$FI)
df$FI_z <- FI_sc$z

if("mainstem" %in% names(df)){
    if(is.factor(df$mainstem)){
        df$mainstem <- droplevels(df$mainstem)
    } else if(is.character(df$mainstem)){
        df$mainstem <- factor(trimws(df$mainstem))
    } else {
        mainstem_num <- suppressWarnings(as.numeric(df$mainstem))
        if(sum(is.finite(mainstem_num)) > 0){
            uniq_ms <- sort(unique(mainstem_num[is.finite(mainstem_num)]))
            if(length(uniq_ms) <= 5 && all(abs(uniq_ms - round(uniq_ms)) < 1e-8)){
                df$mainstem <- factor(mainstem_num)
            } else {
                df$mainstem <- mainstem_num
            }
        }
    }
}

RAW_ENV_VARS <- c("t", "do", "ph", "tds", "tsm", "cond", "nh4", "no2", "no3", "si", "po4", "chla", "river_width")
cand_env <- RAW_ENV_VARS[RAW_ENV_VARS %in% names(df)]

cov_z_names <- character(0)
for(v in cand_env){
    x <- df[[v]]
    if(sum(is.finite(x)) > 0){
        if(IMPUTE_MISSING_ENV){
            med <- median(x[is.finite(x)], na.rm = TRUE)
            x[!is.finite(x)] <- med
        }
        sc <- scale_z(x)
        v_z <- paste0(v, "_z")
        df[[v_z]] <- sc$z
        cov_z_names <- c(cov_z_names, v_z)
    }
}

# -------------------------
# 9) Main loop (brms modelling)
# -------------------------
all_hop_tbls <- list()

for(model_type in MODEL_TYPES){
    
    message("\n==================================================")
    message(">>> Running Pipeline for: ", model_type, " MODEL <<<")
    message("==================================================")
    
    OUT_DIR <- paste0("FI_FFG_", model_type)
    MODEL_DIR <- file.path(OUT_DIR, "models")
    DIAG_DIR  <- file.path(OUT_DIR, "diagnostics")
    PPC_DIR   <- file.path(OUT_DIR, "ppc")
    LOO_DIR   <- file.path(OUT_DIR, "loo")
    
    dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
    dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
    dir.create(DIAG_DIR,  showWarnings = FALSE, recursive = TRUE)
    dir.create(PPC_DIR,   showWarnings = FALSE, recursive = TRUE)
    dir.create(LOO_DIR,   showWarnings = FALSE, recursive = TRUE)
    
    safe_write_tsv(metrics_df, file.path(OUT_DIR, "Network_metrics_by_site_FFG.tsv"))
    
    build_formula <- function(has_river, has_mainstem = FALSE, has_river_width_total = FALSE){
        if(model_type == "TOTAL"){
            rhs_terms <- c("FI_z")
            if(has_mainstem) rhs_terms <- c(rhs_terms, "mainstem")
            if(has_river_width_total) rhs_terms <- c(rhs_terms, "river_width_z")
            return(as.formula(paste0("y ~ ", paste(rhs_terms, collapse = " + "))))
        } else {
            rhs <- paste(c("FI_z", cov_z_names), collapse = " + ")
            if(USE_RIVER_RE_DIRECT && has_river){
                return(as.formula(paste0("y ~ ", rhs, " + (1|river)")))
            } else {
                return(as.formula(paste0("y ~ ", rhs)))
            }
        }
    }
    
    fit_one_metric <- function(metric){
        metric_col <- metric
        if(RESPONSE_MODE != "raw"){
            if(RESPONSE_MODE == "devnull") metric_col <- paste0(metric, "_devnull")
            if(RESPONSE_MODE == "ses")     metric_col <- paste0(metric, "_SES")
        }
        
        need <- c("Site","river","FI","FI_z", metric_col)
        if(model_type == "TOTAL") need <- c(need, "mainstem", "river_width_z")
        if(model_type == "DIRECT") need <- c(need, cov_z_names)
        
        cols <- intersect(need, names(df))
        d0 <- df[, cols, drop = FALSE]
        
        if(!(metric_col %in% names(d0))){
            return(list(ok = FALSE, metric = metric, reason = paste0("metric missing: ", metric_col)))
        }
        
        ok_idx <- is.finite(d0[[metric_col]]) & is.finite(d0$FI_z)
        if(model_type == "TOTAL"){
            if("mainstem" %in% names(d0)){
                if(is.factor(d0$mainstem) || is.character(d0$mainstem)){
                    ok_idx <- ok_idx & !is.na(d0$mainstem) & as.character(d0$mainstem) != ""
                } else {
                    ok_idx <- ok_idx & is.finite(suppressWarnings(as.numeric(d0$mainstem)))
                }
            }
            if("river_width_z" %in% names(d0)) ok_idx <- ok_idx & is.finite(d0$river_width_z)
        }
        if(model_type == "DIRECT" && !IMPUTE_MISSING_ENV){
            ok_idx <- ok_idx & apply(d0[, cov_z_names, drop = FALSE], 1, function(r) all(is.finite(r)))
        }
        d0 <- d0[ok_idx, , drop = FALSE]
        
        if(nrow(d0) < 12){
            return(list(ok = FALSE, metric = metric, reason = "too few rows"))
        }
        
        fam_tag <- if(RESPONSE_MODE != "raw") {
            "gaussian"
        } else if(is_prop_metric(metric)) {
            "beta"
        } else if(is_count_metric(metric)) {
            "negbin"
        } else {
            "gaussian"
        }
        
        y_raw <- suppressWarnings(as.numeric(d0[[metric_col]]))
        
        if(fam_tag == "beta"){
            y <- squeeze01(y_raw, nrow(d0))
            family_brms <- Beta(link = "logit")
            family_tmb  <- glmmTMB::beta_family(link = "logit")
        } else if(fam_tag == "negbin"){
            y <- pmax(round(y_raw), 0)
            family_brms <- negbinomial(link = "log")
            family_tmb  <- glmmTMB::nbinom2(link = "log")
        } else {
            y <- y_raw
            family_brms <- gaussian()
            family_tmb  <- gaussian()
        }
        
        d1 <- d0
        d1$y <- y
        
        has_river <- ("river" %in% names(d1)) &&
            any(!is.na(d1$river)) &&
            length(unique(d1$river[!is.na(d1$river)])) >= 2
        
        has_mainstem_total <- identical(model_type, "TOTAL") && ("mainstem" %in% names(d1)) && has_valid_term(d1$mainstem)
        has_river_width_total <- identical(model_type, "TOTAL") && ("river_width_z" %in% names(d1)) && has_valid_term(d1$river_width_z)
        
        fml <- build_formula(
            has_river = has_river,
            has_mainstem = has_mainstem_total,
            has_river_width_total = has_river_width_total
        )
        
        has_river_in_model <- identical(model_type, "DIRECT") && USE_RIVER_RE_DIRECT && has_river
        
        model_path <- file.path(MODEL_DIR, paste0("model_", metric, "_", fam_tag, "_", RESPONSE_MODE, "_", RUN_MODE, ".rds"))
        diag_summary_path <- file.path(DIAG_DIR, paste0("diag_summary_", metric, "_", fam_tag, "_", RESPONSE_MODE, "_", RUN_MODE, ".tsv"))
        loo_prefix <- file.path(LOO_DIR, paste0("loo_", metric, "_", fam_tag, "_", RESPONSE_MODE, "_", RUN_MODE))
        
        if(CAN_TRY_BRMS){
            fit <- NULL
            fit_ok <- TRUE
            
            if(file.exists(model_path)){
                fit <- readRDS(model_path)
            } else {
                pri <- make_priors(fam_tag, has_river_in_model)
                fit <- tryCatch({
                    brm(
                        formula = fml,
                        data = d1,
                        family = family_brms,
                        prior = pri,
                        chains = CHAINS,
                        iter = ITER,
                        warmup = WARMUP,
                        cores = CORES,
                        backend = "rstan",
                        seed = STAN_SEED,
                        control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TDEPTH),
                        refresh = 0
                    )
                }, error = function(e){
                    fit_ok <<- FALSE
                    NULL
                })
                if(fit_ok && !is.null(fit)) saveRDS(fit, model_path)
            }
            
            if(fit_ok && !is.null(fit)){
                if(SAVE_DIAGNOSTICS){
                    dg <- diagnose_brms_fit(fit, metric, fam_tag, MAX_TDEPTH)
                    safe_write_tsv(dg$summary, diag_summary_path)
                }
                
                if(RUN_LOO_COMPARE){
                    loo_out <- save_loo_outputs(fit, metric, loo_prefix)
                    fit <- loo_out$fit
                    saveRDS(fit, model_path)
                }
                
                # Use tryCatch to avoid failure when Stan returns no posterior draws.
                draws <- tryCatch({
                    as_draws_df(fit)
                }, error = function(e) NULL)
                
                if(is.null(draws)) {
                    return(list(ok = FALSE, metric = metric, reason = "Model contains no samples (Stan crash or iter=warmup)"))
                }
                # ----------------------------------------------------
                
                if(!("b_FI_z" %in% names(draws))){
                    return(list(ok = FALSE, metric = metric, reason = "FI coef missing"))
                }
                
                b <- draws[["b_FI_z"]]
                pd_pos <- mean(b > 0)
                hop <- max(pd_pos, 1 - pd_pos)
                ci <- quantile(b, c(0.025, 0.5, 0.975), na.rm = TRUE)
                
                FI_seq_raw <- seq(min(d1$FI, na.rm = TRUE), max(d1$FI, na.rm = TRUE), length.out = FI_GRID_N)
                FI_seq_z   <- (FI_seq_raw - FI_sc$mu) / FI_sc$sd
                
                newdata <- data.frame(FI_z = FI_seq_z, stringsAsFactors = FALSE)
                if(model_type == "TOTAL"){
                    if(has_mainstem_total) newdata$mainstem <- ref_value_for_prediction(d1$mainstem)
                    if(has_river_width_total) newdata$river_width_z <- 0
                }
                if(model_type == "DIRECT"){
                    for(cv in cov_z_names) newdata[[cv]] <- 0
                }
                if(model_type == "DIRECT" && has_river_in_model){
                    newdata$river <- d1$river[which(!is.na(d1$river))[1]]
                }
                
                ep <- posterior_epred(fit, newdata = newdata, re_formula = NA)
                pred_df <- data.frame(
                    FI = FI_seq_raw,
                    mu = as.numeric(apply(ep, 2, mean)),
                    lo = as.numeric(apply(ep, 2, quantile, probs = 0.05)),
                    hi = as.numeric(apply(ep, 2, quantile, probs = 0.95))
                )
                
                keep_id <- sample(seq_len(nrow(ep)), min(NDRAW_SPAGHETTI, nrow(ep)))
                ep_long <- tidyr::pivot_longer(
                    as.data.frame(ep[keep_id, , drop = FALSE]) %>% mutate(draw = row_number()),
                    cols = -draw,
                    names_to = "k",
                    values_to = "mu_draw"
                )
                ep_long$k  <- as.integer(gsub("\\D+", "", ep_long$k))
                ep_long$FI <- FI_seq_raw[ep_long$k]
                
                return(list(
                    ok = TRUE,
                    method = "brms",
                    metric = metric,
                    fam = fam_tag,
                    fit = fit,
                    ci = ci,
                    pd_pos = pd_pos,
                    hop = hop,
                    pred_df = pred_df,
                    ep_long = ep_long
                ))
            }
        }
        
        if(!ALLOW_FALLBACK_GLMMTMB){
            return(list(ok = FALSE, metric = metric, reason = "brms failed and fallback disabled"))
        }
        
        fit2 <- tryCatch({
            glmmTMB(formula = fml, data = d1, family = family_tmb)
        }, error = function(e) NULL)
        
        if(is.null(fit2)){
            return(list(ok = FALSE, metric = metric, reason = "all models failed"))
        }
        
        sm <- summary(fit2)
        coefs <- sm$coefficients$cond
        est <- coefs["FI_z","Estimate"]
        se  <- coefs["FI_z","Std. Error"]
        pd_pos <- pnorm(est / se)
        ci <- c(est - 1.96 * se, est, est + 1.96 * se)
        
        FI_seq_raw <- seq(min(d1$FI, na.rm = TRUE), max(d1$FI, na.rm = TRUE), length.out = FI_GRID_N)
        newdata <- data.frame(FI_z = (FI_seq_raw - FI_sc$mu) / FI_sc$sd, stringsAsFactors = FALSE)
        if(model_type == "TOTAL"){
            if(has_mainstem_total) newdata$mainstem <- ref_value_for_prediction(d1$mainstem)
            if(has_river_width_total) newdata$river_width_z <- 0
        }
        if(model_type == "DIRECT"){
            for(cv in cov_z_names) newdata[[cv]] <- 0
        }
        if(model_type == "DIRECT" && has_river_in_model){
            newdata$river <- d1$river[which(!is.na(d1$river))[1]]
        }
        
        pr_link <- predict(fit2, newdata = newdata, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
        eta_mu_raw <- as.numeric(pr_link$fit)
        eta_se     <- as.numeric(pr_link$se.fit)
        
        pred_df <- data.frame(
            FI = FI_seq_raw,
            mu = inv_link_mu(eta_mu_raw, fam_tag),
            lo = inv_link_mu(eta_mu_raw - 1.645 * eta_se, fam_tag),
            hi = inv_link_mu(eta_mu_raw + 1.645 * eta_se, fam_tag)
        )
        
        return(list(
            ok = TRUE,
            method = "glmmTMB",
            metric = metric,
            fam = fam_tag,
            fit = fit2,
            ci = ci,
            pd_pos = pd_pos,
            hop = max(pd_pos, 1 - pd_pos),
            pred_df = pred_df,
            ep_long = NULL
        ))
    }
    
    METRICS_BASE <- c(
        "Genus_richness",
        "Modularity", "Connectance_dir", "Clustering_und",
        "Prop_Basal", "Prop_Predators", "Prop_BasalConsumers",
        "Prop_Omnivores",
        "Degree_skewness", "Mean_TL", "Trophic_incoherence",
        "Mean_Generality", "Omnivory", "Nestedness_NODF", "Niche_overlap",
        "Inter_module_link_ratio",
        "Mean_PC_Consumers",
        "Mean_PC_TopPredators",
        "Prop_CrossModulePredators_k2",
        "Prop_CrossModulePredators_k3",
        "Mean_PreyModuleRichness_Predators"
    )
    
    if(USE_LINK_DENSITY)     METRICS_BASE <- unique(c(METRICS_BASE, "LinkDensity"))
    if(ADD_DISTANCE_METRICS) METRICS_BASE <- unique(c(METRICS_BASE, "Dist_between_nodes"))
    
    METRICS_TO_MODEL <- unique(METRICS_BASE)
    
    results <- list()
    for(m in METRICS_TO_MODEL){
        message("---- Fitting metric: ", m, " ----")
        results[[m]] <- fit_one_metric(m)
    }
    
    hop_tbl <- bind_rows(lapply(results, function(res){
        if(is.null(res) || !isTRUE(res$ok)){
            return(data.frame(Metric = res$metric, ok = FALSE, note = res$reason))
        }
        data.frame(
            Metric = res$metric,
            ok = TRUE,
            method = res$method,
            family = res$fam,
            FI_slope_median = as.numeric(res$ci[2]),
            FI_slope_CI_lo  = as.numeric(res$ci[1]),
            FI_slope_CI_hi  = as.numeric(res$ci[3]),
            pd_pos = as.numeric(res$pd_pos),
            HOP = as.numeric(res$hop),
            note = ifelse(res$method == "glmmTMB", "HOP is pseudo", "HOP from posterior"),
            stringsAsFactors = FALSE
        )
    }))
    
    hop_tbl$Model_Type <- model_type
    all_hop_tbls[[model_type]] <- hop_tbl
    safe_write_tsv(hop_tbl, file.path(OUT_DIR, paste0("FI_model_HOP_table_", RESPONSE_MODE, ".tsv")))
    
    # ---------- Plotting ----------
    AR <- 1.5
    COLS <- 3
    H_UNIT <- 4
    
    plot_resp <- function(res){
        ggplot(res$pred_df, aes(x = FI)) +
            geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#93C5FD", alpha = 0.28) +
            { if(!is.null(res$ep_long)) geom_line(data = res$ep_long, aes(y = mu_draw, group = draw), color = "#60A5FA", alpha = 0.10, linewidth = 0.35) } +
            geom_line(aes(y = mu), color = "black", linewidth = 1.05) +
            theme_bw(base_size = 10) +
            labs(
                title = paste0(res$metric, " | FI effect"),
                subtitle = paste0(
                    "mode=", RESPONSE_MODE,
                    " | pd(FI>0)=", sprintf("%.3f", res$pd_pos),
                    " | HOP=", sprintf("%.3f", res$hop)
                ),
                x = "FI (raw)",
                y = paste0("Predicted ", res$metric)
            ) +
            theme(
                panel.grid.minor = element_blank(),
                aspect.ratio = 1 / AR
            )
    }
    
    blank_panel <- function(metric, reason = NULL){
        ggplot() +
            theme_void() +
            annotate("text", x = 0, y = 0, label = if(!is.null(reason)) paste0("FAILED: ", reason) else "FAILED", size = 4) +
            labs(title = metric) +
            theme(
                plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
                aspect.ratio = 1 / AR
            )
    }
    
    PLOT_ORDER <- c(
        "Genus_richness",
        "Connectance_dir", "LinkDensity", "Nestedness_NODF",
        "Modularity", "Dist_between_nodes", "Clustering_und", "Niche_overlap",
        "Inter_module_link_ratio", "Mean_PC_Consumers", "Mean_PC_TopPredators",
        "Prop_CrossModulePredators_k2", "Prop_CrossModulePredators_k3", "Mean_PreyModuleRichness_Predators",
        "Prop_Basal", "Prop_Predators", "Prop_BasalConsumers",
        "Degree_skewness", "Mean_TL", "Mean_Generality", "Omnivory", "Prop_Omnivores", "Trophic_incoherence"
    )
    PLOT_METRICS <- intersect(PLOT_ORDER, METRICS_TO_MODEL)
    
    make_panel_list <- function(metric_vec){
        lapply(metric_vec, function(mn){
            if(!is.null(results[[mn]]) && isTRUE(results[[mn]]$ok)){
                plot_resp(results[[mn]])
            } else {
                blank_panel(mn, if(!is.null(results[[mn]])) results[[mn]]$reason else "not fitted")
            }
        })
    }
    
    if(length(results) > 0 && length(PLOT_METRICS) > 0){
        GRID_W <- (H_UNIT * AR) * COLS
        NROW_ALL <- ceiling(length(PLOT_METRICS) / COLS)
        
        out_all <- file.path(OUT_DIR, paste0("FI_partial_effects_BayesCurves_GRID_ALL_", RESPONSE_MODE, ".pdf"))
        tmp_all <- open_pdf_tmp(out_all, width = GRID_W, height = H_UNIT * NROW_ALL)
        print(wrap_plots(make_panel_list(PLOT_METRICS), ncol = COLS))
        close_pdf_tmp(tmp_all, out_all)
    }
}