#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phylolm)
  library(ggplot2)
  library(jsonlite)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))]
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

default_master_csv <- file.path(
  project_root, "output", "runs", "mixed_cellpose_yolo_full_dataset_v1_bgclean",
  "integrated_multinode_phylo", "integrated_master_observed.csv"
)
default_tree_file <- file.path(dirname(project_root), "Desmognathus_TE", "results", "phylogeny", "processed_phylogeny.nwk")
default_output_dir <- file.path(
  project_root, "output", "runs", "mixed_cellpose_yolo_full_dataset_v1_bgclean",
  "integrated_multinode_phylo"
)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

master_csv <- arg_value("--master-csv", default_master_csv)
tree_file <- arg_value("--tree-file", default_tree_file)
output_dir <- arg_value("--output-dir", default_output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

standardize_species <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("_", " ", x)
  x <- sub("^Desmognathus\\s+", "", x, ignore.case = TRUE)
  x <- sub("^D[.]?\\s*", "", x)
  trimws(x)
}

display_species <- function(x) {
  sprintf("D. %s", standardize_species(x))
}

zscore <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) {
    return(rep(0, length(x)))
  }
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

safe_log10 <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log10(x[ok])
  out
}

as_boolish <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  chr <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(chr))
  out[chr %in% c("TRUE", "T")] <- TRUE
  out[chr %in% c("FALSE", "F")] <- FALSE
  if (all(!is.na(out) | is.na(chr) | chr == "")) {
    return(as.logical(out))
  }
  x
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

table_rows <- function(df, cols) {
  out <- character(0)
  for (i in seq_len(nrow(df))) {
    row <- df[i, , drop = FALSE]
    cells <- vapply(cols, function(col) sprintf("<td>%s</td>", row[[col]][[1]]), character(1))
    out <- c(out, sprintf("<tr>%s</tr>", paste(cells, collapse = "")))
  }
  paste(out, collapse = "\n")
}

bind_rows_fill <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) {
    return(data.frame())
  }
  all_cols <- unique(unlist(lapply(dfs, colnames), use.names = FALSE))
  aligned <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_cols, colnames(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })
  do.call(rbind, aligned)
}

master_df <- read.csv(master_csv, stringsAsFactors = FALSE, check.names = FALSE)
for (col in colnames(master_df)) {
  master_df[[col]] <- as_boolish(master_df[[col]])
}
master_df$species_key <- standardize_species(master_df$species)
master_df$species_display <- display_species(master_df$species)

tree <- read.tree(tree_file)
tree$tip.label <- standardize_species(tree$tip.label)
master_df <- master_df[master_df$species_key %in% tree$tip.label, , drop = FALSE]
tree <- keep.tip(tree, intersect(tree$tip.label, master_df$species_key))
tree <- collapse.singles(tree)
master_df <- master_df[match(tree$tip.label, master_df$species_key), , drop = FALSE]
row.names(master_df) <- master_df$species_key
tree$tip.label <- master_df$species_key

if (is.null(tree$edge.length) || anyNA(tree$edge.length) || all(tree$edge.length == 0)) {
  tree <- compute.brlen(tree, method = "Grafen")
}

master_df$species_display <- display_species(master_df$species_key)

continuous_traits <- list(
  genome_size_pg = list(label = "Genome size (pg)", log_transform = TRUE, family = "core"),
  morph_nucleus_area_um2 = list(label = "Nucleus area (um^2)", log_transform = TRUE, family = "core"),
  morph_cell_area_um2 = list(label = "Cell area (um^2)", log_transform = TRUE, family = "core"),
  morph_nc_ratio = list(label = "N:C area ratio", log_transform = TRUE, family = "core"),
  order_pielou = list(label = "TE evenness (Pielou)", log_transform = FALSE, family = "te"),
  ltr_line_logratio = list(label = "LTR:LINE balance", log_transform = FALSE, family = "te"),
  ectopic_mean_ratio = list(label = "Ectopic recombination index", log_transform = TRUE, family = "ectopic"),
  weighted_te_divergence_p90 = list(label = "TE divergence p90", log_transform = TRUE, family = "te"),
  body_size_proxy_mm = list(label = "Body size proxy (mm)", log_transform = TRUE, family = "life_history"),
  aquaticity_index = list(label = "Aquaticity index", log_transform = FALSE, family = "life_history"),
  ltr_history_age_central_mya = list(label = "LTR history age (mya)", log_transform = TRUE, family = "ltr_history")
)

discrete_traits <- list(
  development_mode = list(label = "Development mode"),
  microhabitat_class = list(label = "Microhabitat class")
)

phylo_impute_continuous <- function(values_raw, tree, log_transform = FALSE) {
  values_raw <- as.numeric(values_raw)
  transformed <- if (log_transform) safe_log10(values_raw) else values_raw
  obs <- which(is.finite(transformed))
  miss <- which(!is.finite(transformed))
  n_all <- length(transformed)
  filled_transformed <- transformed
  fill_source <- rep("missing", n_all)
  fill_source[obs] <- "observed"
  support <- rep(NA_real_, n_all)
  cond_sd <- rep(NA_real_, n_all)

  if (length(obs) >= 3 && length(miss) > 0) {
    V <- vcv.phylo(tree)
    Voo <- V[obs, obs, drop = FALSE]
    one <- matrix(1, nrow = length(obs), ncol = 1)
    y_obs <- matrix(transformed[obs], ncol = 1)
    Vinv_one <- solve(Voo, one)
    Vinv_y <- solve(Voo, y_obs)
    mu <- as.numeric(solve(t(one) %*% Vinv_one, t(one) %*% Vinv_y))

    Vmo <- V[miss, obs, drop = FALSE]
    centered_obs <- matrix(transformed[obs] - mu, ncol = 1)
    pred_mean <- as.numeric(mu + Vmo %*% solve(Voo, centered_obs))
    Vmm <- V[miss, miss, drop = FALSE]
    cond_cov <- Vmm - Vmo %*% solve(Voo, V[obs, miss, drop = FALSE])
    pred_sd <- sqrt(pmax(diag(cond_cov), 0))

    filled_transformed[miss] <- pred_mean
    fill_source[miss] <- "phylofill"
    cond_sd[miss] <- pred_sd
    support[miss] <- 1 / (1 + pred_sd)
  }

  filled_raw <- if (log_transform) 10^filled_transformed else filled_transformed

  cv_mae <- NA_real_
  cv_rmse <- NA_real_
  if (length(obs) >= 5) {
    loo_pred <- rep(NA_real_, length(obs))
    for (i in seq_along(obs)) {
      target <- obs[[i]]
      train <- setdiff(obs, target)
      if (length(train) < 3) {
        next
      }
      V <- vcv.phylo(tree)
      Vtt <- V[train, train, drop = FALSE]
      one <- matrix(1, nrow = length(train), ncol = 1)
      y_train <- matrix(transformed[train], ncol = 1)
      Vinv_one <- solve(Vtt, one)
      Vinv_y <- solve(Vtt, y_train)
      mu <- as.numeric(solve(t(one) %*% Vinv_one, t(one) %*% Vinv_y))
      Vmo <- V[target, train, drop = FALSE]
      loo_pred[[i]] <- as.numeric(mu + Vmo %*% solve(Vtt, matrix(transformed[train] - mu, ncol = 1)))
    }
    ok <- is.finite(loo_pred)
    if (any(ok)) {
      truth <- transformed[obs][ok]
      if (log_transform) {
        truth <- 10^truth
        pred_raw <- 10^loo_pred[ok]
      } else {
        pred_raw <- loo_pred[ok]
      }
      cv_mae <- mean(abs(pred_raw - truth))
      cv_rmse <- sqrt(mean((pred_raw - truth)^2))
    }
  }

  anc_fit <- NULL
  root_state <- NA_real_
  if (sum(is.finite(filled_transformed)) >= 3) {
    anc_fit <- ace(filled_transformed, tree, type = "continuous", method = "ML")
    if (length(anc_fit$ace) >= 1) {
      root_state <- anc_fit$ace[[1]]
      if (log_transform) {
        root_state <- 10^root_state
      }
    }
  }

  list(
    filled_raw = filled_raw,
    filled_transformed = filled_transformed,
    fill_source = fill_source,
    support = support,
    cond_sd = cond_sd,
    cv_mae = cv_mae,
    cv_rmse = cv_rmse,
    n_observed = length(obs),
    n_imputed = sum(fill_source == "phylofill"),
    root_state = root_state,
    anc_fit = anc_fit
  )
}

predict_discrete_tip <- function(target_idx, observed_idx, values, dist_mat) {
  d <- dist_mat[target_idx, observed_idx]
  d <- d[is.finite(d)]
  if (!length(d)) {
    return(NULL)
  }
  ordered <- names(sort(d, method = "radix"))
  weights <- 1 / (as.numeric(d[ordered]) + 1e-6)^2
  weights <- weights / sum(weights)
  states <- values[ordered]
  weight_by_state <- tapply(weights, states, sum)
  weight_by_state <- sort(weight_by_state, decreasing = TRUE)
  list(
    predicted = names(weight_by_state)[[1]],
    support = as.numeric(weight_by_state[[1]])
  )
}

phylo_impute_discrete <- function(values_raw, tree) {
  values_raw <- trimws(as.character(values_raw))
  values_raw[values_raw == ""] <- NA_character_
  obs <- which(!is.na(values_raw))
  miss <- which(is.na(values_raw))
  filled <- values_raw
  fill_source <- rep("missing", length(values_raw))
  fill_source[obs] <- "observed"
  support <- rep(NA_real_, length(values_raw))
  dist_mat <- cophenetic.phylo(tree)

  if (length(obs) >= 2 && length(miss) > 0) {
    for (idx in miss) {
      pred <- predict_discrete_tip(tree$tip.label[[idx]], tree$tip.label[obs], setNames(values_raw[obs], tree$tip.label[obs]), dist_mat)
      if (is.null(pred)) {
        next
      }
      filled[[idx]] <- pred$predicted
      fill_source[[idx]] <- "phylofill"
      support[[idx]] <- pred$support
    }
  }

  cv_accuracy <- NA_real_
  if (length(obs) >= 5) {
    cv_pred <- rep(NA_character_, length(obs))
    for (i in seq_along(obs)) {
      target <- obs[[i]]
      train <- setdiff(obs, target)
      pred <- predict_discrete_tip(tree$tip.label[[target]], tree$tip.label[train], setNames(values_raw[train], tree$tip.label[train]), dist_mat)
      if (is.null(pred)) {
        next
      }
      cv_pred[[i]] <- pred$predicted
    }
    ok <- !is.na(cv_pred)
    if (any(ok)) {
      cv_accuracy <- mean(cv_pred[ok] == values_raw[obs][ok])
    }
  }

  anc_fit <- NULL
  root_state <- NA_character_
  root_support <- NA_real_
  if (sum(!is.na(filled)) >= 3 && length(unique(filled[!is.na(filled)])) >= 2) {
    states <- sort(unique(filled[!is.na(filled)]))
    tip_codes <- match(filled, states)
    names(tip_codes) <- tree$tip.label
    anc_fit <- ace(tip_codes, tree, type = "discrete", method = "ML", model = "ER")
    if (!is.null(anc_fit$lik.anc) && nrow(anc_fit$lik.anc) >= 1) {
      root_probs <- anc_fit$lik.anc[1, ]
      root_idx <- which.max(root_probs)
      root_state <- states[[root_idx]]
      root_support <- as.numeric(root_probs[[root_idx]])
    }
  }

  list(
    filled = filled,
    fill_source = fill_source,
    support = support,
    cv_accuracy = cv_accuracy,
    n_observed = length(obs),
    n_imputed = sum(fill_source == "phylofill"),
    root_state = root_state,
    root_support = root_support,
    anc_fit = anc_fit
  )
}

plot_trait_tree <- function(values_raw, tree, title, outfile, log_transform = FALSE) {
  analysis_values <- if (log_transform) safe_log10(values_raw) else as.numeric(values_raw)
  anc_fit <- ace(analysis_values, tree, type = "continuous", method = "ML")
  all_states <- c(analysis_values, anc_fit$ace)
  edge_vals <- rowMeans(cbind(all_states[tree$edge[, 1]], all_states[tree$edge[, 2]]))
  palette <- colorRampPalette(c("#f7f1e5", "#c58c63", "#7a3e1d"))(128)
  bins <- cut(
    edge_vals,
    breaks = seq(min(all_states, na.rm = TRUE), max(all_states, na.rm = TRUE), length.out = length(palette) + 1),
    include.lowest = TRUE,
    labels = FALSE
  )
  edge_cols <- palette[bins]

  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
  par(mar = c(2, 1, 2.8, 1))
  plot(tree, edge.color = edge_cols, show.tip.label = TRUE, cex = 0.55, main = title)
  par(mar = c(4, 2, 2.8, 5))
  y_breaks <- seq(min(all_states, na.rm = TRUE), max(all_states, na.rm = TRUE), length.out = length(palette) + 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(y_breaks))
  rect(0, y_breaks[-length(y_breaks)], 1, y_breaks[-1], col = palette, border = NA)
  axis(4)
  mtext(side = 4, line = 3, if (log_transform) paste0(title, " (log10)") else title)
}

plot_discrete_trait <- function(values, tree, ace_fit, title) {
  values <- as.character(values)
  states <- sort(unique(values[!is.na(values)]))
  palette <- setNames(rainbow(length(states), s = 0.55, v = 0.85), states)
  plot(tree, show.tip.label = TRUE, cex = 0.55, main = title)
  tiplabels(pch = 21, bg = palette[values], col = "#4d3f32", cex = 1.2)
  if (!is.null(ace_fit) && !is.null(ace_fit$lik.anc)) {
    nodelabels(pie = ace_fit$lik.anc, piecol = unname(palette), cex = 0.22)
  }
  legend("topleft", legend = names(palette), pt.bg = unname(palette), pch = 21, cex = 0.55, bty = "n")
}

imputation_summary_rows <- list()
continuous_node_rows <- list()
discrete_node_rows <- list()

for (trait_id in names(continuous_traits)) {
  spec <- continuous_traits[[trait_id]]
  result <- phylo_impute_continuous(master_df[[trait_id]], tree, log_transform = spec$log_transform)
  master_df[[paste0(trait_id, "_phylofill")]] <- result$filled_raw
  master_df[[paste0(trait_id, "_phylofill_source")]] <- result$fill_source
  master_df[[paste0(trait_id, "_phylofill_support")]] <- result$support
  imputation_summary_rows[[trait_id]] <- data.frame(
    trait = trait_id,
    label = spec$label,
    type = "continuous",
    observed_n = result$n_observed,
    imputed_n = result$n_imputed,
    total_n = nrow(master_df),
    cv_metric_1 = "mae",
    cv_value_1 = result$cv_mae,
    cv_metric_2 = "rmse",
    cv_value_2 = result$cv_rmse,
    root_state = result$root_state,
    stringsAsFactors = FALSE
  )
  if (!is.null(result$anc_fit) && length(result$anc_fit$ace)) {
    continuous_node_rows[[trait_id]] <- data.frame(
      trait = trait_id,
      node_id = names(result$anc_fit$ace),
      reconstructed_state = as.numeric(result$anc_fit$ace),
      stringsAsFactors = FALSE
    )
  }
}

for (trait_id in names(discrete_traits)) {
  spec <- discrete_traits[[trait_id]]
  result <- phylo_impute_discrete(master_df[[trait_id]], tree)
  master_df[[paste0(trait_id, "_phylofill")]] <- result$filled
  master_df[[paste0(trait_id, "_phylofill_source")]] <- result$fill_source
  master_df[[paste0(trait_id, "_phylofill_support")]] <- result$support
  imputation_summary_rows[[paste0(trait_id, "_disc")]] <- data.frame(
    trait = trait_id,
    label = spec$label,
    type = "discrete",
    observed_n = result$n_observed,
    imputed_n = result$n_imputed,
    total_n = nrow(master_df),
    cv_metric_1 = "accuracy",
    cv_value_1 = result$cv_accuracy,
    cv_metric_2 = "root_support",
    cv_value_2 = result$root_support,
    root_state = result$root_state,
    stringsAsFactors = FALSE
  )
  if (!is.null(result$anc_fit) && !is.null(result$anc_fit$lik.anc)) {
    probs <- as.data.frame(result$anc_fit$lik.anc, stringsAsFactors = FALSE)
    probs$node_id <- row.names(probs)
    probs$trait <- trait_id
    discrete_node_rows[[trait_id]] <- probs
  }
}

imputation_summary <- bind_rows_fill(imputation_summary_rows)
write.csv(imputation_summary, file.path(output_dir, "trait_imputation_summary.csv"), row.names = FALSE)
write.csv(master_df, file.path(output_dir, "integrated_master_phylofill.csv"), row.names = FALSE)

if (length(continuous_node_rows)) {
  write.csv(bind_rows_fill(continuous_node_rows), file.path(output_dir, "continuous_ancestral_nodes.csv"), row.names = FALSE)
}
if (length(discrete_node_rows)) {
  write.csv(bind_rows_fill(discrete_node_rows), file.path(output_dir, "discrete_ancestral_nodes.csv"), row.names = FALSE)
}

family_configs <- list(
  genome_morphology = list(
    label = "Genome-Morphology",
    models = list(
      G_to_N_to_C = list(ns = c("gs"), cs = c("ns")),
      G_to_C_to_N = list(cs = c("gs"), ns = c("cs")),
      G_to_both = list(ns = c("gs"), cs = c("gs")),
      G_to_both_plus_N_to_C = list(ns = c("gs"), cs = c("gs", "ns")),
      G_to_both_plus_C_to_N = list(cs = c("gs"), ns = c("gs", "cs"))
    )
  ),
  te_genome = list(
    label = "TE-Genome",
    models = list(
      ltr_balance_only = list(gs = c("ltr_balance")),
      evenness_only = list(gs = c("te_evenness")),
      additive_load_evenness = list(gs = c("ltr_balance", "te_evenness")),
      mediated_evenness = list(te_evenness = c("ltr_balance"), gs = c("te_evenness"))
    )
  ),
  te_genome_organismal = list(
    label = "TE-Genome-Organismal",
    models = list(
      te_baseline = list(te_evenness = c("ltr_balance"), gs = c("te_evenness")),
      body_size_additive = list(te_evenness = c("ltr_balance"), gs = c("te_evenness", "body_size")),
      aquaticity_additive = list(te_evenness = c("ltr_balance"), gs = c("te_evenness", "aquaticity")),
      organismal_additive = list(te_evenness = c("ltr_balance"), gs = c("te_evenness", "body_size", "aquaticity")),
      aquaticity_confounds_body_and_te = list(body_size = c("aquaticity"), te_evenness = c("ltr_balance", "aquaticity"), gs = c("te_evenness", "body_size"))
    )
  ),
  te_genome_ltr_history = list(
    label = "TE-Genome-LTR History",
    models = list(
      te_baseline = list(te_evenness = c("ltr_balance"), gs = c("te_evenness")),
      history_direct = list(gs = c("ltr_history")),
      history_additive = list(te_evenness = c("ltr_balance"), gs = c("te_evenness", "ltr_history")),
      history_to_balance = list(ltr_balance = c("ltr_history"), te_evenness = c("ltr_balance"), gs = c("te_evenness")),
      history_to_evenness = list(te_evenness = c("ltr_balance", "ltr_history"), gs = c("te_evenness"))
    )
  ),
  te_genome_ectopic = list(
    label = "TE-Genome-Ectopic",
    models = list(
      ectopic_only = list(gs = c("ectopic_index")),
      ltr_to_ectopic = list(ectopic_index = c("ltr_balance"), gs = c("ectopic_index")),
      evenness_and_ectopic = list(ectopic_index = c("ltr_balance"), gs = c("te_evenness", "ectopic_index")),
      full_mechanism = list(ectopic_index = c("ltr_balance", "te_evenness"), gs = c("ltr_balance", "te_evenness", "ectopic_index"))
    )
  ),
  te_genome_ectopic_organismal = list(
    label = "TE-Genome-Ectopic-Organismal",
    models = list(
      ectopic_baseline = list(gs = c("ectopic_index")),
      ectopic_body_size_additive = list(gs = c("ectopic_index", "body_size")),
      te_body_size_baseline = list(te_evenness = c("ltr_balance"), gs = c("te_evenness", "body_size")),
      te_ectopic_body_size = list(te_evenness = c("ltr_balance"), gs = c("te_evenness", "ectopic_index", "body_size")),
      ltr_to_ectopic_body_size = list(ectopic_index = c("ltr_balance"), gs = c("ectopic_index", "body_size"))
    )
  ),
  te_genome_morphology = list(
    label = "TE-Genome-Morphology",
    models = list(
      te_to_genome_to_nucleus_to_cell = list(gs = c("ltr_balance", "te_evenness"), ns = c("gs"), cs = c("ns")),
      te_to_genome_partial_cell = list(gs = c("ltr_balance", "te_evenness"), ns = c("gs"), cs = c("gs", "ns")),
      te_direct_to_genome_nucleus = list(gs = c("ltr_balance"), ns = c("gs"), cs = c("ns")),
      te_evenness_path = list(te_evenness = c("ltr_balance"), gs = c("te_evenness"), ns = c("gs"), cs = c("ns"))
    )
  )
)

col_for <- function(df, base_col, mode) {
  phy_col <- paste0(base_col, "_phylofill")
  if (mode == "phylofill" && phy_col %in% colnames(df)) {
    return(phy_col)
  }
  base_col
}

prepare_family_data <- function(df, family, mode) {
  get_col <- function(x) col_for(df, x, mode)
  analysis_df <- switch(
    family,
    genome_morphology = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ns = zscore(safe_log10(df[[get_col("morph_nucleus_area_um2")]])),
      cs = zscore(safe_log10(df[[get_col("morph_cell_area_um2")]])),
      stringsAsFactors = FALSE
    ),
    te_genome = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ltr_balance = zscore(df[[get_col("ltr_line_logratio")]]),
      te_evenness = zscore(df[[get_col("order_pielou")]]),
      stringsAsFactors = FALSE
    ),
    te_genome_organismal = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ltr_balance = zscore(df[[get_col("ltr_line_logratio")]]),
      te_evenness = zscore(df[[get_col("order_pielou")]]),
      body_size = zscore(safe_log10(df[[get_col("body_size_proxy_mm")]])),
      aquaticity = zscore(as.numeric(df[[get_col("aquaticity_index")]])),
      stringsAsFactors = FALSE
    ),
    te_genome_ltr_history = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ltr_balance = zscore(df[[get_col("ltr_line_logratio")]]),
      te_evenness = zscore(df[[get_col("order_pielou")]]),
      ltr_history = zscore(safe_log10(df[[get_col("ltr_history_age_central_mya")]])),
      stringsAsFactors = FALSE
    ),
    te_genome_ectopic = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ltr_balance = zscore(df[[get_col("ltr_line_logratio")]]),
      te_evenness = zscore(df[[get_col("order_pielou")]]),
      ectopic_index = zscore(safe_log10(df[[get_col("ectopic_mean_ratio")]])),
      stringsAsFactors = FALSE
    ),
    te_genome_ectopic_organismal = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ltr_balance = zscore(df[[get_col("ltr_line_logratio")]]),
      te_evenness = zscore(df[[get_col("order_pielou")]]),
      ectopic_index = zscore(safe_log10(df[[get_col("ectopic_mean_ratio")]])),
      body_size = zscore(safe_log10(df[[get_col("body_size_proxy_mm")]])),
      stringsAsFactors = FALSE
    ),
    te_genome_morphology = data.frame(
      species = df$species_key,
      gs = zscore(safe_log10(df[[get_col("genome_size_pg")]])),
      ltr_balance = zscore(df[[get_col("ltr_line_logratio")]]),
      te_evenness = zscore(df[[get_col("order_pielou")]]),
      ns = zscore(safe_log10(df[[get_col("morph_nucleus_area_um2")]])),
      cs = zscore(safe_log10(df[[get_col("morph_cell_area_um2")]])),
      stringsAsFactors = FALSE
    ),
    stop("Unknown family: ", family)
  )
  analysis_df <- analysis_df[stats::complete.cases(analysis_df), , drop = FALSE]
  analysis_df <- analysis_df[!duplicated(analysis_df$species), , drop = FALSE]
  analysis_df
}

fit_equation <- function(response, predictors, analysis_df, analysis_tree) {
  rhs <- if (length(predictors)) paste(predictors, collapse = " + ") else "1"
  form <- as.formula(sprintf("%s ~ %s", response, rhs))
  fit <- phylolm(form, data = analysis_df, phy = analysis_tree, model = "lambda")
  sm <- summary(fit)
  coef_df <- as.data.frame(sm$coefficients)
  coef_df$term <- row.names(coef_df)
  coef_df$response <- response
  coef_df$lambda <- unname(fit$optpar)
  coef_df$logLik <- as.numeric(logLik(fit))
  coef_df$k <- attr(logLik(fit), "df")
  row.names(coef_df) <- NULL
  list(fit = fit, coefs = coef_df)
}

calc_aicc <- function(logLik_value, k, n) {
  aic <- -2 * logLik_value + 2 * k
  if ((n - k - 1) <= 0) {
    return(aic)
  }
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

run_family_models <- function(family, mode, analysis_df, tree, output_dir) {
  model_specs <- family_configs[[family]]$models
  analysis_tree <- keep.tip(tree, analysis_df$species)
  analysis_df <- analysis_df[match(analysis_tree$tip.label, analysis_df$species), , drop = FALSE]
  row.names(analysis_df) <- analysis_df$species

  model_rows <- list()
  coef_rows <- list()
  for (model_name in names(model_specs)) {
    spec <- model_specs[[model_name]]
    equation_fits <- lapply(names(spec), function(response) fit_equation(response, spec[[response]], analysis_df, analysis_tree))
    names(equation_fits) <- names(spec)
    total_logLik <- sum(vapply(equation_fits, function(x) x$coefs$logLik[1], numeric(1)))
    total_k <- sum(vapply(equation_fits, function(x) x$coefs$k[1], numeric(1)))
    n_species <- nrow(analysis_df)
    aicc <- calc_aicc(total_logLik, total_k, n_species)
    model_rows[[model_name]] <- data.frame(
      family = family,
      mode = mode,
      model = model_name,
      equations = paste(vapply(names(spec), function(response) {
        rhs <- if (length(spec[[response]])) paste(spec[[response]], collapse = " + ") else "1"
        sprintf("%s ~ %s", response, rhs)
      }, character(1)), collapse = " ; "),
      n_species = n_species,
      total_logLik = total_logLik,
      total_k = total_k,
      AICc = aicc,
      stringsAsFactors = FALSE
    )
    coef_df <- do.call(rbind, lapply(equation_fits, function(x) x$coefs))
    coef_df$model <- model_name
    coef_df$family <- family
    coef_df$mode <- mode
    coef_rows[[model_name]] <- coef_df
  }

  model_scores <- do.call(rbind, model_rows)
  model_scores <- model_scores[order(model_scores$AICc), , drop = FALSE]
  model_scores$delta_AICc <- model_scores$AICc - min(model_scores$AICc)
  model_scores$model_weight <- exp(-0.5 * model_scores$delta_AICc)
  model_scores$model_weight <- model_scores$model_weight / sum(model_scores$model_weight)

  coef_df <- do.call(rbind, coef_rows)
  coef_df <- coef_df[, c("family", "mode", "model", "response", "term", "Estimate", "StdErr", "p.value", "lambda")]
  names(coef_df) <- c("family", "mode", "model", "response", "term", "estimate", "std_error", "p_value", "lambda")

  write.csv(model_scores, file.path(output_dir, sprintf("%s_%s_model_scores.csv", family, mode)), row.names = FALSE)
  write.csv(coef_df, file.path(output_dir, sprintf("%s_%s_model_coefficients.csv", family, mode)), row.names = FALSE)

  list(model_scores = model_scores, coef_df = coef_df, analysis_df = analysis_df)
}

panel_summary_rows <- list()
all_model_scores <- list()
all_best_coeffs <- list()

for (family in names(family_configs)) {
  for (mode in c("observed", "phylofill")) {
    analysis_df <- prepare_family_data(master_df, family, mode)
    panel_summary_rows[[paste(family, mode, sep = "_")]] <- data.frame(
      family = family,
      family_label = family_configs[[family]]$label,
      mode = mode,
      n_species = nrow(analysis_df),
      species_list = paste(display_species(analysis_df$species), collapse = "; "),
      stringsAsFactors = FALSE
    )
    if (nrow(analysis_df) < 6) {
      next
    }
    fit <- run_family_models(family, mode, analysis_df, tree, output_dir)
    all_model_scores[[paste(family, mode, sep = "_")]] <- fit$model_scores
    best_model <- fit$model_scores$model[[1]]
    best_coef_df <- subset(fit$coef_df, model == best_model & term != "(Intercept)")
    if (nrow(best_coef_df)) {
      best_coef_df$family_label <- family_configs[[family]]$label
      all_best_coeffs[[paste(family, mode, sep = "_")]] <- best_coef_df
    }
  }
}

panel_summary <- bind_rows_fill(panel_summary_rows)
write.csv(panel_summary, file.path(output_dir, "panel_summary_observed_vs_phylofill.csv"), row.names = FALSE)

model_summary <- if (length(all_model_scores)) bind_rows_fill(all_model_scores) else data.frame()
best_model_summary <- data.frame()
if (nrow(model_summary)) {
  best_model_summary <- bind_rows_fill(
    lapply(split(model_summary, interaction(model_summary$family, model_summary$mode, drop = TRUE)), function(df) df[1, , drop = FALSE])
  )
  best_model_summary$family_label <- unname(vapply(best_model_summary$family, function(x) family_configs[[x]]$label, character(1)))
  write.csv(model_summary, file.path(output_dir, "path_model_summary_all.csv"), row.names = FALSE)
  write.csv(best_model_summary, file.path(output_dir, "path_model_best_models.csv"), row.names = FALSE)
}

best_coeff_summary <- if (length(all_best_coeffs)) bind_rows_fill(all_best_coeffs) else data.frame()
if (nrow(best_coeff_summary)) {
  write.csv(best_coeff_summary, file.path(output_dir, "path_model_best_coefficients.csv"), row.names = FALSE)
}

coverage_plot <- ggplot(panel_summary, aes(x = reorder(family_label, n_species), y = n_species, fill = mode)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  coord_flip() +
  scale_fill_manual(values = c(observed = "#c98d63", phylofill = "#7a3e1d")) +
  labs(title = "Panel coverage before and after phylogenetic filling", x = NULL, y = "Species in panel", fill = NULL) +
  theme_classic(base_size = 12)
ggsave(file.path(output_dir, "figure_1_panel_coverage.png"), coverage_plot, width = 10.5, height = 5.8, dpi = 180)

heatmap_traits <- c(
  "genome_size_pg", "morph_nucleus_area_um2", "morph_cell_area_um2",
  "order_pielou", "ltr_line_logratio", "ectopic_mean_ratio",
  "body_size_proxy_mm", "aquaticity_index", "ltr_history_age_central_mya",
  "development_mode", "microhabitat_class"
)
heatmap_labels <- c(
  genome_size_pg = "Genome",
  morph_nucleus_area_um2 = "Nucleus",
  morph_cell_area_um2 = "Cell",
  order_pielou = "TE evenness",
  ltr_line_logratio = "LTR:LINE",
  ectopic_mean_ratio = "Ectopic",
  body_size_proxy_mm = "Body size",
  aquaticity_index = "Aquaticity",
  ltr_history_age_central_mya = "LTR history",
  development_mode = "Development",
  microhabitat_class = "Microhabitat"
)

heatmap_rows <- list()
for (trait_id in heatmap_traits) {
  source_col <- paste0(trait_id, "_phylofill_source")
  status <- if (source_col %in% colnames(master_df)) master_df[[source_col]] else ifelse(!is.na(master_df[[trait_id]]), "observed", "missing")
  heatmap_rows[[trait_id]] <- data.frame(
    species = factor(master_df$species_display, levels = rev(master_df$species_display)),
    trait = heatmap_labels[[trait_id]],
    status = factor(status, levels = c("observed", "phylofill", "missing")),
    stringsAsFactors = FALSE
  )
}
heatmap_df <- bind_rows_fill(heatmap_rows)

heatmap_plot <- ggplot(heatmap_df, aes(x = trait, y = species, fill = status)) +
  geom_tile(color = "#f7f3ee", linewidth = 0.2) +
  scale_fill_manual(values = c(observed = "#c98d63", phylofill = "#7a3e1d", missing = "#ddd5cb")) +
  labs(title = "Trait coverage across species", x = NULL, y = NULL, fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )
ggsave(file.path(output_dir, "figure_2_trait_status_heatmap.png"), heatmap_plot, width = 10.5, height = 9.5, dpi = 180)

continuous_plot_traits <- c("genome_size_pg", "morph_nucleus_area_um2", "morph_cell_area_um2", "order_pielou", "ectopic_mean_ratio", "body_size_proxy_mm")
png(file.path(output_dir, "figure_3_continuous_trait_ancestry.png"), width = 2200, height = 1800, res = 180)
par(mfrow = c(2, 3))
for (trait_id in continuous_plot_traits) {
  filled_col <- paste0(trait_id, "_phylofill")
  spec <- continuous_traits[[trait_id]]
  plot_trait_tree(master_df[[filled_col]], tree, spec$label, "", log_transform = spec$log_transform)
}
dev.off()

disc_dev <- phylo_impute_discrete(master_df$development_mode_phylofill, tree)
disc_micro <- phylo_impute_discrete(master_df$microhabitat_class_phylofill, tree)

png(file.path(output_dir, "figure_4_discrete_trait_ancestry.png"), width = 1800, height = 950, res = 180)
par(mfrow = c(1, 2), mar = c(1.5, 1.5, 2.5, 1.5))
plot_discrete_trait(master_df$development_mode_phylofill, tree, disc_dev$anc_fit, "Development mode ancestry")
plot_discrete_trait(master_df$microhabitat_class_phylofill, tree, disc_micro$anc_fit, "Microhabitat ancestry")
dev.off()

if (nrow(best_model_summary)) {
  best_model_summary$family_mode <- factor(
    paste(best_model_summary$family_label, ifelse(best_model_summary$mode == "observed", "(observed)", "(phylofilled)")),
    levels = rev(paste(best_model_summary$family_label, ifelse(best_model_summary$mode == "observed", "(observed)", "(phylofilled)")))
  )
  support_plot <- ggplot(best_model_summary, aes(x = model_weight, y = family_mode, color = mode)) +
    geom_vline(xintercept = c(0.25, 0.5, 0.75), color = "#ece5db", linewidth = 0.4) +
    geom_point(size = 3.4) +
    geom_text(aes(label = sprintf("%s (n=%d)", model, n_species)), hjust = -0.03, size = 3.1, family = "serif") +
    scale_color_manual(values = c(observed = "#c98d63", phylofill = "#7a3e1d")) +
    scale_x_continuous(limits = c(0, 1.25)) +
    labs(title = "Best model support by family", x = "AICc model weight", y = NULL, color = NULL) +
    theme_classic(base_size = 12)
  ggsave(file.path(output_dir, "figure_5_best_model_support.png"), support_plot, width = 11.5, height = 6.8, dpi = 180)
}

if (nrow(best_coeff_summary)) {
  best_coeff_summary$edge_label <- paste(best_coeff_summary$response, "<-", best_coeff_summary$term)
  best_coeff_summary$ci_low <- best_coeff_summary$estimate - 1.96 * best_coeff_summary$std_error
  best_coeff_summary$ci_high <- best_coeff_summary$estimate + 1.96 * best_coeff_summary$std_error
  coef_plot <- ggplot(best_coeff_summary, aes(x = estimate, y = edge_label, xmin = ci_low, xmax = ci_high, color = mode)) +
    geom_vline(xintercept = 0, color = "#ddd5cb", linewidth = 0.5) +
    geom_errorbar(width = 0.16, orientation = "y", linewidth = 0.9) +
    geom_point(size = 2.3) +
    facet_wrap(~ family_label, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c(observed = "#c98d63", phylofill = "#7a3e1d")) +
    labs(title = "Best-model coefficients across family panels", x = "Standardized coefficient (95% CI)", y = NULL, color = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file.path(output_dir, "figure_6_best_model_coefficients.png"), coef_plot, width = 11.8, height = 9.0, dpi = 180)
}

summary_json <- list(
  n_species_tree = nrow(master_df),
  n_species_with_current_cellprofiler_updates = 0,
  n_families = length(family_configs),
  n_best_models = nrow(best_model_summary),
  observed_panel_species_max = max(panel_summary$n_species[panel_summary$mode == "observed"]),
  phylofill_panel_species_max = max(panel_summary$n_species[panel_summary$mode == "phylofill"]),
  missing_after_phylofill_continuous = sum(vapply(names(continuous_traits), function(tr) sum(is.na(master_df[[paste0(tr, "_phylofill")]])), numeric(1))),
  missing_after_phylofill_discrete = sum(vapply(names(discrete_traits), function(tr) sum(is.na(master_df[[paste0(tr, "_phylofill")]])), numeric(1)))
)
summary_json$n_species_with_current_cellprofiler_updates <- sum(master_df$current_cellprofiler_bridge, na.rm = TRUE)
write(toJSON(summary_json, pretty = TRUE, auto_unbox = TRUE), file.path(output_dir, "summary.json"))

imputation_html <- imputation_summary
imputation_html$observed_n <- as.integer(imputation_html$observed_n)
imputation_html$imputed_n <- as.integer(imputation_html$imputed_n)
imputation_html$total_n <- as.integer(imputation_html$total_n)
imputation_html$cv_value_1 <- fmt_num(imputation_html$cv_value_1)
imputation_html$cv_value_2 <- fmt_num(imputation_html$cv_value_2)
imputation_html$root_state <- ifelse(is.na(imputation_html$root_state), "NA", as.character(imputation_html$root_state))

panel_html <- panel_summary
panel_html$n_species <- as.integer(panel_html$n_species)

best_model_html <- best_model_summary
if (nrow(best_model_html)) {
  best_model_html$AICc <- fmt_num(best_model_html$AICc)
  best_model_html$delta_AICc <- fmt_num(best_model_html$delta_AICc)
  best_model_html$model_weight <- fmt_num(best_model_html$model_weight)
}

html <- sprintf(
'<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Integrated Desmognathus Phylogenetic Path Report</title>
  <style>
    :root { --ink:#241d17; --muted:#66584a; --line:#d9ccbf; --accent:#7a3e1d; --paper:#fffdfa; }
    html, body { margin: 0; padding: 0; background: #f7f3ee; color: var(--ink); }
    body { font: 16px/1.65 Georgia, serif; }
    main { max-width: 1060px; margin: 0 auto; padding: 40px 28px 56px; background: var(--paper); box-shadow: 0 8px 30px rgba(33, 23, 14, 0.06); }
    h1, h2 { font-weight: 600; }
    h1 { margin: 0 0 10px; font-size: 2.2rem; line-height: 1.2; text-align: center; }
    h2 { margin: 34px 0 14px; font-size: 1.2rem; border-bottom: 1px solid var(--line); padding-bottom: 6px; }
    p { margin: 0 0 14px; }
    .eyebrow { text-align: center; text-transform: uppercase; letter-spacing: 0.14em; color: var(--muted); font-size: 0.75rem; margin-bottom: 10px; }
    .subhead { text-align: center; color: var(--muted); margin-bottom: 24px; }
    .summary { border-top: 2px solid #cdb8a3; border-bottom: 2px solid #cdb8a3; padding: 18px 0 6px; margin-bottom: 26px; }
    .links { margin: 18px 0 8px; }
    .links a { display: inline-block; margin: 0 10px 10px 0; color: var(--accent); text-decoration: none; border-bottom: 1px solid #ccb39e; }
    .figure { margin: 30px 0 36px; }
    .figure img { width: 100%%; height: auto; display: block; border: 1px solid var(--line); background: white; }
    .figcaption { margin-top: 10px; font-size: 0.93rem; color: #43372d; }
    .table-wrap { overflow-x: auto; margin: 12px 0 24px; }
    table { width: 100%%; border-collapse: collapse; }
    th, td { padding: 8px 10px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; font-size: 0.9rem; }
    th { font-size: 0.82rem; letter-spacing: 0.03em; text-transform: uppercase; color: var(--muted); }
    .note { color: var(--muted); font-size: 0.92rem; }
    .mono { font-family: "SFMono-Regular", Menlo, Consolas, monospace; font-size: 0.9em; }
  </style>
</head>
<body>
  <main>
    <header>
      <div class="eyebrow">Integrated Comparative Reconstruction</div>
      <h1>Full-node phylogenetic integration across genome, morphology, TE, ectopic recombination, LTR history, and organismal traits</h1>
      <p class="subhead">Current bg-clean CellProfiler summaries merged with the broader Desmognathus trait workspace and extended with phylogenetic filling</p>
    </header>

    <section class="summary">
      <p>This report integrates %d tree species across %d path-analysis families. Current bg-clean CellProfiler values replaced the legacy genome and morphology estimates for %d species, and missing values across the broader trait matrix were then filled phylogenetically so observed-only and phylofilled analyses could be compared side by side.</p>
      <p class="note">Continuous traits were filled with Brownian-motion conditional expectations on the phylogeny and discrete traits were filled with phylogenetic distance-weighted state voting, while ancestral internal-node states were reconstructed separately for visualization.</p>
      <div class="links">
        <a href="summary.json">summary json</a>
        <a href="integrated_master_observed.csv">observed master</a>
        <a href="integrated_master_phylofill.csv">phylofilled master</a>
        <a href="trait_imputation_summary.csv">trait imputation summary</a>
        <a href="panel_summary_observed_vs_phylofill.csv">panel summary</a>
        <a href="path_model_best_models.csv">best model summary</a>
        <a href="path_model_best_coefficients.csv">best coefficients</a>
      </div>
    </section>

    <section class="figure">
      <img src="figure_1_panel_coverage.png" alt="Panel coverage before and after phylogenetic filling">
      <div class="figcaption"><strong>Figure 1.</strong> Species overlap in each major analysis family before and after phylogenetic filling. The observed-only mode uses only directly measured values, whereas the phylofilled mode retains observed values and supplements missing traits with phylogenetically inferred tip estimates.</div>
    </section>

    <section class="figure">
      <img src="figure_2_trait_status_heatmap.png" alt="Trait coverage heatmap">
      <div class="figcaption"><strong>Figure 2.</strong> Trait availability matrix across the full salamander backbone. Each tile indicates whether a given species-trait combination was directly observed, phylogenetically filled, or still unavailable.</div>
    </section>

    <section class="figure">
      <img src="figure_3_continuous_trait_ancestry.png" alt="Continuous trait ancestry maps">
      <div class="figcaption"><strong>Figure 3.</strong> Continuous-trait ancestral reconstructions across the phylogeny for representative morphology, genome, TE, ectopic, and body-size traits after tip filling.</div>
    </section>

    <section class="figure">
      <img src="figure_4_discrete_trait_ancestry.png" alt="Discrete trait ancestral reconstructions">
      <div class="figcaption"><strong>Figure 4.</strong> Ancestral-state reconstructions for discrete life-history and habitat traits. Tip colors show the final tip-state matrix used for filled analyses and node pies summarize marginal ancestral-state support.</div>
    </section>

    <section class="figure">
      <img src="figure_5_best_model_support.png" alt="Best model support by family">
      <div class="figcaption"><strong>Figure 5.</strong> Best-supported path model in each analysis family under observed-only and phylofilled datasets. Label text gives the winning model and the number of species retained in that family-mode combination.</div>
    </section>

    <section class="figure">
      <img src="figure_6_best_model_coefficients.png" alt="Best model coefficients across families">
      <div class="figcaption"><strong>Figure 6.</strong> Standardized coefficients from the best-supported model in each family. Comparing observed-only and phylofilled estimates helps show whether the larger filled panels change the inferred direction or magnitude of the main paths.</div>
    </section>

    <section>
      <h2>Table 1. Trait Imputation Summary</h2>
      <div class="table-wrap">
        <table>
          <tr><th>Trait</th><th>Type</th><th>Observed</th><th>Imputed</th><th>Total</th><th>CV metric 1</th><th>Value</th><th>CV metric 2</th><th>Value</th><th>Root state</th></tr>
          %s
        </table>
      </div>
    </section>

    <section>
      <h2>Table 2. Panel Sizes</h2>
      <div class="table-wrap">
        <table>
          <tr><th>Family</th><th>Mode</th><th>Species</th><th>Species list</th></tr>
          %s
        </table>
      </div>
    </section>

    <section>
      <h2>Table 3. Best Path Model Per Family</h2>
      <div class="table-wrap">
        <table>
          <tr><th>Family</th><th>Mode</th><th>Best model</th><th>AICc</th><th>Delta AICc</th><th>Weight</th><th>Species</th><th>Equations</th></tr>
          %s
        </table>
      </div>
    </section>
  </main>
</body>
</html>',
  nrow(master_df),
  length(family_configs),
  sum(master_df$current_cellprofiler_bridge, na.rm = TRUE),
  table_rows(imputation_html[, c("label", "type", "observed_n", "imputed_n", "total_n", "cv_metric_1", "cv_value_1", "cv_metric_2", "cv_value_2", "root_state"), drop = FALSE], c("label", "type", "observed_n", "imputed_n", "total_n", "cv_metric_1", "cv_value_1", "cv_metric_2", "cv_value_2", "root_state")),
  table_rows(panel_html[, c("family_label", "mode", "n_species", "species_list"), drop = FALSE], c("family_label", "mode", "n_species", "species_list")),
  if (nrow(best_model_html)) table_rows(best_model_html[, c("family_label", "mode", "model", "AICc", "delta_AICc", "model_weight", "n_species", "equations"), drop = FALSE], c("family_label", "mode", "model", "AICc", "delta_AICc", "model_weight", "n_species", "equations")) else ""
)

writeLines(html, file.path(output_dir, "index.html"))
message("Wrote integrated multi-node report: ", file.path(output_dir, "index.html"))
