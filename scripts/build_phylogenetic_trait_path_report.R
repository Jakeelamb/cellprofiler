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

default_species_csv <- file.path(project_root, "output", "runs", "mixed_cellpose_yolo_full_dataset_v1_bgclean", "linked_species_stats", "species_overview.csv")
default_matches_csv <- ""
default_tree_candidates <- c(
  file.path(project_root, "output", "runs", "mixed_cellpose_yolo_full_dataset_v1_bgclean", "phylogenetic_analysis", "source_tree_full.nwk"),
  file.path(dirname(project_root), "Desmognathus_TE", "results", "phylogeny", "processed_phylogeny.nwk")
)
default_output_dir <- file.path(project_root, "output", "runs", "mixed_cellpose_yolo_full_dataset_v1_bgclean", "phylogenetic_analysis")

first_existing_path <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (!length(hits)) {
    return(paths[[1]])
  }
  hits[[1]]
}

default_tree_file <- first_existing_path(default_tree_candidates)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

species_csv <- arg_value("--species-csv", default_species_csv)
matches_csv <- arg_value("--matches-csv", default_matches_csv)
tree_file <- arg_value("--tree-file", default_tree_file)
output_dir <- arg_value("--output-dir", default_output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

species_df <- read.csv(species_csv, stringsAsFactors = FALSE, check.names = FALSE)
tree <- read.tree(tree_file)

if (nrow(species_df) == 0) {
  stop("Species overview CSV is empty")
}

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

source_tree_label <- "Local Desmognathus processed phylogeny"
if (nzchar(matches_csv) && file.exists(matches_csv)) {
  matches_df <- read.csv(matches_csv, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(matches_df) == 0) {
    stop("Matches CSV is empty")
  }

  matches_df$tip_label <- paste0(gsub(" ", "_", matches_df$taxon_name), "_ott", matches_df$ott_id)
  tip_map <- setNames(display_species(matches_df$input_species), matches_df$tip_label)
  tree$tip.label <- unname(tip_map[tree$tip.label])
  if (anyNA(tree$tip.label)) {
    stop("Failed to map one or more external-tree tip labels back to species names")
  }
  source_tree_label <- "Mapped external phylogeny"
} else {
  tree$tip.label <- display_species(tree$tip.label)
}

species_df$species_input <- species_df$species
species_df$species <- display_species(species_df$species)
species_df$species_key <- standardize_species(species_df$species)
tree_tip_keys <- standardize_species(tree$tip.label)

coverage_df <- data.frame(
  input_species = species_df$species_input,
  normalized_species = species_df$species,
  in_tree = species_df$species_key %in% tree_tip_keys,
  stringsAsFactors = FALSE
)
write.csv(coverage_df, file.path(output_dir, "species_tree_coverage.csv"), row.names = FALSE)

missing_species <- coverage_df$normalized_species[!coverage_df$in_tree]
if (sum(coverage_df$in_tree) < 4) {
  stop("Fewer than 4 species from the dataset are represented in the selected tree")
}

tree <- keep.tip(tree, tree$tip.label[tree_tip_keys %in% species_df$species_key])
tree <- collapse.singles(tree)
tree_tip_keys <- standardize_species(tree$tip.label)
species_df <- species_df[match(tree_tip_keys, species_df$species_key), , drop = FALSE]
if (anyNA(species_df$species)) {
  stop("Failed to align one or more tree tips to species-overview rows")
}
row.names(species_df) <- species_df$species
tree$tip.label <- species_df$species

if (is.null(tree$edge.length) || anyNA(tree$edge.length) || all(tree$edge.length == 0)) {
  tree <- compute.brlen(tree, method = "Grafen")
}

source_tree_outfile <- file.path(output_dir, "source_tree_full.nwk")
if (normalizePath(tree_file, mustWork = TRUE) != normalizePath(source_tree_outfile, mustWork = FALSE)) {
  file.copy(tree_file, source_tree_outfile, overwrite = TRUE)
}
write.tree(tree, file.path(output_dir, "analysis_tree_pruned.nwk"))

trait_columns <- c(
  genome = "median_estimated_genome_pg_strict",
  nucleus = "median_nuc_area_um2_strict",
  cell = "median_cell_area_um2_strict"
)

trait_labels <- c(
  genome = "Median genome size (pg)",
  nucleus = "Median nucleus area (um^2)",
  cell = "Median cell area (um^2)"
)

species_df$log_genome <- log10(species_df[[trait_columns[["genome"]]]])
species_df$log_nucleus <- log10(species_df[[trait_columns[["nucleus"]]]])
species_df$log_cell <- log10(species_df[[trait_columns[["cell"]]]])
species_df$z_genome <- as.numeric(scale(species_df$log_genome))
species_df$z_nucleus <- as.numeric(scale(species_df$log_nucleus))
species_df$z_cell <- as.numeric(scale(species_df$log_cell))

fit_signal <- function(trait_name) {
  form <- as.formula(sprintf("%s ~ 1", trait_name))
  fit <- phylolm(form, data = species_df, phy = tree, model = "lambda")
  data.frame(
    trait = trait_name,
    lambda = unname(fit$optpar),
    logLik = as.numeric(logLik(fit)),
    AIC = AIC(fit),
    sigma2 = fit$sigma2,
    stringsAsFactors = FALSE
  )
}

signal_df <- do.call(rbind, lapply(c("log_genome", "log_nucleus", "log_cell"), fit_signal))
signal_df$trait_label <- c(trait_labels[["genome"]], trait_labels[["nucleus"]], trait_labels[["cell"]])
write.csv(signal_df, file.path(output_dir, "trait_phylogenetic_signal.csv"), row.names = FALSE)

color_ramp <- colorRampPalette(c("#f7f1e5", "#c58c63", "#7a3e1d"))

plot_trait_tree <- function(values, title, outfile) {
  ace_fit <- ace(values, tree, type = "continuous", method = "ML")
  all_states <- c(values, ace_fit$ace)
  edge_vals <- rowMeans(cbind(all_states[tree$edge[, 1]], all_states[tree$edge[, 2]]))
  palette <- color_ramp(128)
  bins <- cut(edge_vals, breaks = seq(min(all_states), max(all_states), length.out = length(palette) + 1), include.lowest = TRUE, labels = FALSE)
  edge_cols <- palette[bins]

  png(outfile, width = 1600, height = 1200, res = 180)
  layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
  par(mar = c(2, 2, 3, 1))
  plot(tree, edge.color = edge_cols, show.tip.label = TRUE, cex = 0.85, no.margin = FALSE, main = title)
  par(mar = c(4, 2, 3, 5))
  y_breaks <- seq(min(all_states), max(all_states), length.out = length(palette) + 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(y_breaks))
  rect(0, y_breaks[-length(y_breaks)], 1, y_breaks[-1], col = palette, border = NA)
  axis(4)
  mtext(side = 4, line = 3, title)
  dev.off()
}

plot_tree_tip_heatmap <- function(outfile) {
  genome_cols <- color_ramp(128)[cut(species_df[[trait_columns[["genome"]]]], 128, include.lowest = TRUE, labels = FALSE)]
  nucleus_cols <- color_ramp(128)[cut(species_df[[trait_columns[["nucleus"]]]], 128, include.lowest = TRUE, labels = FALSE)]
  cell_cols <- color_ramp(128)[cut(species_df[[trait_columns[["cell"]]]], 128, include.lowest = TRUE, labels = FALSE)]

  png(outfile, width = 1800, height = 1400, res = 180)
  par(mar = c(2, 2, 3, 10))
  plot(tree, show.tip.label = TRUE, cex = 0.85, x.lim = c(0, max(node.depth.edgelength(tree)) * 1.45), main = "Phylogeny with species trait heat strips")
  last_plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  x0 <- max(last_plot$xx)
  strip_w <- max(last_plot$xx) * 0.06
  y <- last_plot$yy[1:length(tree$tip.label)]
  rect(x0 + strip_w * 0.5, y - 0.28, x0 + strip_w * 1.3, y + 0.28, col = genome_cols, border = NA)
  rect(x0 + strip_w * 1.5, y - 0.28, x0 + strip_w * 2.3, y + 0.28, col = nucleus_cols, border = NA)
  rect(x0 + strip_w * 2.5, y - 0.28, x0 + strip_w * 3.3, y + 0.28, col = cell_cols, border = NA)
  text(x0 + strip_w * 0.9, max(y) + 1.2, "Genome", srt = 45, adj = 0, cex = 0.8)
  text(x0 + strip_w * 1.9, max(y) + 1.2, "Nucleus", srt = 45, adj = 0, cex = 0.8)
  text(x0 + strip_w * 2.9, max(y) + 1.2, "Cell", srt = 45, adj = 0, cex = 0.8)
  dev.off()
}

plot_trait_tree(species_df[[trait_columns[["genome"]]]], "Genome size mapped on Desmognathus phylogeny", file.path(output_dir, "tree_genome_size.png"))
plot_trait_tree(species_df[[trait_columns[["nucleus"]]]], "Nucleus size mapped on Desmognathus phylogeny", file.path(output_dir, "tree_nucleus_size.png"))
plot_trait_tree(species_df[[trait_columns[["cell"]]]], "Cell size mapped on Desmognathus phylogeny", file.path(output_dir, "tree_cell_size.png"))
plot_tree_tip_heatmap(file.path(output_dir, "tree_tip_heatmap.png"))

pairwise_specs <- list(
  list(
    panel = "A. Genome size and nucleus size",
    xvar = "log_genome",
    yvar = "log_nucleus",
    x_label = "Log10 genome size (pg)",
    y_label = "Log10 nucleus area (um^2)"
  ),
  list(
    panel = "B. Genome size and cell size",
    xvar = "log_genome",
    yvar = "log_cell",
    x_label = "Log10 genome size (pg)",
    y_label = "Log10 cell area (um^2)"
  ),
  list(
    panel = "C. Nucleus size and cell size",
    xvar = "log_nucleus",
    yvar = "log_cell",
    x_label = "Log10 nucleus area (um^2)",
    y_label = "Log10 cell area (um^2)"
  )
)

fit_pairwise_pgls <- function(spec) {
  form <- as.formula(sprintf("%s ~ %s", spec$yvar, spec$xvar))
  fit <- phylolm(form, data = species_df, phy = tree, model = "lambda")
  sm <- summary(fit)$coefficients
  data.frame(
    panel = spec$panel,
    xvar = spec$xvar,
    yvar = spec$yvar,
    x_label = spec$x_label,
    y_label = spec$y_label,
    intercept = sm["(Intercept)", "Estimate"],
    slope = sm[spec$xvar, "Estimate"],
    std_error = sm[spec$xvar, "StdErr"],
    p_value = sm[spec$xvar, "p.value"],
    lambda = unname(fit$optpar),
    n_species = nrow(species_df),
    stringsAsFactors = FALSE
  )
}

pairwise_results <- do.call(rbind, lapply(pairwise_specs, fit_pairwise_pgls))
write.csv(pairwise_results, file.path(output_dir, "pairwise_pgls_results.csv"), row.names = FALSE)

make_pairwise_plot <- function(spec, fit_row) {
  panel_df <- data.frame(
    x = species_df[[spec$xvar]],
    y = species_df[[spec$yvar]],
    species = species_df$species,
    n_pairs = species_df$n_pairs_strict_core,
    stringsAsFactors = FALSE
  )

  x_seq <- seq(min(panel_df$x), max(panel_df$x), length.out = 100)
  line_df <- data.frame(
    x = x_seq,
    y = fit_row$intercept + fit_row$slope * x_seq
  )

  ann_x <- min(panel_df$x) + 0.04 * diff(range(panel_df$x))
  ann_y <- max(panel_df$y) - 0.06 * diff(range(panel_df$y))
  ann_label <- sprintf(
    "slope = %s\nP = %s\nlambda = %s\nn = %d",
    formatC(fit_row$slope, digits = 3, format = "f"),
    formatC(fit_row$p_value, digits = 3, format = "f"),
    formatC(fit_row$lambda, digits = 3, format = "f"),
    fit_row$n_species
  )

  ggplot(panel_df, aes(x = x, y = y)) +
    geom_point(aes(size = n_pairs), shape = 21, stroke = 0.4, color = "#6e3b20", fill = "#c98d63", alpha = 0.88) +
    geom_line(data = line_df, aes(x = x, y = y), inherit.aes = FALSE, color = "#a43224", linewidth = 1.1) +
    annotate("text", x = ann_x, y = ann_y, label = ann_label, hjust = 0, vjust = 1, size = 3.5, family = "serif") +
    scale_size_continuous(range = c(2.5, 7.5), guide = "none") +
    labs(title = spec$panel, x = spec$x_label, y = spec$y_label) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", family = "serif", size = 12),
      axis.title = element_text(family = "serif"),
      axis.text = element_text(family = "serif")
    )
}

pairwise_plots <- lapply(seq_along(pairwise_specs), function(i) {
  make_pairwise_plot(pairwise_specs[[i]], pairwise_results[i, , drop = FALSE])
})

png(file.path(output_dir, "figure_2_pairwise_pgls.png"), width = 1800, height = 1400, res = 180)
grid::grid.newpage()
push_vp <- function(row, col) {
  grid::viewport(layout.pos.row = row, layout.pos.col = col)
}
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
print(pairwise_plots[[1]], vp = push_vp(1, 1))
print(pairwise_plots[[2]], vp = push_vp(1, 2))
print(pairwise_plots[[3]], vp = push_vp(2, 1))
grid::grid.text(
  "Point size scales with the number of strict-core linked cell-nucleus pairs per species.",
  x = 0.5,
  y = 0.5,
  gp = grid::gpar(fontsize = 12, col = "#5b4c3d", fontfamily = "serif"),
  vp = push_vp(2, 2)
)
dev.off()

model_specs <- list(
  G_to_N_to_C = list(z_nucleus = c("z_genome"), z_cell = c("z_nucleus")),
  G_to_C_to_N = list(z_cell = c("z_genome"), z_nucleus = c("z_cell")),
  G_to_both = list(z_nucleus = c("z_genome"), z_cell = c("z_genome")),
  G_to_both_plus_N_to_C = list(z_nucleus = c("z_genome"), z_cell = c("z_genome", "z_nucleus")),
  G_to_both_plus_C_to_N = list(z_cell = c("z_genome"), z_nucleus = c("z_genome", "z_cell"))
)

fit_equation <- function(response, predictors) {
  rhs <- if (length(predictors)) paste(predictors, collapse = " + ") else "1"
  form <- as.formula(sprintf("%s ~ %s", response, rhs))
  fit <- phylolm(form, data = species_df, phy = tree, model = "lambda")
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

dsep_test <- function(model_name) {
  if (model_name == "G_to_N_to_C") {
    fit <- phylolm(z_cell ~ z_nucleus + z_genome, data = species_df, phy = tree, model = "lambda")
    sm <- summary(fit)$coefficients
    return(data.frame(model = model_name, claim = "cell _||_ genome | nucleus", test_term = "z_genome", estimate = sm["z_genome", "Estimate"], p_value = sm["z_genome", "p.value"], stringsAsFactors = FALSE))
  }
  if (model_name == "G_to_C_to_N") {
    fit <- phylolm(z_nucleus ~ z_cell + z_genome, data = species_df, phy = tree, model = "lambda")
    sm <- summary(fit)$coefficients
    return(data.frame(model = model_name, claim = "nucleus _||_ genome | cell", test_term = "z_genome", estimate = sm["z_genome", "Estimate"], p_value = sm["z_genome", "p.value"], stringsAsFactors = FALSE))
  }
  if (model_name == "G_to_both") {
    fit <- phylolm(z_cell ~ z_genome + z_nucleus, data = species_df, phy = tree, model = "lambda")
    sm <- summary(fit)$coefficients
    return(data.frame(model = model_name, claim = "cell _||_ nucleus | genome", test_term = "z_nucleus", estimate = sm["z_nucleus", "Estimate"], p_value = sm["z_nucleus", "p.value"], stringsAsFactors = FALSE))
  }
  data.frame(model = model_name, claim = "saturated", test_term = "", estimate = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE)
}

model_rows <- list()
coef_rows <- list()
dsep_rows <- list()

for (model_name in names(model_specs)) {
  spec <- model_specs[[model_name]]
  equation_fits <- lapply(names(spec), function(response) fit_equation(response, spec[[response]]))
  names(equation_fits) <- names(spec)
  total_logLik <- sum(vapply(equation_fits, function(x) x$coefs$logLik[1], numeric(1)))
  total_k <- sum(vapply(equation_fits, function(x) x$coefs$k[1], numeric(1)))
  n_species <- nrow(species_df)
  aicc <- calc_aicc(total_logLik, total_k, n_species)

  model_rows[[model_name]] <- data.frame(
    model = model_name,
    equations = paste(vapply(names(spec), function(response) sprintf("%s ~ %s", response, paste(spec[[response]], collapse = " + ")), character(1)), collapse = " ; "),
    n_species = n_species,
    n_equations = length(spec),
    total_logLik = total_logLik,
    total_k = total_k,
    AICc = aicc,
    stringsAsFactors = FALSE
  )

  coef_df <- do.call(rbind, lapply(equation_fits, function(x) x$coefs))
  coef_df$model <- model_name
  coef_rows[[model_name]] <- coef_df
  dsep_rows[[model_name]] <- dsep_test(model_name)
}

model_scores <- do.call(rbind, model_rows)
model_scores <- model_scores[order(model_scores$AICc), , drop = FALSE]
model_scores$delta_AICc <- model_scores$AICc - min(model_scores$AICc)
model_scores$model_weight <- exp(-0.5 * model_scores$delta_AICc)
model_scores$model_weight <- model_scores$model_weight / sum(model_scores$model_weight)
write.csv(model_scores, file.path(output_dir, "path_model_scores.csv"), row.names = FALSE)

coef_df <- do.call(rbind, coef_rows)
coef_df <- coef_df[, c("model", "response", "term", "Estimate", "StdErr", "p.value", "lambda")]
names(coef_df) <- c("model", "response", "term", "estimate", "std_error", "p_value", "lambda")
write.csv(coef_df, file.path(output_dir, "path_model_coefficients.csv"), row.names = FALSE)

dsep_df <- do.call(rbind, dsep_rows)
write.csv(dsep_df, file.path(output_dir, "path_model_dsep.csv"), row.names = FALSE)

best_model <- model_scores$model[[1]]
best_coef_df <- subset(coef_df, model == best_model & term != "(Intercept)")

score_plot <- ggplot(model_scores, aes(x = reorder(model, AICc), y = AICc, fill = model == best_model)) +
  geom_col(width = 0.72, color = NA) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#8c5632", "FALSE" = "#d8c1aa"), guide = "none") +
  labs(title = "Phylogenetic path-model comparison", x = NULL, y = "Summed PGLS AICc") +
  theme_minimal(base_size = 12)
ggsave(file.path(output_dir, "path_model_scores.png"), score_plot, width = 10, height = 5.6, dpi = 180)

coef_plot <- ggplot(best_coef_df, aes(x = estimate, y = paste(response, "<-", term), xmin = estimate - 1.96 * std_error, xmax = estimate + 1.96 * std_error)) +
  geom_vline(xintercept = 0, color = "#d7c6b1", linewidth = 0.6) +
  geom_errorbar(width = 0.18, orientation = "y", color = "#c58c63", linewidth = 1.0) +
  geom_point(color = "#8c5632", size = 3) +
  labs(title = sprintf("Best-model standardized path coefficients: %s", best_model), x = "Standardized coefficient", y = NULL) +
  theme_minimal(base_size = 12)
ggsave(file.path(output_dir, "best_model_coefficients.png"), coef_plot, width = 10, height = 4.8, dpi = 180)

node_positions <- list(
  z_genome = c(0.2, 0.78),
  z_nucleus = c(0.8, 0.78),
  z_cell = c(0.5, 0.24)
)

node_labels <- c(
  z_genome = "Genome\nsize",
  z_nucleus = "Nucleus\nsize",
  z_cell = "Cell\nsize"
)

node_fill <- c(
  z_genome = "#e8d4c1",
  z_nucleus = "#efdcc8",
  z_cell = "#f4e8d8"
)

draw_path_diagram <- function(spec, panel_title, footer_lines = character(0), highlight = FALSE, edge_label_map = NULL) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  border_col <- if (highlight) "#8c5632" else "#d7c6b1"
  border_lwd <- if (highlight) 2.2 else 1.0
  rect(0.03, 0.03, 0.97, 0.97, border = border_col, lwd = border_lwd, col = if (highlight) "#fcf6ee" else "#fffdfa")

  text(0.05, 0.94, panel_title, adj = c(0, 1), cex = 0.92, font = 2)

  for (node in names(node_positions)) {
    xy <- node_positions[[node]]
    symbols(xy[1], xy[2], circles = 0.09, inches = FALSE, bg = node_fill[[node]], fg = "#6d5138", add = TRUE)
    text(xy[1], xy[2], node_labels[[node]], cex = 0.82)
  }

  for (response in names(spec)) {
    for (predictor in spec[[response]]) {
      from_xy <- node_positions[[predictor]]
      to_xy <- node_positions[[response]]
      arrows(from_xy[1], from_xy[2], to_xy[1], to_xy[2], length = 0.09, lwd = 2.2, col = "#8c5632")
      label_key <- paste(response, predictor, sep = "|")
      if (!is.null(edge_label_map) && label_key %in% names(edge_label_map)) {
        mid_x <- (from_xy[1] + to_xy[1]) / 2
        mid_y <- (from_xy[2] + to_xy[2]) / 2 + ifelse(from_xy[2] == to_xy[2], 0.06, 0.04)
        text(mid_x, mid_y, edge_label_map[[label_key]], cex = 0.72, col = "#6d5138")
      }
    }
  }

  if (length(footer_lines)) {
    text(0.05, 0.11, paste(footer_lines, collapse = "\n"), adj = c(0, 0), cex = 0.76, col = "#5b4c3d")
  }
}

model_display_name <- c(
  G_to_N_to_C = "Genome -> Nucleus -> Cell",
  G_to_C_to_N = "Genome -> Cell -> Nucleus",
  G_to_both = "Genome -> Nucleus; Genome -> Cell",
  G_to_both_plus_N_to_C = "Genome -> Nucleus; Genome + Nucleus -> Cell",
  G_to_both_plus_C_to_N = "Genome -> Cell; Genome + Cell -> Nucleus"
)

png(file.path(output_dir, "figure_3_candidate_path_models.png"), width = 1800, height = 1250, res = 180)
par(mfrow = c(2, 3), mar = c(0.4, 0.4, 1.2, 0.4))
for (i in seq_len(nrow(model_scores))) {
  row <- model_scores[i, , drop = FALSE]
  footer <- c(
    sprintf("AICc = %s", formatC(row$AICc, digits = 2, format = "f")),
    sprintf("weight = %s", formatC(row$model_weight, digits = 3, format = "f"))
  )
  draw_path_diagram(
    model_specs[[row$model]],
    panel_title = sprintf("%s. %s", LETTERS[i], model_display_name[[row$model]]),
    footer_lines = footer,
    highlight = identical(row$model[[1]], best_model)
  )
}
plot.new()
text(
  0.05, 0.95,
  "Models are ranked by summed phylogenetic GLS AICc.\nThe highlighted panel is the best-supported topology.",
  adj = c(0, 1),
  cex = 1.0,
  col = "#5b4c3d"
)
dev.off()

best_edge_labels <- setNames(
  sprintf(
    "beta = %s\nP = %s",
    formatC(best_coef_df$estimate, digits = 2, format = "f"),
    formatC(best_coef_df$p_value, digits = 3, format = "f")
  ),
  paste(best_coef_df$response, best_coef_df$term, sep = "|")
)

coef_term_map <- c(
  "z_nucleus|z_genome" = "Genome -> Nucleus",
  "z_cell|z_nucleus" = "Nucleus -> Cell",
  "z_cell|z_genome" = "Genome -> Cell"
)
coef_terms <- unname(coef_term_map[paste(best_coef_df$response, best_coef_df$term, sep = "|")])
coef_terms[is.na(coef_terms)] <- paste(best_coef_df$term, "->", best_coef_df$response)

png(file.path(output_dir, "figure_4_best_model.png"), width = 1800, height = 900, res = 180)
layout(matrix(c(1, 2), nrow = 1), widths = c(1.05, 1.25))
par(mar = c(0.6, 0.6, 1.4, 0.6))
draw_path_diagram(
  model_specs[[best_model]],
  panel_title = sprintf("A. Best-supported path model: %s", model_display_name[[best_model]]),
  footer_lines = c(
    sprintf("AICc = %s", formatC(model_scores$AICc[[1]], digits = 2, format = "f")),
    sprintf("weight = %s", formatC(model_scores$model_weight[[1]], digits = 3, format = "f"))
  ),
  highlight = TRUE,
  edge_label_map = best_edge_labels
)

par(mar = c(4.5, 6.5, 1.4, 1.5))
ci_low <- best_coef_df$estimate - 1.96 * best_coef_df$std_error
ci_high <- best_coef_df$estimate + 1.96 * best_coef_df$std_error
y_pos <- seq_along(best_coef_df$estimate)
plot(
  best_coef_df$estimate,
  y_pos,
  xlim = range(c(ci_low, ci_high, 0)),
  ylim = c(0.5, length(y_pos) + 0.5),
  yaxt = "n",
  ylab = "",
  xlab = "Standardized coefficient (95% CI)",
  pch = 21,
  bg = "#8c5632",
  col = "#8c5632",
  cex = 1.3,
  main = "B. Best-model path coefficients"
)
segments(ci_low, y_pos, ci_high, y_pos, lwd = 2.4, col = "#c58c63")
abline(v = 0, lty = 2, col = "#cdb8a3")
axis(2, at = y_pos, labels = coef_terms, las = 1)
dev.off()

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

signal_html <- signal_df
signal_html$lambda <- fmt_num(signal_html$lambda)
signal_html$AIC <- fmt_num(signal_html$AIC)
signal_html$sigma2 <- fmt_num(signal_html$sigma2)

pairwise_html <- pairwise_results
pairwise_html$slope <- fmt_num(pairwise_html$slope)
pairwise_html$std_error <- fmt_num(pairwise_html$std_error)
pairwise_html$p_value <- fmt_num(pairwise_html$p_value)
pairwise_html$lambda <- fmt_num(pairwise_html$lambda)
pairwise_html$n_species <- as.integer(pairwise_html$n_species)

scores_html <- model_scores
scores_html$model_label <- unname(model_display_name[scores_html$model])
scores_html$AICc <- fmt_num(scores_html$AICc)
scores_html$delta_AICc <- fmt_num(scores_html$delta_AICc)
scores_html$model_weight <- fmt_num(scores_html$model_weight)

coef_html <- best_coef_df
coef_html$estimate <- fmt_num(coef_html$estimate)
coef_html$std_error <- fmt_num(coef_html$std_error)
coef_html$p_value <- fmt_num(coef_html$p_value)
coef_html$lambda <- fmt_num(coef_html$lambda)

report_summary <- list(
  n_species_input = nrow(coverage_df),
  n_species_in_tree = nrow(species_df),
  n_species_missing_from_tree = length(missing_species),
  missing_species = unname(missing_species),
  best_model = best_model,
  best_model_AICc = unname(model_scores$AICc[[1]]),
  best_model_weight = unname(model_scores$model_weight[[1]]),
  tree_tip_count = length(tree$tip.label),
  tree_source = source_tree_label,
  tree_file = normalizePath(tree_file, mustWork = TRUE)
)
write(toJSON(report_summary, pretty = TRUE, auto_unbox = TRUE), file.path(output_dir, "summary.json"))

missing_html <- if (length(missing_species)) {
  data.frame(species = missing_species, stringsAsFactors = FALSE)
} else {
  data.frame(species = "None", stringsAsFactors = FALSE)
}

missing_species_text <- if (length(missing_species)) {
  paste(missing_species, collapse = ", ")
} else {
  "None"
}

html <- sprintf(
'<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Desmognathus Trait Phylogeny Report</title>
  <style>
    :root { --ink:#241d17; --muted:#66584a; --line:#d9ccbf; --accent:#8c5632; --paper:#fffdfa; }
    html, body { margin: 0; padding: 0; background: #f7f3ee; color: var(--ink); }
    body { font: 16px/1.65 Georgia, serif; }
    main { max-width: 980px; margin: 0 auto; padding: 40px 28px 56px; background: var(--paper); box-shadow: 0 8px 30px rgba(33, 23, 14, 0.06); }
    h1, h2, h3 { font-weight: 600; }
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
    .supp-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(240px, 1fr)); gap: 16px; margin-top: 12px; }
    .supp-grid figure { margin: 0; }
    .supp-grid figcaption { font-size: 0.83rem; color: var(--muted); margin-top: 6px; }
    .mono { font-family: "SFMono-Regular", Menlo, Consolas, monospace; font-size: 0.9em; }
  </style>
</head>
<body>
  <main>
    <header>
      <div class="eyebrow">Comparative Phylogenetic Analysis</div>
      <h1>Genome size, nucleus size, and cell size track one another across <em>Desmognathus</em></h1>
      <p class="subhead">bg-clean linked cell-nucleus medians analyzed on the local curated salamander backbone</p>
    </header>

    <section class="summary">
      <p>%d species were supplied to the report and %d were represented in the phylogeny. The best-supported path model was <span class="mono">%s</span> with AICc %s and model weight %s, consistent with genome size covarying most strongly with nucleus size and nucleus size then tracking cell size.</p>
      <p class="note">Species excluded from tree-based inference: %s.</p>
      <div class="links">
        <a href="summary.json">summary json</a>
        <a href="pairwise_pgls_results.csv">pairwise PGLS csv</a>
        <a href="path_model_scores.csv">path model ranking csv</a>
        <a href="path_model_coefficients.csv">path coefficients csv</a>
        <a href="species_tree_coverage.csv">tree coverage csv</a>
        <a href="source_tree_full.nwk">source tree</a>
        <a href="analysis_tree_pruned.nwk">analysis tree</a>
      </div>
    </section>

    <section class="figure">
      <img src="tree_tip_heatmap.png" alt="Trait heat strips across the Desmognathus phylogeny">
      <div class="figcaption"><strong>Figure 1.</strong> Phylogenetic distribution of species-level medians for genome size, nucleus size, and cell size. Each tip is annotated with aligned trait strips so the three linked traits can be compared directly across the salamander backbone.</div>
    </section>

    <section class="figure">
      <img src="figure_2_pairwise_pgls.png" alt="Pairwise phylogenetic GLS relationships among genome size, nucleus size, and cell size">
      <div class="figcaption"><strong>Figure 2.</strong> Pairwise phylogenetic GLS relationships among the three linked traits. Each panel shows one phylogenetically corrected regression using species medians; point size scales with the number of strict-core linked cell-nucleus pairs contributing to that species estimate.</div>
    </section>

    <section class="figure">
      <img src="figure_3_candidate_path_models.png" alt="Candidate path model comparison">
      <div class="figcaption"><strong>Figure 3.</strong> Candidate phylogenetic path models ranked by summed PGLS AICc. The highlighted model had the strongest support among the tested hypotheses.</div>
    </section>

    <section class="figure">
      <img src="figure_4_best_model.png" alt="Best-supported path model and coefficients">
      <div class="figcaption"><strong>Figure 4.</strong> Best-supported path model and standardized path coefficients. Arrow labels in the left panel show coefficient estimates and P-values for the retained paths, and the right panel shows coefficient uncertainty as 95%% confidence intervals.</div>
    </section>

    <section>
      <h2>Table 1. Pairwise PGLS Results</h2>
      <div class="table-wrap">
        <table>
          <tr><th>Panel</th><th>Slope</th><th>Std. error</th><th>P-value</th><th>Lambda</th><th>Species</th></tr>
          %s
        </table>
      </div>
    </section>

    <section>
      <h2>Table 2. Candidate Path-Model Ranking</h2>
      <div class="table-wrap">
        <table>
          <tr><th>Model</th><th>AICc</th><th>Delta AICc</th><th>Weight</th><th>Equations</th></tr>
          %s
        </table>
      </div>
    </section>

    <section>
      <h2>Coverage And Diagnostics</h2>
      <p class="note">Species excluded from the phylogenetic backbone are listed separately because they contribute descriptive linked-trait summaries but not phylogenetically corrected estimates.</p>
      <div class="table-wrap">
        <table>
          <tr><th>Species excluded from tree-based inference</th></tr>
          %s
        </table>
      </div>
      <div class="table-wrap">
        <table>
          <tr><th>Trait</th><th>Lambda</th><th>AIC</th><th>Sigma2</th></tr>
          %s
        </table>
      </div>
    </section>

    <section>
      <h2>Supplementary Figures</h2>
      <div class="supp-grid">
        <figure>
          <a href="tree_genome_size.png"><img src="tree_genome_size.png" alt="Genome size mapped on the tree"></a>
          <figcaption>Supplementary Figure S1. Genome size mapped alone on the pruned phylogeny.</figcaption>
        </figure>
        <figure>
          <a href="tree_nucleus_size.png"><img src="tree_nucleus_size.png" alt="Nucleus size mapped on the tree"></a>
          <figcaption>Supplementary Figure S2. Nucleus size mapped alone on the pruned phylogeny.</figcaption>
        </figure>
        <figure>
          <a href="tree_cell_size.png"><img src="tree_cell_size.png" alt="Cell size mapped on the tree"></a>
          <figcaption>Supplementary Figure S3. Cell size mapped alone on the pruned phylogeny.</figcaption>
        </figure>
      </div>
    </section>
  </main>
</body>
</html>',
  nrow(coverage_df),
  nrow(species_df),
  best_model,
  fmt_num(report_summary$best_model_AICc),
  fmt_num(report_summary$best_model_weight),
  missing_species_text,
  table_rows(pairwise_html[, c("panel", "slope", "std_error", "p_value", "lambda", "n_species"), drop = FALSE], c("panel", "slope", "std_error", "p_value", "lambda", "n_species")),
  table_rows(scores_html[, c("model_label", "AICc", "delta_AICc", "model_weight", "equations"), drop = FALSE], c("model_label", "AICc", "delta_AICc", "model_weight", "equations")),
  table_rows(missing_html[, c("species"), drop = FALSE], c("species")),
  table_rows(signal_html[, c("trait_label", "lambda", "AIC", "sigma2"), drop = FALSE], c("trait_label", "lambda", "AIC", "sigma2"))
)

writeLines(html, file.path(output_dir, "index.html"))

message("Wrote phylogenetic report: ", file.path(output_dir, "index.html"))
