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

scores_html <- model_scores
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

html <- sprintf(
'<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Phylogenetic Trait and Path Report</title>
  <style>
    :root { --bg:#f3ecdf; --panel:#fffcf8; --ink:#2e261d; --muted:#70604d; --line:#d7c6b1; --accent:#8c5632; }
    body { margin:0; color:var(--ink); background:linear-gradient(180deg,#f8f3eb 0%%, var(--bg) 100%%); font:15px/1.55 Georgia, serif; }
    .wrap { max-width:1520px; margin:0 auto; padding:24px; }
    .hero, .card, table { background:var(--panel); border:1px solid var(--line); border-radius:18px; box-shadow:0 12px 28px rgba(43,31,19,0.06); }
    .hero, .card { padding:18px; }
    .grid, .fig-grid { display:grid; gap:14px; margin:18px 0 24px; }
    .grid { grid-template-columns:repeat(auto-fit, minmax(220px, 1fr)); }
    .fig-grid { grid-template-columns:repeat(auto-fit, minmax(320px, 1fr)); }
    .metric-label { color:var(--muted); text-transform:uppercase; letter-spacing:.08em; font-size:12px; }
    .metric-value { margin-top:10px; color:var(--accent); font-size:28px; }
    .links a { display:inline-block; margin:0 10px 10px 0; padding:8px 12px; border-radius:999px; background:#efe1d2; border:1px solid #dfc9b5; color:var(--ink); text-decoration:none; }
    img { width:100%%; height:auto; display:block; border-radius:12px; border:1px solid #e4d7ca; background:#faf5ef; }
    table { width:100%%; border-collapse:collapse; overflow:hidden; margin-bottom:24px; }
    th, td { padding:8px 10px; border-bottom:1px solid var(--line); text-align:left; font-size:13px; }
    th { background:#f1e3d5; }
    .note { color:var(--muted); }
  </style>
</head>
<body>
  <div class="wrap">
    <section class="hero">
      <h1>Phylogenetic Trait and Path Report</h1>
      <p>Species-level medians from the bg-clean linked cell-nucleus dataset are mapped onto the existing local <em>Desmognathus</em> backbone tree used in the companion phylogeny workflow. This report therefore reuses that curated salamander tree rather than the incomplete OpenTree fallback.</p>
      <p class="note">Traits analyzed here are median strict-core species values for genome size, nucleus size, and cell size. Path models are compared as summed phylogenetic generalized least squares regressions with Pagel\'s lambda estimated per equation. Species absent from the local backbone are listed below and excluded from tree-based inference.</p>
      <div class="links">
        <a href="summary.json">summary json</a>
        <a href="species_tree_coverage.csv">tree coverage</a>
        <a href="trait_phylogenetic_signal.csv">trait phylogenetic signal</a>
        <a href="path_model_scores.csv">path model scores</a>
        <a href="path_model_coefficients.csv">path coefficients</a>
        <a href="path_model_dsep.csv">path d-sep checks</a>
        <a href="source_tree_full.nwk">source tree</a>
        <a href="analysis_tree_pruned.nwk">analysis tree</a>
      </div>
    </section>
    <section class="grid">
      <div class="card"><div class="metric-label">Input species</div><div class="metric-value">%d</div></div>
      <div class="card"><div class="metric-label">Species in tree</div><div class="metric-value">%d</div></div>
      <div class="card"><div class="metric-label">Excluded from tree</div><div class="metric-value">%d</div></div>
      <div class="card"><div class="metric-label">Best path model</div><div class="metric-value" style="font-size:22px">%s</div></div>
      <div class="card"><div class="metric-label">Best model AICc</div><div class="metric-value">%s</div></div>
      <div class="card"><div class="metric-label">Best model weight</div><div class="metric-value">%s</div></div>
    </section>
    <section class="fig-grid">
      <div class="card"><a href="tree_tip_heatmap.png"><img src="tree_tip_heatmap.png" alt="Trait heatmap tree"></a></div>
      <div class="card"><a href="tree_genome_size.png"><img src="tree_genome_size.png" alt="Genome tree"></a></div>
      <div class="card"><a href="tree_nucleus_size.png"><img src="tree_nucleus_size.png" alt="Nucleus tree"></a></div>
      <div class="card"><a href="tree_cell_size.png"><img src="tree_cell_size.png" alt="Cell tree"></a></div>
      <div class="card"><a href="path_model_scores.png"><img src="path_model_scores.png" alt="Path model scores"></a></div>
      <div class="card"><a href="best_model_coefficients.png"><img src="best_model_coefficients.png" alt="Best model coefficients"></a></div>
    </section>
    <section>
      <h2>Species excluded from phylogenetic inference</h2>
      <table>
        <tr><th>Species</th></tr>
        %s
      </table>
    </section>
    <section>
      <h2>Trait phylogenetic signal</h2>
      <table>
        <tr><th>Trait</th><th>Lambda</th><th>AIC</th><th>Sigma2</th></tr>
        %s
      </table>
    </section>
    <section>
      <h2>Path model comparison</h2>
      <table>
        <tr><th>Model</th><th>Equations</th><th>AICc</th><th>Delta AICc</th><th>Weight</th></tr>
        %s
      </table>
    </section>
    <section>
      <h2>Best-model standardized coefficients</h2>
      <table>
        <tr><th>Response</th><th>Predictor</th><th>Estimate</th><th>Std error</th><th>p-value</th><th>Lambda</th></tr>
        %s
      </table>
    </section>
  </div>
</body>
</html>',
  nrow(coverage_df),
  nrow(species_df),
  length(missing_species),
  best_model,
  fmt_num(report_summary$best_model_AICc),
  fmt_num(report_summary$best_model_weight),
  table_rows(missing_html[, c("species"), drop = FALSE], c("species")),
  table_rows(signal_html[, c("trait_label", "lambda", "AIC", "sigma2"), drop = FALSE], c("trait_label", "lambda", "AIC", "sigma2")),
  table_rows(scores_html[, c("model", "equations", "AICc", "delta_AICc", "model_weight"), drop = FALSE], c("model", "equations", "AICc", "delta_AICc", "model_weight")),
  table_rows(coef_html[, c("response", "term", "estimate", "std_error", "p_value", "lambda"), drop = FALSE], c("response", "term", "estimate", "std_error", "p_value", "lambda"))
)

writeLines(html, file.path(output_dir, "index.html"))

message("Wrote phylogenetic report: ", file.path(output_dir, "index.html"))
