# =============================================================================
# visualize_results.R
# Generates plots from the most recent benchmark run.
# Can be sourced standalone to regenerate plots without re-running benchmarks:
#
#   source("benchmarks/setup.R")     # sets results_dir to latest automatically
#   source("benchmarks/visualize_results.R")
# =============================================================================

library(ggplot2)

# If sourced standalone, pick the most recent results folder
if (!exists("results_dir")) {
  dirs       <- list.dirs("benchmarks/results", recursive = FALSE)
  results_dir <- dirs[which.max(file.info(dirs)$mtime)]
  cat("Auto-detected latest results:", results_dir, "\n\n")
  load(file.path(results_dir, "benchmark_data.RData"))
} else if (!exists("lizards_bench")) {
  load(file.path(results_dir, "benchmark_data.RData"))
}

plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

theme_bench <- theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        plot.subtitle = element_text(size = 10, color = "grey40"))

version_colours <- c("original" = "#E41A1C", "new" = "#4DAF4A")

# Helper: extract bench::mark timings into a tidy data frame (nanoseconds -> s)
bench_to_df <- function(bench_obj) {
  data.frame(
    time_s  = c(as.numeric(bench_obj$time[[1]]) / 1e9,
                as.numeric(bench_obj$time[[2]]) / 1e9),
    version = factor(c(rep("original", length(bench_obj$time[[1]])),
                       rep("new",      length(bench_obj$time[[2]]))),
                     levels = c("original", "new"))
  )
}

# Helper: extract microbenchmark timings into a tidy data frame
mb_to_df <- function(mb_obj) {
  data.frame(
    time_s  = mb_obj$time / 1e9,
    version = factor(mb_obj$expr, levels = c("original", "new"))
  )
}

# ====== Plot 1: Lizards boxplot ======
cat("Plot 1: Lizards comparison...\n")
tryCatch({
  df <- mb_to_df(lizards_bench)
  n  <- nrow(lizards_bench) / 2
  p  <- ggplot(df, aes(x = version, y = time_s, fill = version)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
    scale_fill_manual(values = version_colours) +
    labs(title    = "Lizards Dataset: Original vs New",
         subtitle = paste0(n, " evaluations each (interleaved)"),
         x = NULL, y = "Time (seconds)") +
    theme_bench
  ggsave(file.path(plots_dir, "01_lizards_comparison.png"), p,
         width = 7, height = 5, dpi = 120)
  cat("  Saved.\n")
}, error = function(e) cat("  Error:", conditionMessage(e), "\n"))

# ====== Plot 2: Endometrial boxplot ======
cat("Plot 2: Endometrial comparison...\n")
tryCatch({
  df <- mb_to_df(endo_bench)
  n  <- nrow(endo_bench) / 2
  p  <- ggplot(df, aes(x = version, y = time_s, fill = version)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
    scale_fill_manual(values = version_colours) +
    labs(title    = "Endometrial Dataset: Original vs New",
         subtitle = paste0(n, " evaluations each (interleaved)"),
         x = NULL, y = "Time (seconds)") +
    theme_bench
  ggsave(file.path(plots_dir, "02_endometrial_comparison.png"), p,
         width = 7, height = 5, dpi = 120)
  cat("  Saved.\n")
}, error = function(e) cat("  Error:", conditionMessage(e), "\n"))

# ====== Plot 3: MultipleFeatures brglmFit boxplot ======
cat("Plot 3: MultipleFeatures (brglmFit) comparison...\n")
tryCatch({
  df <- bench_to_df(mf_bench)
  n  <- length(mf_bench$time[[1]])
  p  <- ggplot(df, aes(x = version, y = time_s, fill = version)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
    scale_fill_manual(values = version_colours) +
    labs(title    = "MultipleFeatures (brglmFit): Original vs New",
         subtitle = paste0(n, " iterations each (interleaved via bench::mark)"),
         x = NULL, y = "Time (seconds)") +
    theme_bench
  ggsave(file.path(plots_dir, "03_mf_brglmfit_comparison.png"), p,
         width = 7, height = 5, dpi = 120)
  cat("  Saved.\n")
}, error = function(e) cat("  Error:", conditionMessage(e), "\n"))

# ====== Plot 4: MultipleFeatures mdyplFit boxplot ======
cat("Plot 4: MultipleFeatures (mdyplFit) comparison...\n")
tryCatch({
  df <- bench_to_df(mdypl_bench)
  n  <- length(mdypl_bench$time[[1]])
  p  <- ggplot(df, aes(x = version, y = time_s, fill = version)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4) +
    scale_fill_manual(values = version_colours) +
    labs(title    = "MultipleFeatures (mdyplFit): Original vs New",
         subtitle = paste0(n, " iterations each (interleaved via bench::mark)"),
         x = NULL, y = "Time (seconds)") +
    theme_bench
  ggsave(file.path(plots_dir, "04_mf_mdyplfit_comparison.png"), p,
         width = 7, height = 5, dpi = 120)
  cat("  Saved.\n")
}, error = function(e) cat("  Error:", conditionMessage(e), "\n"))

# ====== Plot 5: Speedup summary bar chart ======
cat("Plot 5: Speedup summary...\n")
tryCatch({
  clean_names <- gsub("\\.elapsed$", "", names(speedups))
  speedup_df  <- data.frame(
    Test    = factor(clean_names, levels = clean_names),
    Speedup = as.numeric(speedups),
    Size    = c("Small", "Medium", "Large", "Large")
  )
  p <- ggplot(speedup_df, aes(x = Test, y = Speedup, fill = Size)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.2fx", Speedup)), vjust = -0.4, size = 3.5) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "red", linewidth = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, max(speedup_df$Speedup) * 1.18) +
    labs(title    = "Speedup: New vs Original",
         subtitle = "Dashed line = no improvement (1.0×).  Higher is better.",
         x = NULL, y = "Speedup factor") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          legend.position = "top")
  ggsave(file.path(plots_dir, "05_speedup_summary.png"), p,
         width = 8, height = 5, dpi = 120)
  cat("  Saved.\n")
}, error = function(e) cat("  Error:", conditionMessage(e), "\n"))

# ====== Plot 6: Lizards density ======
cat("Plot 6: Lizards distribution...\n")
tryCatch({
  df      <- mb_to_df(lizards_bench)
  medians <- aggregate(time_s ~ version, df, median)
  p <- ggplot(df, aes(x = time_s, fill = version)) +
    geom_density(alpha = 0.55) +
    geom_vline(data = medians, aes(xintercept = time_s, colour = version),
               linetype = "dashed", linewidth = 0.9) +
    scale_fill_manual(values   = version_colours) +
    scale_colour_manual(values = version_colours) +
    labs(title    = "Execution Time Distribution — Lizards",
         subtitle = "Dashed lines = medians",
         x = "Time (seconds)", y = "Density") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  ggsave(file.path(plots_dir, "06_lizards_distribution.png"), p,
         width = 7, height = 5, dpi = 120)
  cat("  Saved.\n")
}, error = function(e) cat("  Error:", conditionMessage(e), "\n"))

cat("\nAll plots saved to:", plots_dir, "\n\n")
cat("Files generated:\n")
cat("  01_lizards_comparison.png\n")
cat("  02_endometrial_comparison.png\n")
cat("  03_mf_brglmfit_comparison.png\n")
cat("  04_mf_mdyplfit_comparison.png\n")
cat("  05_speedup_summary.png\n")
cat("  06_lizards_distribution.png\n")