# Fixed Benchmark Visualization Script
# This version properly handles the benchmark data structure

library(ggplot2)
library(microbenchmark)

# Find the most recent results directory
results_base <- "benchmarks/results"
dirs <- list.dirs(results_base, recursive = FALSE)
latest_dir <- dirs[which.max(file.info(dirs)$mtime)]

cat("Loading benchmark data from:", latest_dir, "\n")

# Load the benchmark data
load(file.path(latest_dir, "benchmark_data.RData"))

# Create plots directory
plots_dir <- file.path(latest_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

cat("\n=== Checking loaded data ===\n")
cat("lizards_bench dimensions:", nrow(lizards_bench), "rows\n")
cat("endo_bench dimensions:", nrow(endo_bench), "rows\n")
cat("Speedups vector:", names(speedups), "\n")
cat("Speedup values:", speedups, "\n\n")

# ====== Plot 1: Microbenchmark comparison (Lizards) ======
cat("Creating Plot 1: Lizards comparison...\n")
tryCatch({
  png(file.path(plots_dir, "01_lizards_comparison.png"), 
      width = 800, height = 600, res = 100)
  
  # Manual boxplot approach (more reliable than autoplot)
  lizards_df <- data.frame(
    time_seconds = lizards_bench$time / 1e9,
    version = lizards_bench$expr
  )
  
  p1 <- ggplot(lizards_df, aes(x = version, y = time_seconds, fill = version)) +
    geom_boxplot() +
    labs(title = "Lizards Dataset: Original vs New Implementation",
         subtitle = "100 evaluations each",
         x = "Version",
         y = "Time (seconds)") +
    scale_fill_manual(values = c("original" = "#E41A1C", "new" = "#4DAF4A")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  print(p1)
  dev.off()
  cat("Plot 1 saved\n")
}, error = function(e) {
  cat("Error in Plot 1:", conditionMessage(e), "\n")
  dev.off()
})

# ====== Plot 2: Microbenchmark comparison (Endometrial) ======
cat("Creating Plot 2: Endometrial comparison...\n")
tryCatch({
  png(file.path(plots_dir, "02_endometrial_comparison.png"), 
      width = 800, height = 600, res = 100)
  
  endo_df <- data.frame(
    time_seconds = endo_bench$time / 1e9,
    version = endo_bench$expr
  )
  
  p2 <- ggplot(endo_df, aes(x = version, y = time_seconds, fill = version)) +
    geom_boxplot() +
    labs(title = "Endometrial Dataset: Original vs New Implementation",
         subtitle = "50 evaluations each",
         x = "Version",
         y = "Time (seconds)") +
    scale_fill_manual(values = c("original" = "#E41A1C", "new" = "#4DAF4A")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  
  print(p2)
  dev.off()
  cat("Plot 2 saved\n")
}, error = function(e) {
  cat("Error in Plot 2:", conditionMessage(e), "\n")
  dev.off()
})

# ====== Plot 3: Speedup summary bar chart ======
cat("Creating Plot 3: Speedup summary...\n")
tryCatch({
  # Clean up speedup names if they have .elapsed suffix
  clean_names <- gsub("\\.elapsed$", "", names(speedups))
  
  speedup_df <- data.frame(
    Test = factor(clean_names, levels = clean_names),
    Speedup = as.numeric(speedups),
    Category = c("Small", "Medium", "Large", "Large (MDYPL)")
  )
  
  png(file.path(plots_dir, "03_speedup_summary.png"), 
      width = 800, height = 600, res = 100)
  
  p3 <- ggplot(speedup_df, aes(x = Test, y = Speedup, fill = Category)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = sprintf("%.2fx", Speedup)), 
              vjust = -0.5, size = 4) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 0.8) +
    labs(title = "Performance Improvement: New vs Original Implementation",
         subtitle = "Higher is better (baseline = 1.0)",
         x = "Test Case",
         y = "Speedup Factor") +
    scale_fill_brewer(palette = "Set2") +
    ylim(0, max(speedup_df$Speedup) * 1.15) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  
  print(p3)
  dev.off()
  cat("Plot 3 saved\n")
}, error = function(e) {
  cat("Error in Plot 3:", conditionMessage(e), "\n")
  dev.off()
})

# ====== Plot 4: Absolute timing comparison ======
cat("Creating Plot 4: Absolute timing...\n")
tryCatch({
  timing_df <- data.frame(
    Test = rep(c("Lizards", "Endometrial", "MultipleFeatures", "MDYPL"), each = 2),
    Version = rep(c("Original", "New"), 4),
    Time = c(
      median(lizards_bench$time[lizards_bench$expr == "original"]) / 1e9,
      median(lizards_bench$time[lizards_bench$expr == "new"]) / 1e9,
      median(endo_bench$time[endo_bench$expr == "original"]) / 1e9,
      median(endo_bench$time[endo_bench$expr == "new"]) / 1e9,
      as.numeric(time_mf_orig["elapsed"]),
      as.numeric(time_mf_new["elapsed"]),
      as.numeric(time_mdypl_orig["elapsed"]),
      as.numeric(time_mdypl_new["elapsed"])
    )
  )
  
  png(file.path(plots_dir, "04_absolute_timing.png"), 
      width = 800, height = 600, res = 100)
  
  p4 <- ggplot(timing_df, aes(x = Test, y = Time, fill = Version)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_text(aes(label = sprintf("%.3fs", Time)), 
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 3) +
    labs(title = "Absolute Execution Time Comparison",
         subtitle = "Lower is better",
         x = "Test Case",
         y = "Time (seconds)") +
    scale_fill_manual(values = c("Original" = "#E41A1C", "New" = "#4DAF4A")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  
  print(p4)
  dev.off()
  cat("Plot 4 saved\n")
}, error = function(e) {
  cat("Error in Plot 4:", conditionMessage(e), "\n")
  dev.off()
})

# ====== Plot 5: Distribution comparison (Lizards) ======
cat("Creating Plot 5: Lizards distribution...\n")
tryCatch({
  lizards_df <- data.frame(
    Time = lizards_bench$time / 1e9,
    Version = lizards_bench$expr
  )
  
  medians <- aggregate(Time ~ Version, lizards_df, median)
  
  png(file.path(plots_dir, "05_lizards_distribution.png"), 
      width = 800, height = 600, res = 100)
  
  p5 <- ggplot(lizards_df, aes(x = Time, fill = Version)) +
    geom_density(alpha = 0.6) +
    geom_vline(data = medians,
               aes(xintercept = Time, color = Version),
               linetype = "dashed", size = 1) +
    labs(title = "Execution Time Distribution: Lizards Dataset",
         subtitle = "Dashed lines show median values",
         x = "Time (seconds)",
         y = "Density") +
    scale_fill_manual(values = c("original" = "#E41A1C", "new" = "#4DAF4A")) +
    scale_color_manual(values = c("original" = "#E41A1C", "new" = "#4DAF4A")) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
  
  print(p5)
  dev.off()
  cat("Plot 5 saved\n")
}, error = function(e) {
  cat("Error in Plot 5:", conditionMessage(e), "\n")
  dev.off()
})

cat("\n=== Visualization Complete ===\n")
cat("All plots saved to:", plots_dir, "\n\n")
cat("Generated files:\n")
cat("  01_lizards_comparison.png - Boxplot comparison\n")
cat("  02_endometrial_comparison.png - Boxplot comparison\n")
cat("  03_speedup_summary.png - Speedup bar chart\n")
cat("  04_absolute_timing.png - Absolute time comparison\n")
cat("  05_lizards_distribution.png - Distribution densities\n")