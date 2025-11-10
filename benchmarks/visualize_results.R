# Benchmark visualization script
# Run this after the main benchmark script completes

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

# ====== Plot 1: Microbenchmark comparison (Lizards) ======
png(file.path(plots_dir, "01_lizards_comparison.png"), 
    width = 800, height = 600, res = 100)
autoplot(lizards_bench) + 
  labs(title = "Lizards Dataset: Original vs New Implementation",
       subtitle = "100 evaluations each",
       y = "Time (seconds)") +
  theme_minimal(base_size = 12)
dev.off()

# ====== Plot 2: Microbenchmark comparison (Endometrial) ======
png(file.path(plots_dir, "02_endometrial_comparison.png"), 
    width = 800, height = 600, res = 100)
autoplot(endo_bench) + 
  labs(title = "Endometrial Dataset: Original vs New Implementation",
       subtitle = "50 evaluations each",
       y = "Time (seconds)") +
  theme_minimal(base_size = 12)
dev.off()

# ====== Plot 3: Speedup summary bar chart ======
speedup_df <- data.frame(
  Test = factor(names(speedups), levels = names(speedups)),
  Speedup = as.numeric(speedups),
  Category = c("Small", "Medium", "Large", "Large (MDYPL)")
)

png(file.path(plots_dir, "03_speedup_summary.png"), 
    width = 800, height = 600, res = 100)
ggplot(speedup_df, aes(x = Test, y = Speedup, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = sprintf("%.2fx", Speedup)), 
            vjust = -0.5, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(title = "Performance Improvement: New vs Original Implementation",
       subtitle = "Higher is better",
       x = "Test Case",
       y = "Speedup Factor") +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

# ====== Plot 4: Absolute timing comparison ======
timing_df <- data.frame(
  Test = rep(c("Lizards", "Endometrial", "MultipleFeatures", "MDYPL"), each = 2),
  Version = rep(c("Original", "New"), 4),
  Time = c(
    median(lizards_bench$time[lizards_bench$expr == "original"]) / 1e9,
    median(lizards_bench$time[lizards_bench$expr == "new"]) / 1e9,
    median(endo_bench$time[endo_bench$expr == "original"]) / 1e9,
    median(endo_bench$time[endo_bench$expr == "new"]) / 1e9,
    time_mf_orig["elapsed"],
    time_mf_new["elapsed"],
    time_mdypl_orig["elapsed"],
    time_mdypl_new["elapsed"]
  )
)

png(file.path(plots_dir, "04_absolute_timing.png"), 
    width = 800, height = 600, res = 100)
ggplot(timing_df, aes(x = Test, y = Time, fill = Version)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.2fs", Time)), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(title = "Absolute Execution Time Comparison",
       subtitle = "Lower is better",
       x = "Test Case",
       y = "Time (seconds)") +
  scale_fill_manual(values = c("Original" = "#E41A1C", "New" = "#4DAF4A")) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")
dev.off()

# ====== Plot 5: Distribution comparison (Lizards) ======
lizards_df <- data.frame(
  Time = lizards_bench$time / 1e9,  # Convert to seconds
  Version = lizards_bench$expr
)

png(file.path(plots_dir, "05_lizards_distribution.png"), 
    width = 800, height = 600, res = 100)
ggplot(lizards_df, aes(x = Time, fill = Version)) +
  geom_density(alpha = 0.6) +
  geom_vline(data = aggregate(Time ~ Version, lizards_df, median),
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
dev.off()

# ====== Summary statistics table ======
summary_stats <- rbind(
  data.frame(
    Dataset = "Lizards",
    Version = c("Original", "New"),
    Median = c(
      median(lizards_bench$time[lizards_bench$expr == "original"]) / 1e9,
      median(lizards_bench$time[lizards_bench$expr == "new"]) / 1e9
    ),
    Mean = c(
      mean(lizards_bench$time[lizards_bench$expr == "original"]) / 1e9,
      mean(lizards_bench$time[lizards_bench$expr == "new"]) / 1e9
    ),
    SD = c(
      sd(lizards_bench$time[lizards_bench$expr == "original"]) / 1e9,
      sd(lizards_bench$time[lizards_bench$expr == "new"]) / 1e9
    )
  ),
  data.frame(
    Dataset = "Endometrial",
    Version = c("Original", "New"),
    Median = c(
      median(endo_bench$time[endo_bench$expr == "original"]) / 1e9,
      median(endo_bench$time[endo_bench$expr == "new"]) / 1e9
    ),
    Mean = c(
      mean(endo_bench$time[endo_bench$expr == "original"]) / 1e9,
      mean(endo_bench$time[endo_bench$expr == "new"]) / 1e9
    ),
    SD = c(
      sd(endo_bench$time[endo_bench$expr == "original"]) / 1e9,
      sd(endo_bench$time[endo_bench$expr == "new"]) / 1e9
    )
  )
)

write.csv(summary_stats, 
          file.path(plots_dir, "summary_statistics.csv"),
          row.names = FALSE)

cat("\nAll plots saved to:", plots_dir, "\n")
cat("\nGenerated plots:\n")
cat("  1. Lizards microbenchmark comparison\n")
cat("  2. Endometrial microbenchmark comparison\n")
cat("  3. Overall speedup summary\n")
cat("  4. Absolute timing comparison\n")
cat("  5. Lizards timing distribution\n")
cat("  6. Summary statistics (CSV)\n")