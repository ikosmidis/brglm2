# Load required packages
library(brglm2)
library(microbenchmark)

# Load the dataset
data("MultipleFeatures", package = "brglm2")

# Set up the data and formulas (based on the package documentation)
# Get the fou (Fourier) and kar (Karhunen-Loève) feature names
vars <- grep("fou|kar", names(MultipleFeatures), value = TRUE)

# Identify training rows
training <- which(MultipleFeatures$training)

# Center the features on the training set
MultipleFeatures[training, vars] <- scale(MultipleFeatures[training, vars], 
                                           scale = FALSE)

# Create the full model formula (predicting digit 7)
full_mf_fm <- formula(paste("I(digit == 7) ~", paste(vars, collapse = " + ")))

# ===== PROFILING STARTS HERE =====

# 1. Detailed profiling with Rprof 
print("\n === Rprof Function Timing (ordered by time) ===")
Rprof("brglm_profile.out", interval = 0.01, memory.profiling = TRUE)
fit <- glm(full_mf_fm, 
           data = MultipleFeatures, 
           family = binomial(),
           method = brglmFit,
           subset = training, 
           maxit = 200)
Rprof(NULL)

# Get summary showing functions by total time
prof_summary <- summaryRprof("brglm_profile.out", memory = "both")
print("Top functions by total time:")
print(head(prof_summary$by.total, 20))
print("\nTop functions by self time (excluding called functions):")
print(head(prof_summary$by.self, 20))

# 2. Visual profiling with profvis
print("\n === Memory/Time Profiling with profvis ===")
library(profvis)

p <- profvis({
  fit <- glm(full_mf_fm, 
             data = MultipleFeatures, 
             family = binomial(),
             method = brglmFit,
             subset = training, 
             maxit = 200)
}, interval = 0.01)  # Smaller interval for more detail

print(p)