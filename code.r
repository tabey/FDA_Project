# Install required packages
required_packages <- c("quantmod", "xts", "zoo", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(quantmod)
library(xts)
library(zoo)
library(dplyr)
library(fda)
library(robust)

# Function to get monthly returns
get_monthly_returns <- function(ticker) {
  data <- getSymbols(ticker, src = "yahoo",
                     from = "1999-12-31",
                     to = "2026-03-01",
                     auto.assign = FALSE)
  
  close_prices <- Cl(data)
  
  monthly_prices <- to.monthly(close_prices, indexAt = "lastof", OHLC = FALSE)
  
  monthly_returns <- ROC(monthly_prices, type = "discrete")
  
  return(monthly_returns)
}

tickers <- c("XLY","XLP","XLE","XLF","XLV","XLI","XLB","XLK","XLU")

# Fetch monthly return data
monthly_returns_list <- lapply(tickers, get_monthly_returns)

# Combine into one dataframe
combined_monthly_returns <- do.call(merge, monthly_returns_list)
colnames(combined_monthly_returns) <- tickers

combined_monthly_returns <- na.omit(combined_monthly_returns)

# Format dates
dates <- as.Date(format(index(combined_monthly_returns), "%Y-%m-%d"))

# Convert to matrix
X <- coredata(combined_monthly_returns)

# Covariance matrix
cov_mat <- cov(X)
# cov_mat <- covRob(X)$cov
cov_inv <- solve(cov_mat)

n <- nrow(X)

# Mahalanobis distance matrix
D <- matrix(0, n, n)

for (i in 1:n) {
  for (j in 1:n) {
    diff <- X[i,] - X[j,]
    D[i,j] <- t(diff) %*% cov_inv %*% diff
  }
}

# IQR
gamma <- 1 / IQR(as.vector(D))

# RBF Kernel
K <- exp(-gamma * D)

# Convert to dataframe
S <- data.frame(K)
rownames(S) <- dates
colnames(S) <- dates

# Export CSV
# write.csv(S, "sim.csv")

# Time index
time_points <- 1:nrow(K)

# Create basis
nbasis <- 50
rangeval <- c(min(time_points), max(time_points))
basis <- create.bspline.basis(rangeval = rangeval,
                              nbasis = nbasis,
                              norder = 4)

plot(basis)

# Grid search over lambda values
lambda_vec <- 10^seq(-10, 10, by = 0.5)
gcv_mean <- numeric(length(lambda_vec))

for (i in seq_along(lambda_vec)) {
  fdPar_i <- fdPar(basis, Lfdobj = 2, lambda = lambda_vec[i])
  smooth_i <- smooth.basis(time_points, K, fdPar_i)
  gcv_mean[i] <- mean(smooth_i$gcv)  # average GCV across all curves
}
smooth_i$gcv
# Select optimal lambda
lambda_opt <- lambda_vec[which.min(gcv_mean)]

# Plot GCV curve
plot(log10(lambda_vec), gcv_mean, type = "b",
     xlab = "log10(lambda)", ylab = "Mean GCV",
     main = "GCV Optimization")
abline(v = log10(lambda_opt), col = "red", lty = 2)
lambda_opt

# Final smooth with optimal lambda
fdPar_opt <- fdPar(basis, Lfdobj = 2, lambda = lambda_opt)
smoothb <- smooth.basis(time_points, K, fdPar_opt)
fd_obj <- smoothb$fd

# Individual curves
plotfit.fd(K, time_points, smooth.basis(time_points, K, fdPar_opt)$fd, type='l')

plot(fd_obj, lwd=0.5, xlab = "Time", ylab = "Similarity",
     main = "Smoothed Similarity Curves (GCV-optimized)")

# Mean / Variance
meanfd  = mean.fd(fd_obj)
stdfd = std.fd(fd_obj)

lines(meanfd, lwd=4, lty=1, col=1)
lines(meanfd-stdfd, lwd=4, lty=2, col=1)
lines(meanfd+stdfd, lwd=4, lty=2, col=1)

# Covariance Surface
time        = seq(1,314,length=100)
logprecvar.bifd = var.fd(fd_obj)

logprecvar_mat  = eval.bifd(time, time,
                            logprecvar.bifd)

persp(time, time, logprecvar_mat,
      theta=-45, phi=25, r=3, expand = 0.5,
      ticktype='detailed',
      xlab="Month",
      ylab="Month",
      zlab="Variance")

contour(time, time, logprecvar_mat,
        col=terrain.colors(12),
        xlab="Month",
        ylab="Month",
        lwd=2,
        labcex=1)

# The diagonal of the covariance matrix
variances <- diag(logprecvar_mat)

# Standard deviations
std_devs <- sqrt(variances)

# Outer product
outer_std <- outer(std_devs, std_devs)

# Correlation surface
correlation_mat <- logprecvar_mat / outer_std

persp(time, time, correlation_mat,
      theta=-45, phi=25, r=3, expand = 0.5,
      ticktype='detailed',
      xlab="Month",
      ylab="Month",
      zlab="Correlation")

contour(time, time, correlation_mat,
        col=terrain.colors(12),
        xlab="Month",
        ylab="Month",
        lwd=2,
        labcex=1)

# Maybe outliers here later
library(fdaoutlier)
tt <- seq(1,314,length=314)
lgp <- eval.fd(tt,fd_obj)

fbplot(lgp, method="BD2")
fbplot(lgp, method="MBD")
fbplot(lgp, method="Both")

fboxplot(data = smoothb, plot.type = "bivariate", ylim=c(-3,3), xlim=c(-5,8),
         type = "bag", projmethod="PCAproj")


foutliers(smoothb, method = "robMah")
# doesn't work?
# foutliers(smoothb, method = "lrt")
# foutliers(smoothb, method = "depth.trim")
# foutliers(smoothb, method = "depth.pond")
# foutliers(smoothb, method = "HUoutliers")

# FPCA
par(mfrow=c(2,2))
nharm = 4
pcalist = pca.fd(fd_obj, nharm, centerfns = TRUE)
plot(pcalist)

par(mfrow=c(1,1))
plot(pcalist$harmonics, lwd=2)

# Rotation (useless?)
# varmx <- varmx.pca.fd(pcalist)
# plot(varmx)
# plot(varmx$harmonics)

# --- HYPOTHESIS TESTING: Recession vs. Expansion ---

# 2. Define NBER Recession Periods
nber_recessions <- list(
  c(as.Date("2001-03-01"), as.Date("2001-11-01")),
  c(as.Date("2007-12-01"), as.Date("2009-06-01")),
  c(as.Date("2020-02-01"), as.Date("2020-04-01"))
)

is_recession <- function(date) {
  for (period in nber_recessions) {
    if (date >= period[1] && date <= period[2]) {
      return(TRUE)
    }
  }
  return(FALSE)
}

group_labels <- sapply(dates, is_recession)

# 3. Split CURVES (columns), not rows
# Each column of K is a curve representing one month's similarity profile
recession_col_idx <- which(group_labels)
expansion_col_idx <- which(!group_labels)

cat("Recession months:", length(recession_col_idx), "\n")
cat("Expansion months:", length(expansion_col_idx), "\n")

if (length(recession_col_idx) < 5 || length(expansion_col_idx) < 5) {
  warning("Not enough data points in one group.")
} else {
  # Extract curves (columns) for each group
  K_recession <- K[, recession_col_idx, drop = FALSE]  # Keep all rows, select columns
  K_expansion <- K[, expansion_col_idx, drop = FALSE]
  
  # Verify dimensions match time_points
  cat("K_recession dimensions:", dim(K_recession), "\n")
  cat("time_points length:", length(time_points), "\n")
  
  # Re-smooth using SAME basis and lambda (critical for fair comparison)
  fdPar_rec <- fdPar(basis, Lfdobj = 2, lambda = lambda_opt)
  smooth_rec <- smooth.basis(time_points, K_recession, fdPar_rec)
  fd_rec <- smooth_rec$fd
  
  fdPar_exp <- fdPar(basis, Lfdobj = 2, lambda = lambda_opt)
  smooth_exp <- smooth.basis(time_points, K_expansion, fdPar_exp)
  fd_exp <- smooth_exp$fd
  
  # 4. Calculate Test Statistic
  eval_grid <- seq(1, 314, length.out = 100)
  
  mean_rec <- eval.fd(eval_grid, mean.fd(fd_rec))
  mean_exp <- eval.fd(eval_grid, mean.fd(fd_exp))
  
  diff_curve <- mean_rec - mean_exp
  step_size <- eval_grid[2] - eval_grid[1]
  stat_obs <- sum(diff_curve^2) * step_size
  
  # 5. Permutation Test
  n_perm <- 1000
  stat_perm <- numeric(n_perm)
  
  # Combine all curves (columns)
  all_curves <- cbind(K_recession, K_expansion)
  n_total <- ncol(all_curves)
  n_rec <- ncol(K_recession)
  
  set.seed(123)
  
  cat("Running Permutation Test...\n")
  for (i in 1:n_perm) {
    perm_indices <- sample(n_total)
    perm_rec <- all_curves[, perm_indices[1:n_rec], drop = FALSE]
    perm_exp <- all_curves[, perm_indices[(n_rec+1):n_total], drop = FALSE]
    
    # Re-smooth permuted groups
    smooth_perm_rec <- smooth.basis(time_points, perm_rec, fdPar_rec)
    smooth_perm_exp <- smooth.basis(time_points, perm_exp, fdPar_exp)
    
    mean_perm_rec <- eval.fd(eval_grid, mean.fd(smooth_perm_rec$fd))
    mean_perm_exp <- eval.fd(eval_grid, mean.fd(smooth_perm_exp$fd))
    
    diff_perm <- mean_perm_rec - mean_perm_exp
    stat_perm[i] <- sum(diff_perm^2) * step_size
  }
  
  p_value <- (sum(stat_perm >= stat_obs) + 1) / (n_perm + 1)
  
  cat("\n--- Hypothesis Test Results ---\n")
  cat("Observed Statistic:", stat_obs, "\n")
  cat("P-value:", p_value, "\n")
  
  if (p_value < 0.05) {
    cat("Result: Significant difference (p < 0.05)\n")
  } else {
    cat("Result: No significant difference\n")
  }
  
  # Plot
  plot(eval_grid, mean_exp, type = "l", col = "blue", lwd = 2,
       xlab = "Time Index", ylab = "Mean Similarity",
       main = "Mean Similarity Curves: Recession vs Expansion")
  lines(eval_grid, mean_rec, col = "red", lwd = 2)
  legend("topright", legend = c("Expansion", "Recession"), col = c("blue", "red"), lwd = 2)
}

# --- BLOCK PERMUTATION TEST ---

# Parameters
block_size <- 12  # Try 3, 6, or 12 months. 6 is a good starting point for monthly data.
n_blocks <- ceiling(n_total / block_size) # Number of blocks needed

# Ensure we cover the whole dataset (circular padding if necessary)
# We will create a circular buffer to handle the wrap-around
K_circular <- cbind(K, K[, 1:(block_size-1)]) # Append first few cols to end

# Function to get a single block permutation
get_block_permutation_stat <- function(K_data, block_size, n_rec, time_points, basis, lambda_opt, eval_grid) {
  n_cols <- ncol(K_data)
  n_blocks <- ceiling(n_cols / block_size)
  
  # Create indices for blocks (start position of each block)
  # We use a circular approach: if a block goes past the end, it wraps to the start
  block_starts <- seq(1, n_cols, by = block_size)
  
  # If the last block is incomplete, we still treat it as a block
  # But for shuffling, we need to handle the wrap-around carefully.
  # Simpler approach: Just shuffle the blocks of the original data, ignoring the circular padding for the shuffle logic,
  # but handling the wrap-around when extracting.
  
  # Let's use a simpler "Non-overlapping" block shuffle for clarity, 
  # but pad the end to make it divisible if needed.
  
  # Pad K to be divisible by block_size
  remainder <- n_cols %% block_size
  if (remainder != 0) {
    pad_cols <- block_size - remainder
    K_padded <- cbind(K_data, K_data[, 1:pad_cols]) # Wrap around padding
  } else {
    K_padded <- K_data
  }
  
  n_padded <- ncol(K_padded)
  n_blocks_padded <- n_padded / block_size
  
  # Create block indices
  block_indices <- matrix(1:n_padded, nrow = block_size, ncol = n_blocks_padded)
  # Transpose so each column is a block
  block_indices <- t(block_indices) 
  
  # Shuffle the blocks (rows of block_indices)
  shuffled_block_rows <- sample(nrow(block_indices))
  shuffled_blocks <- block_indices[shuffled_block_rows, ]
  
  # Reconstruct the matrix
  # Flatten the shuffled blocks back into a vector of column indices
  new_order <- as.vector(shuffled_blocks)
  
  # Extract the permuted matrix (only the first n_cols columns to remove padding)
  K_perm <- K_padded[, new_order[1:n_cols], drop = FALSE]
  
  # Now assign labels: First n_rec columns are "Recession", rest "Expansion"
  # (This mimics the original logic: we are testing if the FIRST n_rec blocks 
  # happen to be different from the rest, under random block ordering)
  
  K_perm_rec <- K_perm[, 1:n_rec, drop = FALSE]
  K_perm_exp <- K_perm[, (n_rec+1):ncol(K_perm), drop = FALSE]
  
  # Re-smooth
  smooth_perm_rec <- smooth.basis(time_points, K_perm_rec, fdPar(basis, Lfdobj = 2, lambda = lambda_opt))
  smooth_perm_exp <- smooth.basis(time_points, K_perm_exp, fdPar(basis, Lfdobj = 2, lambda = lambda_opt))
  
  # Calculate Statistic
  mean_perm_rec <- eval.fd(eval_grid, mean.fd(smooth_perm_rec$fd))
  mean_perm_exp <- eval.fd(eval_grid, mean.fd(smooth_perm_exp$fd))
  diff_perm <- mean_perm_rec - mean_perm_exp
  stat_perm_val <- sum(diff_perm^2) * step_size
  
  return(stat_perm_val)
}

# Run the Block Permutation Test
n_perm_block <- 1000
stat_perm_block <- numeric(n_perm_block)

cat("Running BLOCK Permutation Test (Block Size =", block_size, ")...\n")

for (i in 1:n_perm_block) {
  stat_perm_block[i] <- get_block_permutation_stat(
    K, block_size, n_rec, time_points, basis, lambda_opt, eval_grid
  )
}

# Calculate P-value
p_value_block <- (sum(stat_perm_block >= stat_obs) + 1) / (n_perm_block + 1)

cat("\n--- Block Permutation Test Results ---\n")
cat("Block Size:", block_size, "months\n")
cat("Observed Statistic:", stat_obs, "\n")
cat("Block P-value:", p_value_block, "\n")

if (p_value_block < 0.05) {
  cat("Result: Significant difference even with block permutation (Robust).\n")
} else {
  cat("Result: Not significant with block permutation (Dependence may be driving the result).\n")
}

# Optional: Compare distributions
par(mfrow=c(1,2))
hist(stat_perm, main="Standard Permutation", col="lightblue", border="white")
abline(v=stat_obs, col="red", lwd=2)
hist(stat_perm_block, main=paste("Block Permutation (L=", block_size, ")", sep=""), col="lightgreen", border="white")
abline(v=stat_obs, col="red", lwd=2)
