# Install required packages
required_packages <- c("quantmod", "xts", "zoo", "dplyr", "tidyr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(quantmod)
library(xts)
library(zoo)
library(dplyr)
library(fda)

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

nbasis <- 20
rangeval <- c(min(time_points), max(time_points))

basis <- create.fourier.basis(rangeval = rangeval,
                              nbasis = nbasis)


# Smooth each sector return series
fd_obj <- smooth.basis(time_points,
                       K,
                       basis)$fd

plot(fd_obj, xlab="Time", ylab="Similarity",
     main="Smoothed Similarity Curves")


# Mean / Variance
plot(mean.fd(fd_obj))
plot(std.fd(fd_obj))


# FPCA
nharm = 4
pcalist = pca.fd(fd_obj, nharm, centerfns = TRUE)
plot(pcalist)