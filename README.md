# Analysis of functional market similarity curves

This project constructs a pairwise similarity structure across time using monthly returns of State Street Sector ETFs, then applies functional data analysis (FDA) techniques to explore structural patterns, identify regime changes, and test whether market similarity behaves differently during recessions versus expansions.

## Sectors Covered

| Ticker | Sector |
|--------|--------|
| XLY | Consumer Discretionary |
| XLP | Consumer Staples |
| XLE | Energy |
| XLF | Financial |
| XLV | Health Care |
| XLI | Industrial |
| XLB | Materials |
| XLK | Technology |
| XLU | Utilities |

## Time Period

2000-01 through 2026-02 (monthly frequency)

## Methodology

### 1. Data Acquisition & Returns

Monthly discrete returns are computed from adjusted closing prices fetched via `quantmod` from Yahoo Finance (source range: 1999-12-31 to 2026-03-01). The nine sector ETF return series are merged and NA-removed into a single matrix **X** of dimension *n × 9*, where *n* is the number of months.

### 2. Similarity Matrix Construction

A Mahalanobis distance matrix **D** is computed over all pairs of months:

$$D[i,j] = (X[i] − X[j])^T Σ^{-1} (X[i] − X[j])$$

where Σ is the sample covariance matrix of returns. An RBF (Gaussian) kernel is then applied:

$$K[i,j] = exp(−γ · D[i,j])$$

with $γ = 1 / IQR(D)$

Each column of $K$ is interpreted as a functional observation: a similarity curve for a given month, describing how similar that month is to every other month across the 26-year span.

### 3. Functional Smoothing (GCV-Optimized)

All similarity curves are smoothed using B-spline basis expansions:

- **Basis**: 50 B-spline basis functions, order 4
- **Roughness penalty**: Second derivative (Lfdobj = 2)
- **Smoothing parameter selection**: Grid search over λ ∈ {10⁻¹⁰, …, 10¹⁰} (step 0.5), minimizing mean Generalized Cross-Validation (GCV) across all curves

The optimal λ is selected and applied uniformly to all curves.

### 4. Descriptive Functional Statistics

- **Mean function** and **pointwise ±1 SD bands** via `mean.fd()` and `std.fd()`
- **Covariance surface** via `var.fd()`, visualized as perspective and contour plots
- **Correlation surface** derived by normalizing the covariance surface by the outer product of standard deviations

### 5. Functional Principal Component Analysis (FPCA)

- 4 principal components extracted via `pca.fd()`
- **Changepoint detection** on PC1 scores using the `changepoint` package:
  - Method: SegNeigh (segment neighborhood)
  - Penalty: BIC
  - Maximum changepoints (Q): 7
  - Detected changepoint dates are mapped back to the original time index

### 6. Hypothesis Testing: Recession vs. Expansion

**NBER recession periods used:**

| Start | End |
|-------|-----|
| 2001-03 | 2001-11 |
| 2007-12 | 2009-06 |
| 2020-02 | 2020-04 |

Curves (columns of **K**) are split into recession and expansion groups, then re-smoothed using the **same** basis and optimal λ to ensure fair comparison.

**Test statistic:** L² distance between group mean curves, approximated via numerical integration on a 100-point grid.

Two permutation tests are performed (1000 permutations each):

| Test | Description |
|------|-------------|
| **Standard permutation** | Randomly reassigns curves to groups, ignoring temporal dependence |
| **Block permutation** | Shuffles non-overlapping blocks of various size to preserve local autocorrelation structure |

A significant result under standard permutation and block permutation suggests significant recession effects.