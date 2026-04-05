Analysis of functional market similarity curves
===============================================

This code uses monthly returns of State Street Sector ETFs as a foundation for constructing a pairwise similarity structure across time.

## Sectors covered:
- XLY - Consumer Discretionary Sector
- XLP - Consumer Staples Sector
- XLE - Energy Sector
- XLF - Financial Sector
- XLV - Health Care Sector
- XLI - Industrial Sector
- XLB - Materials Sector
- XLK - Technology Sector
- XLU - Utilities Sector

## Time period:
- 2000 01 - 2026 02


## Interpretation:
Each row of the constructed matrix `K` is treated as a functional observation of a similarity curve for the given month showing how similar it is to all other months across the 26-year span.