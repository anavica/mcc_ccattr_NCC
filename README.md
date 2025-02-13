# The burden of heat-related mortality attributable to recent human-induced climate change 

Data and `R` code to support:

Vicedo-Cabrera, Scovronick N, Sera F,..., Gasparrini A. The burden of heat-related mortality attributable to recent human-induced climate change. *Nature Climate Change*. 2021;11(6):492-500 ([DOI: 10.1038/s41558-021-01058-x](https://doi.org/10.1038/s41558-021-01058-x)). 

The work was supported by the Medical Research Council-UK (Grant ID: MR/M022625/1), the Natural Environment Research Council UK (Grant ID: NE/R009384/1), and the European Union’s Horizon 2020 Project Exhaustion (Grant ID: 820655).

This code has been implemented to reproduce the main analysis applied in the
study referenced above. It illustrates a simplified example using temperature-mortality data of 10 hypothetical
cities and the corresponding simulated temperature for 10 climate models defining the
two scenarios.

The analysis is divided in three stages: (1) 1st stage analysis (time-series
analysis), (2) Meta-regression (to extract the BLUPs and MMTs),
(3) Quantification of impacts under the two scenarios (factual / counterfactual)
