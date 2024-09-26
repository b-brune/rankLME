
# rankLME

The R-package `rankLME` implements a robust rank-based approach to the estimation of mixed effects models with random slopes. The methodology is presented in "A Rank-Based Estimation Method for Mixed Effects Models in the Presence of Outlying Data" by Barbara Brune, Irene Ortner and Peter Filzmoser (2024). 

The core function in the package is the function `ranklme`.

The package can be installed by running the following code:

```
if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("b-brune/rankLME")
```

In addition, the folder `inst` contains code to reproduce the data examples and the structure of the simulation studies presented in the paper.

