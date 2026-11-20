<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_dark.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_light.png">
  <img alt="fdaPDE logo" src="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_light.png">
</picture>
</p>

Material for the short course “Physics-Informed Statistical Learning for Spatial and Functional Data”.

__Lectures: Laura Maria Sangalli <sup>1</sup>, Alessandro Palummo <sup>1</sup>__.

&nbsp;&nbsp;&nbsp;&nbsp;_<sup>1</sup> MOX, Department of Mathematics, Politecnico di Milano_

# Getting started

The notebooks use the R interface to `fdaPDE` for spatial regression, generalized regression, quantile regression, density estimation, and functional principal component analysis.

## Colab instructions

Open a notebook using one of the links below. To run the R notebooks locally, follow the installation instructions in the next section.

* [Spatial regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/SRPDE_2D.ipynb)
* [Generalized regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/GSRPDE_2D.ipynb)
* [Quantile regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/QSRPDE_2D.ipynb)
* [Density estimation](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/DEPDE_2D.ipynb)
* [Functional Principal Component Analysis](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/fPCA_2D.ipynb)

Google Colab does not officially support R kernels. To execute these notebooks, use a local Jupyter installation with R support.

## Installation instructions

Install R and Jupyter Notebook for your operating system. The following R commands are the same on Linux, macOS, and Windows. Run them in the R console:

```r
install.packages(c(
  "IRkernel", "devtools", "dplyr", "fields", "ggmap", "ggplot2",
  "latex2exp", "leafsync", "mapview", "raster", "RTriangle", "sf",
  "viridis", "Rcpp", "RcppEigen", "R6"
))

devtools::install_github("fdaPDE/fdaPDE-R")
IRkernel::installspec(user = TRUE, name = "ir", displayname = "R")
```

To run the notebooks, open a terminal in the `R/notebook/` folder and start Jupyter:

```sh
jupyter notebook
```

Select the **R** kernel when you open a notebook.
