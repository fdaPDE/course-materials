<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_dark.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_light.png">
  <img alt="fdaPDE logo" src="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_light.png">
</picture>
</p>

Material for the short course “Physics-Informed Statistical Learning for Spatial and Functional Data” held as part of the STOR-i masterclass program of the Lancaster University, Lancaster (UK), 24-25 November 2026.

__Lectures: Laura Maria Sangalli <sup>1</sup>, Alessandro Palummo <sup>1</sup>__.

&nbsp;&nbsp;&nbsp;&nbsp;_<sup>1</sup> MOX, Department of Mathematics, Politecnico di Milano_

# Getting started

`fdaPDE` is a C++ library interfacing with R, one of the most widely used langauges for data analysis. The version of `fdaPDE` used in this course is not yet available as an official CRAN package, and a development version will be considered here.

## Colab instructions

To run the notebooks on [Google Colab](https://colab.google/), just open one of the link below. If necessary, grant permissions to Google to trust the notebook.

**Be sure to run the setup cell at the beginning of each notebook to download the data and install all the required dependencies**.

* [Spatial regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/SRPDE_2D.ipynb)
* [Generalized regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/GSRPDE_2D.ipynb)
* [Quantile regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/QSRPDE_2D.ipynb)
* [Density estimation](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/DEPDE_2D.ipynb)
* [Functional Principal Component Analysis](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Lancaster]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/notebook/fPCA_2D.ipynb)

## Installation instructions

`fdaPDE` requires Rcpp, RcppEigen, R6 and a C++20-compatible compiler to be installed (for instance `gcc14.2`). The notebook contained in this respository further require additional graphical dependencies. Install R and Jupyter Notebook for your operating system, then execute from the R console

```r
install.packages(c(
  "IRkernel", "devtools", "dplyr", "fields", "ggmap", "ggplot2",
  "latex2exp", "leafsync", "mapview", "raster", "RTriangle", "sf",
  "viridis", "Rcpp", "RcppEigen", "R6"
))

devtools::install_github("fdaPDE/fdaPDE-R")
IRkernel::installspec(user = TRUE, name = "ir", displayname = "R")
```

Some graphical dependencies require GDAL. If the installation of the requirements fails, make sure that GDAL is installed. 

To run the notebooks, open a terminal in the `courses/\[2026-Lancaster\]-Physics-Informed-Statistical-Learning-for-Spatial-and-Functional-Data/R/` folder and start Jupyter:

```sh
jupyter notebook
```

If your browser does not open Jupyter automatically, inspect the output produced by the command above and copy-paste the link `http://localhost:8888/tree?token=...` into your browser.

Select a notebook from the `notebook/` folder, then select the R kernel by clicking on "Kernel" → "Change kernel..." and then selecting "R".

You are now up and running. Enjoy the course!
