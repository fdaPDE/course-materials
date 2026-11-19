# R notebooks

The five notebooks in `notebook/` follow the code, explanations, and order of the
Limerick R Markdown vignettes. They use the corresponding Lancaster datasets in
`../data/` and the copied plotting helpers in `utils/graphics.R`.

Open Jupyter with `R/notebook/` as the working directory and select the R
(`IRkernel`) kernel. The notebooks require `fdaPDE2` and the packages loaded by
`utils/graphics.R` (`dplyr`, `ggmap`, `ggplot2`, `latex2exp`, `leafsync`,
`mapview`, `raster`, `RTriangle`, `sf`, `viridis`, and `fields`). See the Limerick
course README for `fdaPDE2` installation instructions.

The Limerick data files are not stored in this repository. Their Lancaster
counterparts have different names and, in a few cases, different formats. The
notebooks build boundary segments from coordinates, combine the Switzerland
locations with log rainfall, and load Florida rasters from TIFF files. The
Florida transport figure is drawn directly from those rasters instead of a
pre-rendered map file.
