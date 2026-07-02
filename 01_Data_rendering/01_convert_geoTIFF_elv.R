# install.packages("terra")  # run if not installed
library(terra)

in_tif <- "C:/Users/baldo/Desktop/EcoViz - Copy/EcoViz/data/BGD_KBT/dem_kost25x1.tif"
out_elv <- "C:/Users/baldo/Desktop/EcoViz - Copy/EcoViz/data/BGD_KBT/dem_kost25x1.elv"

r <- rast(in_tif)
m <- as.matrix(r, wide = TRUE)      # matrix: rows = raster rows (top->bottom)
nrows <- nrow(m)
ncols <- ncol(m)

originX <- xmin(r)                  # left-most x
originY <- ymax(r)                  # top-most y
cellsize_x <- res(r)[1]             # pixel width

# set nodata replacement if desired, e.g. nodata_value <- 600
nodata_value <- NA
if (!is.na(nodata_value)) m[is.na(m)] <- nodata_value

con <- file(out_elv, "w")
writeLines(sprintf("%d %d %.7f %.0f %.0f", ncols, nrows, originX, originY, cellsize_x), con)
for (ridx in 1:nrows) {
  rowvals <- m[ridx, ]
  if (!is.na(nodata_value)) rowvals[is.na(rowvals)] <- nodata_value
  writeLines(paste(as.integer(round(rowvals)), collapse = " "), con)
}
close(con)
