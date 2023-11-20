#' Crop and then transform
#'
#' Transforms bbox to raster CRS, then crops raster, then transforms raster to
#' given CRS
#'
#' @param ras SpatRaster
#' @param bbx sf polygon of bbox
#' @param crs_out crs to use. As returned by `terra::crs`
#'
#' @return cropped and transformed SpatRaster
crop_transform <- function(ras, bbx, crs_out){
  bbx1 <- st_transform(bbx, st_crs(ras)) %>% terra::vect()

  bbx2 <- st_transform(bbx, crs_out) %>% terra::vect()

  ras <- terra::crop(ras, bbx1, snap = "out")

  ras <- terra::project(ras, crs_out)

  ras <- terra::crop(ras, bbx2, snap = "out")

  return(ras)
}


extract_transform <- function(ras, pts){
  pts <- terra::project(pts, terra::crs(ras))

  terra::extract(ras, pts)
}
