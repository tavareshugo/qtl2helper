#' Extract markers from cross2 object
#'
#' This function gets both genetic and physical position of all markers in a `cross2`
#' object.
#'
#' @param x An object of class "cross2".
#'
#' @return a tibble
#'
#' @export
#' @rdname get_markers
get_markers <- function(x){
  UseMethod("get_markers")
}

#' @rdname get_markers
#' @export
get_markers.cross2 <- function(x){
  # get genetic map into a tibble
  gmap <- lapply(x$gmap, tibble::enframe, name = "marker", value = "genetic_pos")
  gmap <- dplyr::bind_rows(gmap, .id = "chrom")

  # get physical map into a tibble
  pmap <- lapply(x$pmap, tibble::enframe, name = "marker", value = "physical_pos")
  pmap <- dplyr::bind_rows(pmap, .id = "chrom")

  # return both joined together
  dplyr::full_join(gmap, pmap, by = c("marker", "chrom"))
}
