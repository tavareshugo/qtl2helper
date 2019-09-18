#' Utility functions for qtl2helper
#'
#' These functions are not exported.
#'
#' @rdname qtl2helper_utils


#' @param .tbl a `data.frame` or `tibble` to add the map data to.
#' Should contain a column named "marker".
#' @param .map the marker map, for example from `qtl2::insert_pseudomarkers()`.
#'
#' @rdname qtl2helper_utils
.join_map_by_marker <- function(.tbl, .map){

  # check inputs
  if(!is.list(.map)) stop(".map should be a list as produced by qtl2::insert_pseudomarkers().")

  # convert to tibble
  map_tbl <- lapply(.map, tibble::enframe, name = "marker", value = "pos")

  # bind chromosomes together
  map_tbl <- dplyr::bind_rows(map_tbl, .id = "chrom")

  # reorder columns
  map_tbl <- dplyr::select(map_tbl, marker, chrom, pos)

  # join with table
  .tbl <- dplyr::right_join(map_tbl, .tbl, by = "marker")

  # coerce chromosome to factor (to ensure ordering when plotting)
  .tbl$chrom <- factor(.tbl$chrom, levels = names(.map))

  # check all markers are present
  if(any(!(.tbl$marker %in% map_tbl$marker))){
    warning("Some markers were not included in the output!")
  }

  return(.tbl)
}

