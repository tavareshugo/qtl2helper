#' Tidy R/qtl2 objects
#'
#' Convert the output of `qtl2` functions to a `tibble` in *long* format.
#' These are generally useful for downstream plotting analysis.
#'
#' @param x The object to coerce.
#' @param map An optional list of vectors of marker positions, as produced by
#' qtl2::insert_pseudomarkers() or the `$gmap` element of a `cross2` object.
#'
#' @details
#' At the moment the output from the following `R/qtl2` functions are recognised:
#' - `scan1()`
#' - `scan1coef()`
#'
#'
#' @rdname qtl2_tidiers
#' @export
tidy.scan1 <- function(x, map = NULL){

  # convert to tibble
  x_tbl <- tibble::as_tibble(x, rownames = "marker")

  # convert to long format
  x_tbl <- tidyr::gather(x_tbl, "pheno", "LOD", -marker)

  # add map data
  if(!is.null(map)){
    x_tbl <- .join_map_by_marker(x_tbl, map)
  }

  return(x_tbl)
}


#' @rdname qtl2_tidiers
#' @export
tidy.scan1coef <- function(x, map = NULL){

  # Convert coefficients to tibble
  coefs <- tibble::as_tibble(x, rownames = "marker")

  # reshape table to long format
  coefs <- tidyr::gather(coefs, "coef", "estimate", -marker)

  # Convert SEs to tibble and join to mean effect
  if("SE" %in% names(attributes(x))){
    SEs <- tibble::as_tibble(attr(x, "SE"), rownames = "marker")
    SEs <- tidyr::gather(SEs, "coef", "SE", -marker)

    coefs <- dplyr::full_join(coefs, SEs, by = c("marker", "coef"))
  }

  # add map data
  if(!is.null(map)){
    coefs <- .join_map_by_marker(coefs, map)
  }

  return(tibble::as_tibble(coefs))

}


#' @rdname qtl2_tidiers
#' @export
tidy.calc_genoprob <- function(x, map = NULL){

  # get the marker names as a vector
  markers <- lapply(x, function(i) attr(i, "dimnames")[[3]])
  markers <- unlist(markers)

  # create empty list for storing data
  marker_probs <- vector("list", length(markers))
  names(marker_probs) <- markers

  # pull each marker out and coerce to tibble
  for (chrom in names(x)){
    chrom_probs <- x[[chrom]]

    for (marker in attr(chrom_probs, "dimnames")[[3]]){
      marker_probs[[marker]] <- tibble::as_tibble(chrom_probs[, , marker])
    }
  }

  # bind them into a single tibble
  marker_probs <- dplyr::bind_rows(marker_probs, .id = "marker")

  # reshape to long format
  marker_probs <- tidyr::gather(marker_probs, "genotype", "probability", -marker)

  # remove missing values
  marker_probs <- tidyr::drop_na(marker_probs)

  # add map data
  if(!is.null(map)){
    marker_probs <- .join_map_by_marker(marker_probs, map)
  }

  # return the tibble
  return(marker_probs)

}


#' @param alpha Vector of significance levels (for `scan1perm` objects).
#'
#' @rdname qtl2_tidiers
#' @export
tidy.scan1perm <- function(x, alpha = 0.05){

  # get LOD thresholds
  perm_sum <- qtl2::summary_scan1perm(x, alpha)

  # convert to tibble
  if(is.list(perm_sum)){

    # if X chromosome was present
    perm_sum <- lapply(perm_sum, tibble::as_tibble, rownames = "alpha")
    perm_sum <- dplyr::bind_rows(perm_sum, .id = "chrom_type")

    # reshape to long format
    perm_sum <- tidyr::gather(perm_sum, "pheno", "threshold", -alpha, -chrom_type)

  } else {

    perm_sum <- tibble::as_tibble(perm_sum, rownames = "alpha")

    # reshape to long format
    perm_sum <- tidyr::gather(perm_sum, "pheno", "threshold", -alpha)

  }

  # return the tibble
  return(perm_sum)
}
