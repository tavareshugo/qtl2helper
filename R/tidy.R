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
  x_tbl <- as.data.frame(x) # to remove `scan1` class
  x_tbl <- tibble::as_tibble(x_tbl, rownames = "marker")

  # convert to long format
  x_tbl <- tidyr::pivot_longer(x_tbl, 
                               cols = c(-marker),
                               names_to = "pheno", 
                               values_to = "LOD")

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
  coefs <- as.data.frame(x) # to remove `scan1coef` class
  coefs <- tibble::as_tibble(coefs, rownames = "marker")

  # reshape table to long format
  coefs <- tidyr::pivot_longer(coefs, 
                               cols = c(-marker),
                               names_to = "coef", 
                               values_to = "estimate")

  # Convert SEs to tibble and join to mean effect
  if("SE" %in% names(attributes(x))){
    SEs <- tibble::as_tibble(attr(x, "SE"), rownames = "marker")
    SEs <- tidyr::pivot_longer(SEs, 
                               cols = c(-marker),
                               names_to = "coef", 
                               values_to = "SE")

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
      marker_probs[[marker]] <- tibble::as_tibble(chrom_probs[, , marker],
                                                  rownames = "id")
      # marker_probs[[marker]] <- as.data.frame(chrom_probs[, , marker])
    }
  }

  # bind them into a single tibble
  marker_probs <- dplyr::bind_rows(marker_probs, .id = "marker")
  # marker_probs <- do.call(rbind, marker_probs)

  # reshape to long format
  marker_probs <- tidyr::pivot_longer(marker_probs, 
                                      cols = c(-marker, -id),
                                      names_to = "genotype", 
                                      values_to = "probability")

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
    perm_sum <- lapply(perm_sum, function(i){
      tibble::as_tibble(as.data.frame(i), rownames = "alpha")
    })
    perm_sum <- dplyr::bind_rows(perm_sum, .id = "chrom_type")

    # reshape to long format
    perm_sum <- tidyr::pivot_longer(perm_sum, 
                                    cols = c(-alpha, -chrom_type),
                                    names_to = "pheno", 
                                    values_to = "threshold")

  } else {

    perm_sum <- tibble::as_tibble(as.data.frame(perm_sum), rownames = "alpha")

    # reshape to long format
    perm_sum <- tidyr::pivot_longer(perm_sum, 
                              cols = c(-alpha),
                              names_to = "pheno", 
                              values_to = "threshold")

  }

  # return the tibble
  return(perm_sum)
}
