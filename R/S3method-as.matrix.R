#' Matrix coercion for qtl2 objects
#'
#' Convert the output of `qtl2` functions to `matrix`.
#'
#' @param x The object to coerce.
#' @rdname as.matrix.qtl2
#' @export
as.matrix.viterbi <- function(x){
  alleles <- attr(x, "alleles")

  # bind all matrices
  genos <- Reduce(cbind, x)

  # replace allele index with name
  out <- matrix(alleles[genos], ncol = ncol(genos))
  rownames(out) <- rownames(genos)
  colnames(out) <- colnames(genos)
  return(out)
}