#' Add phenotype data to a cross2 object
#'
#' @param x An object of class "cross2".
#' @param pheno A data.frame with phenotype data. All columns, apart from the one with
#' individual ID, will be converted to numeric values.
#' @param idcol A numeric or character value with the index or name of the column
#' that contains the individual IDs.
#' @param retain_all A logical value indicating whereas all individuals should be
#' retained in the cross2 object (default) or to only retain those that have
#' both genotype and phenotype data.
#'
#' @export
#' @rdname add_pheno
add_pheno <- function(x, pheno, idcol = 1L, retain_all = TRUE){
  UseMethod("add_pheno")
}

#' @export
#' @rdname add_pheno
add_pheno.cross2 <- function(x, pheno, idcol = 1L, retain_all = TRUE){

  #### Format pheno data.frame ####

  # Coerce to data.frame (in case tibble is provided)
  pheno <- as.data.frame(pheno)

  # Get column name, in case it is provided as a character
  if(is.character(idcol)){
    idcol <- which(colnames(pheno) %in% idcol)
    if(length(idcol) != 1) stop("idcol does not seem to be unique!")
  }

  # Get vector of individual IDs in phenotype provided
  ids <- as.character(pheno[[idcol]])

  # Add these to rownames of phenotype table
  rownames(pheno) <- ids

  # If user wants to remove non-genotyped individuals
  if(!retain_all){
    # Retain only those individuals that are in the cross object
    ids <- ids[which(ids %in% rownames(x$geno[[1]]))]

    # Subset phenotype table to retain only those individuals
    pheno <- pheno[ids, -idcol]
  }

  # remove ID column from pheno table
  pheno <- pheno[, -idcol]

  # Convert to matrix using qtl2 helper function
  pheno <- qtl2:::pheno2matrix(pheno)


  #### Add to cross2 object ####

  if(!retain_all){
    # Subset cross2 object to retain only IDs in phenotype file
    x <- x[rownames(pheno), ]
  }

  # Add pheno slot
  x[["pheno"]] <- pheno

  return(x)
}

