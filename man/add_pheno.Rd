% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_phenotypes.R
\name{add_pheno}
\alias{add_pheno}
\alias{add_pheno.cross2}
\title{Add phenotype data to a cross2 object}
\usage{
add_pheno(x, pheno, idcol = 1L, retain_all = TRUE)

\method{add_pheno}{cross2}(x, pheno, idcol = 1L, retain_all = TRUE)
}
\arguments{
\item{x}{An object of class "cross2".}

\item{pheno}{A data.frame with phenotype data. All columns, apart from the one with
individual ID, will be converted to numeric values.}

\item{idcol}{A numeric or character value with the index or name of the column
that contains the individual IDs.}

\item{retain_all}{A logical value indicating whereas all individuals should be
retained in the cross2 object (default) or to only retain those that have
both genotype and phenotype data.}
}
\description{
Add phenotype data to a cross2 object
}
