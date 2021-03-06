% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Function_qqplot.R
\name{qq_plot}
\alias{qq_plot}
\title{Function to create a QQ plot}
\usage{
qq_plot(GM, pvals, QTN_index = c(), trait = "unknown")
}
\arguments{
\item{GM}{genetic map of data with chr and position of each SNP}

\item{pvals}{pvals from gwas results for each SNP}

\item{QTN_index}{posistion of QTN if applicable}

\item{trait}{character value for trait name}
}
\value{
QQ plot
}
\description{
Function to create a QQ plot
}
