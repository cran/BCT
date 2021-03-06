% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ML}
\alias{ML}
\title{Maximum Likelihood}
\usage{
ML(input_data, depth)
}
\arguments{
\item{input_data}{the sequence to be analysed. 
The sequence needs to be a "character" object. See the examples section of the BCT/kBCT functions on how to transform any dataset to a "character" object.}

\item{depth}{maximum memory length.}
}
\value{
returns the natural logarithm of the maximum likelihood.
}
\description{
Computes the logarithm of the likelihood of the observations, maximised over all models and parameters.
}
\examples{
# Computing the maximum likelihood of the gene_s dataset 
# with a maximum depth of 5:
ML(gene_s, 5)
}
\seealso{
\code{\link{BCT}}, \code{\link{kBCT}}
}
