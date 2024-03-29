% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{CTW}
\alias{CTW}
\title{Context Tree Weighting (CTW) algorithm}
\usage{
CTW(input_data, depth, desired_alphabet = NULL, beta = NULL)
}
\arguments{
\item{input_data}{the sequence to be analysed. 
The sequence needs to be a "character" object. See the examples section of the BCT/kBCT functions on how to transform any dataset to a "character" object.}

\item{depth}{maximum memory length.}

\item{desired_alphabet}{set containing the symbols of the process. If not initialised, the default set contains all the unique symbols which appear in the sequence. 
This parameter is needed for the segmentation problem where short segments might not contain all the symbols in the alphabet.}

\item{beta}{hyper-parameter of the model prior. 
Takes values between 0 and 1. If not initialised in the call function, the default value is \ifelse{html}{\out{1-2<sup>-m+1</sup>}}{\eqn{1 - 2^{-m+1}}}, 
where \ifelse{html}{\out{m}}{\eqn{m}} is the size of the alphabet; for more information see \href{https://arxiv.org/pdf/2007.14900.pdf}{Kontoyiannis et al. (2020)}.}
}
\value{
returns the natural logarithm of the prior predictive likelihood of the data.
}
\description{
Computes the prior predictive likelihood of the data given a specific alphabet. This function is used in for change-point point/segmentation problems
}
\examples{
# For the gene_s dataset with a maximum depth of 10 (with dafault value of beta):
CTW(gene_s, 10)

# With the ["0", "1", "2", "3"] alphabet
CTW(gene_s, 10, "0123")

# For custom beta (e.g. 0.8):
CTW(gene_s, 10, ,0.8)
}
\seealso{
\code{\link{BCT}}, \code{\link{kBCT}}
}
