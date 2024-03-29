% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation_useful_functions.R
\name{plot_individual_changepoint_posterior}
\alias{plot_individual_changepoint_posterior}
\title{Plot empirical conditional posterior of the number of change-points.}
\usage{
plot_individual_changepoint_posterior(res, burn, pm, l = NULL)
}
\arguments{
\item{res}{the output obtained from the Metropolis-Hastings algorithms (either from infer_fixed_changepoints or infer_unknown_changepoints).}

\item{burn}{the proportion of the samples discarded as burn-in.}

\item{pm}{the desired range around the MAP location for each change-point location.}

\item{l}{condition on the number of change-points. If not initialised, the function expects as input the results obtained from the infer_fixed_changepoints function.}
}
\value{
plots of the empirical posterior distributions of the change-points given a specific number of change-points.
}
\description{
This function plots the conditional posterior distribution of the change-points locations given a specific number of change-points.
}
\examples{
# Use as an example the el_nino dataset.
# Run the function with l_max = 3 change-points, a maximum depth of 5 and the [0, 1] alphabet.
# The sampler is run for 10000 iterations.

res_unknown <- infer_unknown_changepoints(el_nino, 3, 5, c("01"), 100, fileName = NULL)

# Because l_max = 3 , there can be 0, 1, 2 or 3 changes.
# Let's see the posterior distribution on the number of changes

 plot_changepoint_posterior(res_unknown, 0.2)

# The MAP l is 2. Let's see the distribution of changes given l = 2.

plot_individual_changepoint_posterior(res_unknown, 0.2, 20, 2)

# One can also see the distribution of changes given l = 1. 

plot_individual_changepoint_posterior(res_unknown, 0.2, 500, 1)

# This function can be also used with the infer_fixed_changepoints
# Assume l = 2.

res_fixed <- infer_fixed_changepoints(el_nino, 2, 5, c("01"), 100, fileName = NULL)

# The function is now called without l = 2 as the number of changes is fixed 
# (all sampled vectors have 2 values). 

plot_individual_changepoint_posterior(res_fixed, 0.2, 20)
}
\seealso{
\code{\link{infer_fixed_changepoints}}, \code{\link{infer_unknown_changepoints}}
}
