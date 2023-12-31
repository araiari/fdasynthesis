% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{compute_weights}
\alias{compute_weights}
\title{Compute the weights of each of the neighbors for weighted averaging}
\usage{
compute_weights(D_i, K, alpha_0 = 1, beta_0 = 1)
}
\arguments{
\item{D_i}{vector of length \eqn{N} : distances between one item and all the others}

\item{K}{scalar : number of neighbors to be considered.}

\item{alpha_0}{scalar, with \eqn{\alpha_0 >0} : concentation parameter of the
Dirichlet distribution from which the weights of neighbors are sampled.}

\item{beta_0}{scalar, with \eqn{\beta_0 > 0} : rate parameter
weighting the inverse distances in the weight computation}
}
\value{
: list of two elements:
\itemize{
\item \code{index_k} vector of length \code{K} : indices of the \code{K} closest neighbors
\item \code{p_k}  vector of length \code{K} : computed weights, one for each neighbors
}
}
\description{
This functions computes the weights to be assigned to each neighbors:
it looks at the distances of the k closest elements, sample weights
from a Dirichlet of parameters proportional to such distances, return
indices of the neighbors and weights.
}
