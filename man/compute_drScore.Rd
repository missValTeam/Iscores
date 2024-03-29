% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{compute_drScore}
\alias{compute_drScore}
\title{compute the density ratio score}
\usage{
compute_drScore(object, Z = Z, num.trees.per.proj, num.proj)
}
\arguments{
\item{object}{a crf object.}

\item{Z}{a matrix of candidate points.}

\item{num.trees.per.proj}{an integer, the number of trees per projection.}

\item{num.proj}{an integer specifying the number of projections.}
}
\value{
a numeric value, the DR I-Score.
}
\description{
compute the density ratio score
}
