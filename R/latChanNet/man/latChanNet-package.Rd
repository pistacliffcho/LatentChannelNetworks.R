\name{latChanNet-package}
\alias{latChanNet-package}
\alias{latChanNet}
\docType{package}
\title{
  Latent Network Models for edge and metadata prediction.
}
\description{
  Fits Latent Channel Networks and the Poisson model of Ball, Karrer and Newman (2011). Allows for unknown edges statuses. Augments network with metadata to allow for metadata predictions. 
}
\details{
  Models are built (but not fitted) with \code{makeLatentModel}. Models are fit with the \code{$fit()} method. Fitted parameters can be extracting via the \code{$get_pars()} method. Predictions of both edges and metadata can be done with \code{predict}. Heatmaps of parameters can be plotted with \code{plot_net}. 
}
\author{
Clifford Anderson-Bergman.

Maintainer: Clifford Anderson-Bergman <pistacliffcho@gmail.com>
}
\references{
  Clifford Anderson-Bergman, Phan Nguyen, and Jose Cadena Pico. "Latent Channel Networks", submitted 2019

BKN model:

Brian Ball, Brian Karrer, and Mark EJ Newman. "Efficient and principled method for detecting communities in networks." Physical Review E 84.3 (2011): 036103.
}
\keyword{ network }
