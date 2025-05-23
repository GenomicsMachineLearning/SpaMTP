% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialMultiomicAlignment.R
\name{rigid.transf}
\alias{rigid.transf}
\title{Creates a transformation matrix for rotation and translation}
\usage{
rigid.transf(h = 0, k = 0, alpha = 0)
}
\arguments{
\item{h}{Numeric: offset along x axis}

\item{k}{Numeric: offset along y axis}

\item{alpha}{rotation angle}
}
\description{
Creates a transformation matrix for clockwise rotation by 'alpha' degrees
followed by a translation with an offset of (h, k). Points are assumed to be
centered at (0, 0).
}
