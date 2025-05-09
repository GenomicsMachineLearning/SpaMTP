% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SpatialMultiomicAlignment.R
\name{get_square_coordinates}
\alias{get_square_coordinates}
\title{This function computes the coordinates of a square's four corners based on a given center point and width.}
\usage{
get_square_coordinates(center_x, center_y, width, name)
}
\arguments{
\item{center_x}{Numeric value defining the x-coordinate of the square's center.}

\item{center_y}{Numeric value defining the y-coordinate of the square's center.}

\item{width}{Numeric value indicating the width of the square.}

\item{name}{Character string specifying the label assigned to the square for identification.}
}
\value{
A data frame with columns \code{Selection}, \code{X}, and \code{Y} representing the square's name and corner coordinates.
}
\description{
This function computes the coordinates of a square's four corners based on a given center point and width.
}
\examples{
get_square_coordinates(center_x = 5, center_y = 5, width = 4, name = "MySquare")
}
