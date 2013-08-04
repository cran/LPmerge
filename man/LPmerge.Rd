\name{LPmerge}
\alias{LPmerge}
\title{
Merging linkage maps by linear programming.
}
\description{
LPmerge creates a composite (or consensus) genetic map by merging
linkage maps from different populations.  The software uses linear programming (LP) to efficiently minimize the mean absolute error between the composite map and the linkage maps.  This minimization is performed subject to linear inequality constraints that ensure the ordering of the markers in the linkage maps is preserved.  When marker order is inconsistent between linkage maps, a minimum set of ordinal constraints is deleted to resolve the conflicts.
}
\usage{
LPmerge(Maps, max.interval = 1:3, weights = NULL)
}
\arguments{
  \item{Maps}{
List of \eqn{T} maps for a single linkage group.  Each genetic map is a data frame with two columns: the first column contains the marker names, the second column contains the map position.
}
  \item{max.interval}{
A whole number (\eqn{K} in formula below) specifying the maximum interval size between bins to include in the objective function.  An array of numbers can be passed to test different values (one composite map is produced for each value in the array).
}
  \item{weights}{
Numeric array of length \eqn{T}, used to weight the maps in the
objective function (see details).
}
}
\details{
Map positions in the ith linkage map are denoted by \eqn{y_i}, and composite map positions are denoted by \eqn{x}.  Within linkage map i, the markers are ordered from j = 1 to \eqn{M_i}, and the map distance between the jth and (j+q)th markers is \eqn{y_i(j+q) - y_i(j)}.  Letting \eqn{u(j;i)} denote the composite map bin containing marker j from map i, the corresponding distance in the composite map is \eqn{x(u(j+q;i))-x(u(j;i))}.  The total error across \eqn{T} maps with maximum interval size \eqn{K} is
\deqn{\sum_{i=1}^{T} W_i N^{-1}_i \sum_{q=1}^{K} \sum_{j=1}^{M_i} |x(u(j+q;i))-x(u(j;i))-[y_i(j+q)-y_i(j)]|}
where \eqn{N_i = \sum_{q=1}^{K} \sum_{j=1}^{M_i} 1} is the number of error terms for map i, and \eqn{W_i} are the weights (equal to 1 by default).  At the end of the linkage map, where the sum j+q exceeds \eqn{M_i}, the expression is evaluated as if the map were circular rather than linear.  These "wrap-around" error terms keep the total composite map length commensurate with the average length of the linkage maps.

Linear inequality constraints are used to ensure the marker order in the composite map is consistent with the order in the linkage maps.  When the linkage maps are not consistent, a minimum set of constraints is deleted to resolve the conflicts. See the reference for more details.  The deleted constraints are printed to the standard output.

One way to select the maximum interval size \eqn{K} is based on the principle of minimizing the root mean-squared error (RMSE) between the composite map and the linkage maps.  The RMSE for each linkage map, and the overall mean, is displayed for convenience.  Since the composite map length can vary with \eqn{K}, this is another factor to consider when selecting this parameter.  See the vignette for a detailed example.
}
\value{
A list with length equal to the length of the max.interval parameter.  Each entry in the list is a data frame containing the composite map and the component linkage maps.
}
\references{
Endelman,J.B. and C. Plomion (submitted)
}
\examples{
#50 random markers (out of 200 possible)
m1 <- c(0,sort(sort(runif(200),index.return=TRUE)$ix[1:50]))
m2 <- c(0,sort(sort(runif(200),index.return=TRUE)$ix[1:50]))
Maps <- list()
Maps[[1]] <- data.frame(marker=m1,map=m1)
Maps[[2]] <- data.frame(marker=m2,map=m2)
ans <- LPmerge(Maps,max.interval=1)
}
