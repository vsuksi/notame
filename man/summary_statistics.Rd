% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.R
\name{summary_statistics}
\alias{summary_statistics}
\title{Summary statistics}
\usage{
summary_statistics(object, grouping_cols = NA)
}
\arguments{
\item{object}{a MetaboSet object}

\item{grouping_cols}{character vector, the columns by which grouping should 
be done. Use \code{NA}
to compute statistics without grouping.}
}
\value{
A data frame with the summary statistics.
}
\description{
Computes summary statistics for each feature, possibly grouped by a factor.
The statistics include mean, standard deviation (sd), median,
median absolute deviation (mad), minimum (min), maximum (max)
as well as 25% and 75% quantiles (Q25 & Q75).
}
\examples{
# Group by "Group"
sum_stats <- summary_statistics(example_set)
# Group by Group and Time
sum_stats <- summary_statistics(example_set, 
  grouping_cols = c("Group", "Time"))
# No Grouping
sum_stats <- summary_statistics(example_set, grouping_cols = NA)

}
