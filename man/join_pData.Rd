% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/class_constructor.R
\name{join_pData}
\alias{join_pData}
\title{Join new columns to pheno data}
\usage{
join_pData(object, dframe)
}
\arguments{
\item{object}{a MetaboSet object}

\item{dframe}{a data frame with the new information}
}
\value{
A MetaboSet object with the new information added to pData(object).
}
\description{
Join a new data frame of information to pheno data of a MetaboSet object.
}
\examples{
new_info <- data.frame(
  Sample_ID = sampleNames(example_set),
  BMI = stats::runif(ncol(example_set), 22, 26)
)
with_new_info <- join_pData(example_set, new_info)
colnames(pData(with_new_info))

}
