% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multivariate.R
\name{muvr_analysis}
\alias{muvr_analysis}
\title{Multivariate modelling with minimally biased variable selection (MUVR)}
\usage{
muvr_analysis(
  object,
  y = NULL,
  id = NULL,
  multi_level = FALSE,
  multi_level_var = NULL,
  covariates = NULL,
  static_covariates = NULL,
  all_features = FALSE,
  nRep = 50,
  nOuter = 6,
  nInner = nOuter - 1,
  varRatio = 0.75,
  method = c("PLS", "RF"),
  ...
)
}
\arguments{
\item{object}{a MetaboSet object}

\item{y}{character, column name in pData of the target variable to predict}

\item{id}{character, column name in pData of the subject ID variable in case 
of repeated measurements}

\item{multi_level}{logical, whether multi-level modeling should be applied, 
see Details}

\item{multi_level_var}{character, column name in pData of the variable for 
splitting the data in multi-level modeling}

\item{covariates, static_covariates}{character, column names of pData to use 
as covariates in the model, in addition to molecular features. 
\code{static_covariates} are ignored for non-multi-level models.
For multi-level models, the change in \code{covariates} is computed, while
\code{static_covariates} are taken from the first time point.}

\item{all_features}{logical, should all features be included in the model? 
if FALSE, flagged features are left out}

\item{nRep}{Number of repetitions of double CV, parameter of MUVR}

\item{nOuter}{Number of outer CV loop segments, parameter of MUVR}

\item{nInner}{Number of inner CV loop segments, parameter of MUVR}

\item{varRatio}{Ratio of variables to include in subsequent inner loop 
iteration, parameter of MUVR}

\item{method}{Multivariate method. Supports 'PLS', 'RF' and 'EN'}

\item{...}{other parameters to MUVR2::MUVR2 or MUVR2::MUVR2_EN 
and MUVR2::getVar (when method == "EN")}
}
\value{
A MUVR object. (make this more descriptive)
}
\description{
A wrapper around the MUVR2 (random forest, PLS(-DA)) and MUVR2_EN (elastic 
net) functions from the MUVR2 package.
}
\details{
This function is now using the MUVR2 package, characterized as an 
upgrade extending the original MUVR package by the inclusion of elastic net 
regression (EN) and other functionality not covered by this wrapper. Elastic 
net regression supports covariate adjustment by suppressing regularization 
of specified features from the regularization procedure. Note that this is 
different from simply including covariates such as sex. EN also differs from 
PLS and RF in that no recursive variable elimination is performed, so an 
additional scheme is used to obtain the 'min', 'mid' and 'max' models using
MUVR2::getVar in this wrapper. 

Sex would be entered as a static covariate, since the change in sex is zero 
for all individuals, so computing the change and using that as a covariate 
does not make sense.

Note that there are several more plots available in MUVR2 for inspecting the 
results, notably MUVR2::plotMV, MUVR2::plotStability and MUVR2::plotVIRank. 
Many of these return different plots depending on the model specification.
}
\examples{
# PLS simple
rf_model <- muvr_analysis(drop_qcs(merged_sample), 
  y = "Group", nRep = 2, method = "PLS")

# RF with covariate and repeated measures (not longitudinal)
ex_set <- drop_qcs(example_set)
ex_set$Injection_order \%<>\% as.numeric()
en_model <- muvr_analysis(drop_qcs(ex_set), y = "Group", id = "Subject_ID", 
  nRep = 2, method = "RF", covariates = "Injection_order")

# RF on multilevel variable comparing levels of y
rf_model_ <- muvr_analysis(drop_qcs(example_set), 
  y = "Group", multi_level = TRUE, id = "Subject_ID", 
  multi_level_var = "Time", method = "RF", nRep = 2)

# EN on multilevel variable with covariate and static covariate
ex_set <- drop_qcs(example_set)
example_set$Injection_order \%<>\% as.numeric()
example_set$Group \%<>\% as.numeric()
en_model <- muvr_analysis(drop_qcs(example_set), id = "Subject_ID", 
 multi_level = TRUE, multi_level_var = "Time", 
 covariates = "Injection_order", static_covariates = "Group", 
 method = "EN", nRep = 2)

}
\seealso{
\code{\link[MUVR2]{MUVR2}} \code{\link[MUVR2]{MUVR2_EN}} 
\code{\link[MUVR2]{getVar}} \code{\link[MUVR2]{plotMV}} 
\code{\link[MUVR2]{plotStability}} \code{\link[MUVR2]{plotVIRank}}
}
