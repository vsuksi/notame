
#' \code{notame} package.
#'
#' Provides functionality for untargeted LC-MS metabolomics research as 
#' specified in the associated publication in the 'Metabolomics Data
#' Processing and Data Analysis—Current Best Practices' special issue of the 
#' Metabolites journal (2020). This includes data pretreatment and quality 
#' control, uni- and multivariate analysis as well as quality control 
#' visualizations, feature-wise visualizations and results visualizations. Raw 
#' data preprocessing and functionality related to biological context, such as 
#' pathway analysis, is not included.
#'
#' @details
#' In roughly chronological order, the functionality of notame is as follows. 
#' Please see the vignettes and paper for more information.
#' 
#' Data pretreatment (for reducing unwanted variation and completing the 
#' dataset):
#' \itemize{
#' \item \code{\link{mark_nas}} mark specified values as missing
#' \item \code{\link{flag_detection}} flag features with low detection rate
#' \item \code{\link{flag_contaminants}} flag contaminants based on blanks
#' \item \code{\link{correct_drift}} QC-based cubic spline drift correction
#' \item \code{\link{assess_quality}} assess quality information of features
#' \item \code{\link{quality}} extract quality information of features
#' \item \code{\link{flag_quality}} flag low-quality features
#' \item \code{\link{flag_report}} report of flagged features
#' \item \code{\link{impute_rf}} impute missing values using random forest
#' \item \code{\link{impute_simple}} impute missing values using simple 
#' imputation
#' \item \code{\link{cluster_features}} cluster correlated features originating 
#' from the same metabolite
#' \item \code{\link{align_batches}} align features between batches
#' \item \code{\link{normalize_batches}}, \code{\link{ruvs_qc}} normalize 
#' between batches
#' \item \code{\link{pca_bhattacharyya_dist}}, 
#' \code{\link{perform_repeatability}} evaluate between-batch normalization
#' }
#'
#' Data pretreatment visualizations (for visualization of data pretreatment and 
#' exploring the data at large):
#' \itemize{
#' \item \code{\link{visualizations}} write all relevant data pretreatment 
#' visualizations to pdf
#' \item \code{\link{plot_injection_lm}} p-value histograms relating 
#' injection order to each feature
#' \item \code{\link{plot_sample_boxplots}} sample boxplots relating 
#' abundance to, for example, study group or injection order
#' \item \code{\link{plot_dist_density}} density plot of euclidean 
#' distances between samples
#' \item \code{\link{plot_tsne}}, \code{\link{plot_pca}} t-SNE and PCA 
#' plots colored by, for example, study group or injection order
#' \item \code{\link{plot_tsne_arrows}}, \code{\link{plot_pca_arrows}} t-
#' SNE and PCA plots with arrows representing the change in samples across time 
#' points 
#' \item \code{\link{plot_tsne_hexbin}}, \code{\link{plot_pca_hexbin}}
#' hexbin t-SNE and PCA plots colored by injection order
#' \item \code{\link{plot_dendrogram}} sample dendrogram, hierarchically 
#' clustered using Ward's criterion on Euclidean distances between samples
#' \item \code{\link{plot_sample_heatmap}} heatmap of Euclidean distances 
#' between samples, ordered by hierarchical clustering using Ward's criterion
#' \item \code{\link{plot_pca_loadings}} PCA loadings plot
#' \item \code{\link{plot_sample_suprahex}} suprahex plots of samples 
#' \item \code{\link{plot_quality}} distributions of quality metrics plot 
#' }
#'
#' Transformations:
#' \itemize{
#' \item \code{\link{log}} logarithms
#' \item \code{\link{exponential}} exponential
#' \item \code{\link{pqn_normalization}} probabilistic quotient normalization
#' \item \code{\link{inverse_normalize}} inverse-rank normalization
#' \item \code{\link{scale}} scale
#' }
#' 
#' Univariate analysis (return data.frames):
#' \itemize{
#' \item \code{\link{perform_lm}} linear models 
#' \item \code{\link{perform_logistic}} logistic regression
#' \item \code{\link{perform_lmer}} linear mixed models 
#' \item \code{\link{perform_oneway_anova}} Welch’s ANOVA and classic ANOVA
#' \item \code{\link{perform_lm_anova}} linear models ANOVA table
#' \item \code{\link{perform_t_test}} Welch's and Student's t-test 
#' \item \code{\link{perform_pairwise_t_test}} pairwise and paired t-test
#' \item \code{\link{perform_kruskal_wallis}} Kruskal-Wallis rank-sum test
#' \item \code{\link{perform_mann_whitney}} Mann-Whitney u test
#' \item \code{\link{perform_wilcoxon_signed_rank}} Wilcoxon signed rank test
#' \item \code{\link{perform_pairwise_non_parametric}} pairwise non-parametric 
#' tests
#' \item \code{\link{perform_correlation_tests}} correlations tests
#' \item \code{\link{perform_auc}} area under curve for each subject and feature
#' \item \code{\link{perform_homoscedasticity_tests}} test homoscedasticity
#' \item \code{\link{cohens_d}} Cohen's D
#' \item \code{\link{fold_change}} fold change
#' \item \code{\link{summary_statistics}} summary statistics
#' \item \code{\link{clean_stats_results}} statistics cleaning
#' }
#'
#' Multivariate analysis (return various objects):
#' \itemize{
#' \item \code{\link{muvr_analysis}} multivariate modelling with minimally 
#' biased variable selection (MUVR2) with random forest, PLS(-DA)
#' and elastic net regression
#' \item \code{\link{mixomics_pls}}, \code{\link{mixomics_plsda}} PLS(-DA)
#' \item \code{\link{mixomics_pls_optimize}},
#' \code{\link{mixomics_plsda_optimize}} optimal number of components for
#' PLS(-DA)
#' \item \code{\link{mixomics_spls_optimize}},
#' \code{\link{mixomics_spls_optimize}} optimal number of components 
#' and features in each component for PLS(-DA)
#' \item \code{\link{fit_rf}}, \code{\link{importance_rf}} random forest
#' and feature importance
#'  \item \code{\link{perform_permanova}} PERMANOVA
#' }
#'
#' Feature-wise visualizations (these are often drawn for a subset of 
#' interesting features after analysis):
#' \itemize{
#' \item \code{\link{save_beeswarm_plots}} save beeswarm plots of each feature 
#' \item \code{\link{save_group_boxplots}} save box plots of each feature
#' \item \code{\link{save_scatter_plots}} save scatter plots of each feature 
#' against a set variable
#' \item \code{\link{save_subject_line_plots}} save line plots with mean
#' \item \code{\link{save_group_lineplots}} save line plots with error bars by 
#' group
#' \item \code{\link{save_dc_plots}} save drift correction plots
#' \item \code{\link{save_batch_plots}} save batch correction plots
#' }
#'
#' Results visualizations (save them using \code{\link{save_plot}}):
#' \itemize{
#' \item \code{\link{plot_p_histogram}} histogram of p-values
#' \item \code{\link{volcano_plot}} volcano plot
#' \item \code{\link{manhattan_plot}} manhattan plot
#' \item \code{\link{mz_rt_plot}} m/z vs retention time plot (aka cloud plot)
#' \item \code{\link{plot_effect_heatmap}} heatmap of effects between 
#' variables, such as correlations
#' }
#' 
#' MetaboSet utilities:
#' \itemize{
#' \item \code{\link{read_from_excel}} read formatted Excel files
#' \item \code{\link{construct_metabosets}} construct MetaboSet objects
#' \item \code{\link{write_to_excel}} write results to Excel file
#' \item \code{\link{group_col}} get and set name of the special column for 
#' group labels
#' \item \code{\link{time_col}} get and set the name of the special column for 
#' time points
#' \item \code{\link{subject_col}} get and set the name of the special column 
#' for subject identifiers
#' \item \code{\link{flag}} get and set values in the flag column
#' \item \code{\link{drop_flagged}} drop flagged features
#' \item \code{\link{drop_qcs}} drop QC samples
#' \item \code{\link{join_fData}} join new columns to feature data
#' \item \code{\link{join_pData}} join new columns to pheno data
#' \item \code{\link{combined_data}} retrieve both sample information and 
#' features
#' }
#' 
#' Other utilities:
#' \itemize{
#' \item \code{\link{citations}} show citations
#' \item \code{\link{init_log}} initialize log to a file
#' \item \code{\link{log_text}} log text to the current log file
#' \item \code{\link{finish_log}} finish a log
#' \item \code{\link{save_plot}} save plot to chosen format
#' \item \code{\link{fix_MSMS}} transform the MS/MS output to publication ready
#' \item \code{\link{merge_metabosets}} merge MetaboSet objects together
#' }
#'
#' @name notame-package
#' @references
#' Klåvus et al. (2020). "notame": Workflow for Non-Targeted LC-MS Metabolic 
#' Profiling. Metabolites, 10: 135.
"_PACKAGE"
NULL


#' @rawNamespace import(ggplot2, except = Position)
#' @importFrom utils citation
#' @importFrom Biobase exprs exprs<- phenoData pData pData<- featureData fData 
#' fData<- sampleNames sampleNames<- featureNames featureNames<- assayData 
#' protocolData
#' @importFrom magrittr "%>%" "%<>%"
#' @importClassesFrom Biobase ExpressionSet
#' @import BiocGenerics
#' @import methods
NULL

utils::globalVariables(c('i', '.'))

#' Set default color scales on load
#'
#' @param libname,pkgname default parameters
#' @noRd
.onLoad <- function(libname, pkgname) {
  op <- options()
  op_notame <- list(
    notame.citations = list(
      "Preprocessing and analyses were performed using notame package:" =
      utils::citation("notame"),
      "notame is built on a class from Biobase package:" =
      utils::citation("Biobase"),
      "visualizations in notame are built with ggplot2:" =
      utils::citation("ggplot2")),
    notame.color_scale_con = scale_color_viridis_c(),
    notame.color_scale_dis = scale_color_brewer(palette = "Set1"),
    notame.fill_scale_con = scale_fill_viridis_c(),
    notame.fill_scale_dis = scale_fill_brewer(palette = "Set1"),
    notame.fill_scale_div_con = scale_fill_distiller(palette = "RdBu"),
    notame.fill_scale_div_dis = scale_fill_brewer(palette = "RdBu"),
    notame.shape_scale = scale_shape_manual(values = c(16, 17, 15, 3,
                                                       7, 8, 11, 13)))
  toset <- !(names(op_notame) %in% names(op))
  if (any(toset)) {
    options(op_notame[toset])
  }
  invisible()
}

.add_citation <- function(name, ref) {
  cites <- getOption("notame.citations")
  if (!name %in% names(cites)) {
    cites[[name]] <- ref
    options(notame.citations = cites)
  }
}

#' Show citations
#'
#' This function lists citations for all the major packages used by the notame 
#' functions that have been called during the session. All notame functions 
#' update the list automatically. The citations are taken from the call to 
#' '\code{citation("package")}, and complemented with a brief description of 
#' what the package was used for.
#' NOTE: the citations might not point to the correct paper if the package 
#' authors have not supplied correct citation information for their package.
#' The output is written to the current log file, if specified.
#'
#' @return None, the function is invoked for its side effect.
#'
#' @examples
#'
#' citations()
#' plot_tsne(merged_sample)
#' # Rtsne added to citations
#' citations()
#'
#' @export
citations <- function() {
  cites <- getOption("notame.citations")
  for (i in seq_along(cites)) {
    log_text(names(cites)[i])
    log_text(utils::capture.output(show(cites[[i]])))
  }
}


#' Summary statistics of finite elements
#'
#' These functions first remove non-finite and missing values, then 
#' compute the summary statistic in question. They are helper 
#' functions used for computing quality measurements.
#' 
#' @param x a numeric vector.
#' @param ... other parameters passed to underlying function
#' @return A named, numeric vector with the summary statistic in question.
#'
#' @name finite_helpers
#' @noRd
NULL

finite_sd <- function(x) {
  sd(x[is.finite(x)], na.rm = TRUE)
}

finite_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x[is.finite(x)], na.rm = TRUE)
}

finite_median <- function(x) {
  stats::median(x[is.finite(x)], na.rm = TRUE)
}

finite_min <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  min(x[is.finite(x)], na.rm = TRUE)
}

finite_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x[is.finite(x)], na.rm = TRUE)
}

finite_mad <- function(x) {
  mad(x[is.finite(x)], 
      center = stats::median(x[is.finite(x)], na.rm = TRUE), 
      na.rm = TRUE)
}

finite_quantile <- function(x, ...) {
  unname(stats::quantile(x[is.finite(x)], na.rm = TRUE, ...))
}

# Defaults for NULL values
`%||%` <- function(a, b) {
  if (is.null(a)) {
    b
  } else if (is.na(a)) {
    b
  } else {
    a
  }
}

#' Proportion of NA values in a vector
#'
#' @param x a numeric vector
#'
#' @return A numeric, the proportion of non-missing values in a vector.
#'
#' @examples
#' example_set <- mark_nas(example_set)
#' prop_na(exprs(example_set))
#' 
#' @noRd
prop_na <- function(x) {
  sum(is.na(x)) / length(x)
}

#' Proportion of non-missing values in a vector
#'
#' @param x a numeric vector
#' 
#' @return A numeric, the proportion of non-missing values in vector.
#'
#' @examples
#' example_set <- mark_nas(example_set, value = 0)
#' prop_found(exprs(example_set))
#'
#' @noRd
prop_found <- function(x) {
  sum(!is.na(x)) / length(x)
}

.best_class <- function(x) {
  x <- utils::type.convert(as.character(x), as.is = TRUE)
  if (inherits(x, "numeric")) {
    x <- x
  } else if (length(unique(x)) < length(x) / 4) {
    x <- as.factor(x)
  } else if (is.integer(x)) {
    x <- as.numeric(x)
  } else {
    x <- as.character(x)
  }
  x
}

.best_classes <- function(x) {
  as.data.frame(lapply(x, .best_class), stringsAsFactors = FALSE)
}

.all_unique <- function(x) {
  !any(duplicated(x))
}
