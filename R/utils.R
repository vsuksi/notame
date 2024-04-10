#' @importFrom grDevices dev.off
#' @importFrom stats dist hclust lm na.omit runif
#' @importFrom utils capture.output citation type.convert
#' @importFrom Biobase exprs exprs<- phenoData pData pData<- featureData fData fData<- sampleNames sampleNames<- featureNames featureNames<- assayData protocolData
NULL

utils::globalVariables(c('i', '.'))

#' Set default color scales on load
#'
#' @param libname,pkgname default parameters
#' @rawNamespace import(ggplot2, except = Position)
#' @importFrom grDevices rgb
#' @noRd
.onLoad <- function(libname, pkgname) {
  op <- options()
  op_notame <- list(
    notame.citations = list(
      "Preprocessing and analyses were performed using notame package:" = utils::citation("notame"),
      "notame is built on a class from Biobase package:" = utils::citation("Biobase"),
      "visualizations in notame are built with ggplot2:" = utils::citation("ggplot2")
    ),
    notame.color_scale_con = ggplot2::scale_color_viridis_c(),
    notame.color_scale_dis = ggplot2::scale_color_brewer(palette = "Set1"),
    notame.fill_scale_con = ggplot2::scale_fill_viridis_c(),
    notame.fill_scale_dis = ggplot2::scale_fill_brewer(palette = "Set1"),
    notame.fill_scale_div_con = ggplot2::scale_fill_distiller(palette = "RdBu"),
    notame.fill_scale_div_dis = ggplot2::scale_fill_brewer(palette = "RdBu"),
    notame.shape_scale = ggplot2::scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 13))
  )
  toset <- !(names(op_notame) %in% names(op))
  if (any(toset)) options(op_notame[toset])

  invisible()
}


add_citation <- function(name, ref) {
  cites <- getOption("notame.citations")
  if (!name %in% names(cites)) {
    cites[[name]] <- ref
    options(notame.citations = cites)
  }
}

#' Show citations
#'
#' This function lists citations for all the major packages used by the notame functions that
#' have been called during the session. All notame functions update the list automatically.
#' The citations are taken from the call to \code{citation("package")}, and complemented with
#' a brief description of what the package was used for.
#' NOTE: the citations might not point to the correct paper if the package authors have not
#' supplied correct citation information for their package.
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
    log_text(capture.output(show(cites[[i]])))
  }
}


#' Summary statistics of finite elements
#'
#' These functions first remove non-finite and missing values, then compute the summary statistic in question.
#' They are helper functions used for computing quality measurements.
#' 
#' @param x a numeric vector.
#' @param ... other parameters passed to underlying function
#' @return a named, numeric vector with the summary statistic in question.
#'
#' @name finite_helpers
NULL

#' @export
#' @rdname finite_helpers
finite_sd <- function(x) {
  sd(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x[is.finite(x)], na.rm = TRUE)
}

#' @examples
#' data <- combined_data(example_set)
#' features <- fData(example_set)
#' features$MPA <- sapply(data[, features[, "Feature_ID"]], finite_median)
#'
#' @importFrom stats median
#' @export
#' @rdname finite_helpers
finite_median <- function(x) {
  median(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_min <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  min(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_max <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  max(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_mad <- function(x) {
  mad(x[is.finite(x)], center = median(x[is.finite(x)], na.rm = TRUE), na.rm = TRUE)
}

#' @importFrom stats quantile
#' @export
#' @rdname finite_helpers
finite_quantile <- function(x, ...) {
  unname(quantile(x[is.finite(x)], na.rm = TRUE, ...))
}



# Defaults for NULL values
`%||%` <- function(a, b) {
  suppressWarnings(if (is.null(a)) {
    b
  } else if (is.na(a)) {
    b
  } else {
    a
  })
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
#' @export
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
#' @export
prop_found <- function(x) {
  sum(!is.na(x)) / length(x)
}

best_class <- function(x) {
  x <- type.convert(as.character(x), as.is = TRUE)
  if (is(x, "numeric")) {
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


best_classes <- function(x) {
  as.data.frame(lapply(x, best_class), stringsAsFactors = FALSE)
}

all_unique <- function(x) {
  !any(duplicated(x))
}
