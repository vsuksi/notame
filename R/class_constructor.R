.log_text_if <- function(text, logif) {
  if (logif) {
    log_text(text)
  }
}

.create_sample_col <- function(x, id_prefix, id_column, log_messages) {
  # If id_column is provided, try to change name of the column to "Sample_ID"
  if (!is.null(id_column)) {
    .log_text_if("Checking provided sample ID column", log_messages)
    if (!id_column %in% colnames(x)) {
      .log_text_if(paste0("ID column '", id_column, "' not found"),
                   log_messages)
    } else if (!any(duplicated(x[, id_column])) && !any(is.na(x[, id_column]))){
      x$Sample_ID <- x[, id_column]
      .log_text_if(paste0("Column 'Sample_ID' created from ", id_column),
                   log_messages)
    } else {
      .log_text_if("Provided sample ID column is not valid", log_messages)
    }
  }
  # If Sample_ID is not provided explicitly, it will be created
  if (!"Sample_ID" %in% colnames(x)) {
    x$Sample_ID <- paste0(id_prefix, x$Injection_order)
    .log_text_if(paste("Sample ID generated from injection orders and prefix",
                       id_prefix), 
                 log_messages)
  } else {
    # Fill sample IDs of "QC" identifiers
    x$Sample_ID <- as.character(x$Sample_ID)
    if (any(is.na(x$Sample_ID))) {
      .log_text_if(paste("Missing values found in Sample_ID,",
                         "filling missing sample IDs of QC samples."), 
                   log_messages)
      x$Sample_ID[is.na(x$Sample_ID) & x$QC == "QC"] <- "QC"
    }
    # Fill sample IDs of "Blank" identifiers
    if (any(is.na(x$Sample_ID))) {
      .log_text_if(paste("Missing values found in Sample_ID after filling IDs",
                         "of QCs, trying to detect 'Blank' samples"), 
                  log_messages)
      blank_found <- which(vapply(x, function(y) {
        any(y == "Blank")
      }, logical(1)))
      if (length(blank_found)) {
        x$Sample_ID[is.na(x$Sample_ID) &
        x[[blank_found[1]]] == "Blank"] <- "Blank"
        x$Sample_ID[x$Sample_ID == "Blank"] <- paste0("Blank_",
          seq_len(sum(x$Sample_ID == "Blank")))
      }
      if (any(is.na(x$Sample_ID))) {
        stop("Missing values in Sample_ID after attempting",
             " filling IDs of QCs and Blanks")
      }
    }
    # Add running index to QC identifiers
    if (sum(x$Sample_ID == "QC")) {
      .log_text_if("Adding running index to 'QC' sample IDs", log_messages)
      x$Sample_ID[x$Sample_ID == "QC"] <- paste0("QC_", 
        seq_len(sum(x$Sample_ID == "QC")))
    }
    # Add running index to Blank identifiers
    if (sum(x$Sample_ID == "Blank")) {
      .log_text_if("Adding running index to 'Blank' sample IDs", log_messages)
      x$Sample_ID[x$Sample_ID == "Blank"] <- 
        paste0("Blank_", seq_len(sum(x$Sample_ID == "Blank")))
    }
    # After this, the Sample IDs should be unique
    if (any(duplicated(x$Sample_ID))) {
      stop("Sample_ID is not unique")
    }
  }
  x
}

## Helper function for checking integrity of pheno data
.check_pheno_data <- function(x, id_prefix, 
                             id_column = NULL, log_messages = FALSE) {
  .log_text_if("\nChecking sample information", log_messages)
  # Check that Injection order is included
  if (!"Injection_order" %in% colnames(x)) {
    stop("'Injection_order' not found for the samples.")
  }
  # Check that injection order can be numeric
  if (!is.numeric(x$Injection_order) && !.looks_numeric(x$Injection_order)) {
    stop("'Injection_order' is not numeric and cannot be converted to numeric")
  }
  # No NAs allowed in Injection order
  if (any(is.na(x$Injection_order))) {
    stop("Missing values in Injection_order")
  }
  # Injection order should be unique
  if (length(unique(x$Injection_order)) != nrow(x)) {
    stop("Injection_order is not unique")
  }
  # If QC column is not provided explicitly, attempt to create it
  if (!"QC" %in% colnames(x)) {
    qc_found <- apply(x, 1, function(y) {
      any(grepl("QC", y))
    })
    if (any(qc_found)) {
      x$QC <- ifelse(qc_found, "QC", "Sample")
      .log_text_if(paste("QC column generated from rows containing 'QC'"),
                  log_messages)
    } else {
      warning("QC not found and column can not be generated.", 
              " Please create one before constructing a MetaboSet object.")
    }
  }
  # Create and populate 'Sample_ID' column
  x <- .create_sample_col(x, id_prefix, id_column, log_messages)
  # Select best classes for columns and prepare data.frame
  x <- .best_classes(x)
  rownames(x) <- x$Sample_ID
  x <- as.data.frame(dplyr::select(x, "Sample_ID", dplyr::everything()))
  
  x
}

## Check that the position of the corner row and column is OK
.check_position <- function(x, cc, cr) {
  condition <- is.na(x[cr - 1, cc - 1]) &
    (is.numeric(utils::type.convert(x[cr + 1, cc + 1], as.is = TRUE)) |
      is.na(x[cr + 1, cc + 1])) &
    !is.na(x[cr, cc]) &
    !is.na(x[cr - 1, cc]) &
    !is.na(x[cr, cc - 1])
  if (!condition) {
    stop("The corner row and column coordinates seem to be incorrect!")
  }
}

# Check if a vector can be safely converted to numeric
.looks_numeric <- function(x) {
  stopifnot(is.atomic(x) || is.list(x))
  tryCatch(
    {
      as.numeric(x)
      TRUE
    },
    warning = function(w) FALSE
  )
}

## Check that all abundances look OK
.check_exprs <- function(exprs_, log_messages = FALSE) {
  .log_text_if("Checking that feature abundances only contain numeric values",
               log_messages)
  # Check that all rows are full of numbers
  non_numerics <- exprs_ %>% apply(1, function(x) !.looks_numeric(x))
  if (sum(non_numerics)) {
    stop(paste("Non-numeric values found in the abundances on rows",
               paste(which(non_numerics), collapse = ", ")))
  }
  # Convert to numeric
  exprs_ <- exprs_ %>% apply(2, as.numeric)
    
  exprs_
}

.check_feature_data <- function(feature_data, check_limits = TRUE, 
                                mz_limits = c(10, 2000), rt_limits = c(0, 20),
                                log_messages = FALSE) {
  .log_text_if("\nChecking feature information", log_messages)
  .log_text_if("Checking that feature IDs are unique and not stored as numbers",
               log_messages)
  fid <- feature_data$Feature_ID
  if (any(duplicated(fid))) {
    stop("Feature_ID values are not unique")
  }
  if (any(is.na(fid))) {
    stop("Missing values in Feature IDs")
  }
  fid_num <- withCallingHandlers(
    expr = as.numeric(fid), 
    warning = function(w) tryInvokeRestart("muffleWarning"))
  if (any(!is.na(fid_num))) {
    stop("Numbers are not allowed as feature IDs.", call. = FALSE)
  }
  fid_chr <- withCallingHandlers(
    expr = as.character(fid),
    warning = function(w) tryInvokeRestart("muffleWarning"))
  if (any(grepl("^[[:digit:]]", fid_chr))) {
    stop("Feature IDs can not start with numbers.", call. = FALSE)
  }
  if (check_limits) {
    .log_text_if("Checking that m/z and retention time values are reasonable.",
                 log_messages)
    mz <- feature_data[, .find_mz_rt_cols(feature_data)$mz_col]
    rt <- feature_data[, .find_mz_rt_cols(feature_data)$rt_col]
    if (!(all(mz > mz_limits[1]) && all(mz < mz_limits[2])) ||
        !(all(rt > rt_limits[1]) && all(rt < rt_limits[2]))) {
      stop("Values in m/z or retention time columns are outside limits.")
    }
  }

  feature_data
}

.extract_information <- function(file, sheet, corner_row, corner_column, name) {
  dada <- openxlsx::read.xlsx(file, sheet, colNames = FALSE)
  
  # Define excel column order A-Z, AA - ZZ, AAA - ZZZ
  grid2 <- expand.grid(LETTERS, LETTERS)
  combinations2 <- paste0(grid2$Var2, grid2$Var1)
  grid3 <- expand.grid(LETTERS, combinations2)
  combinations3 <- paste0(grid3$Var2, grid3$Var1)
  excel_columns <- c(LETTERS, combinations2, combinations3)
  
  # If corner coordinates are omitted, try to find them automatically
  if (is.null(corner_row) || is.null(corner_column)) {
    log_text("Detecting corner position")
  }
  corner_row <- corner_row %||% which(!is.na(dada[, 1]))[1]
  corner_column <- corner_column %||% which(!is.na(dada[1, ]))[1]
  # Column can be given as a character
  cc <- ifelse(is.character(corner_column),
    which(excel_columns == corner_column), corner_column)
  cr <- corner_row
  # Check that corner is in the right spot
  .check_position(dada, cc, cr)
  log_text(paste0("Corner detected correctly at row ", cr, ", column ",
                  excel_columns[cc]))
  # Extract sample information
  log_text(paste0("\nExtracting sample information from rows 1 to ", cr,
                  " and columns ", excel_columns[cc + 1], " to ",
                  excel_columns[ncol(dada)]))
  pheno_data <- as.data.frame(t(dada[seq_len(cr), (cc + 1):ncol(dada)]),
                              stringsAsFactors = FALSE)
  log_text(paste0("Replacing spaces in sample information column names with", 
           " underscores (_)"))
  colnames(pheno_data) <- gsub(" ", "_",
                               c(dada[seq_len(cr - 1), cc], "Datafile"))
  
  # If a single mode is given, datafile will indicate the mode
  if (!is.null(name)) {
    colnames(pheno_data)[ncol(pheno_data)] <- paste0(name, "_Datafile")
    log_text(paste0('Naming the last column of sample information "', name,
                    '_Datafile"'))
  } else {
    log_text('Naming the last column of sample information "Datafile"')
  }
  
  # Extract feature information
  log_text(paste0("\nExtracting feature information from rows ", cr + 1, 
                  " to ", nrow(dada), " and columns ", excel_columns[1], 
                  " to ", excel_columns[cc]))
  feature_data <- dada[(cr + 1):nrow(dada), seq_len(cc)]
  colnames(feature_data) <- dada[cr, seq_len(cc)]
  # Extract LC-MS measurements as matrix
  log_text(paste0("\nExtracting feature abundances from rows ", cr + 1, " to ",
                  nrow(dada), " and columns ", excel_columns[cc + 1], " to ",
                  excel_columns[ncol(dada)]))
  exprs_ <- dada[(cr + 1):nrow(dada), (cc + 1):ncol(dada)]
  list("pheno_data" = pheno_data, "feature_data" = feature_data, 
       "exprs_" = exprs_ )
}
  
.prepare_feature_data <- function(name, split_by, feature_data) {
  # If the file only contains one mode, add the mode name as Split column
  if (!is.null(name)) {
    log_text(paste0("Assigning ", name,
                    " as the value of the Split column for each feature"))
    feature_data$Split <- name
    split_by <- "Split"
  } else { # Multiple modes in the file, create Split column to separate modes
    if (!all(split_by %in% colnames(feature_data))) {
      stop(paste0("Couldn't find column(s): ",
                  paste(split_by[!(split_by %in% colnames(feature_data))],
                        collapse = ", ")))
    }
    log_text(paste0("Creating Split column from ",
                    paste0(split_by, collapse = ", ")))
    feature_data <- feature_data %>%
      tidyr::unite("Split", split_by, remove = FALSE)
  }

  # Create feature ID if necessary
  if (!"Feature_ID" %in% colnames(feature_data)) {
    log_text("Feature_ID column not found, creating feature IDs")
    feature_data <- .name_features(feature_data = feature_data)
  }
  # Reorganise columns and change classes
  feature_data <- feature_data %>%
    dplyr::select("Feature_ID", "Split", dplyr::everything()) %>%
    .best_classes() %>%
    dplyr::mutate_if(is.factor, as.character)
  rownames(feature_data) <- feature_data$Feature_ID
  log_text(paste0("Replacing dots (.) in feature information column names with",
                  " underscores (_)"))
  colnames(feature_data) <- gsub("[.]", "_", colnames(feature_data)) %>%
    # remove duplicate underscores
    gsub("_{2,}", "_", .)
    
  feature_data
}

#' Read formatted Excel files
#'
#' Reads data from an Excel file of the following format:
#' \itemize{
#'   \item Left side of the sheet contains information about the features, size 
#' features x feature info columns
#'   \item Top part contains sample information, size sample info variables x 
#' samples
#'   \item The middle contains the actual abundances, size features x samples
#' }
#' This function separates the three parts from the file, and returns them in a 
#' list.
#'
#' @param file path to the Excel file
#' @param sheet the sheet number or name
#' @param id_column character, column name for unique identification of samples
#' @param corner_row integer, the bottom row of sample information,
#' usually contains data file names and feature info column names.
#' If set to NULL, will be detected automatically.
#' @param corner_column integer or character, the corresponding column number 
#' or the column name (letter) in Excel.
#' If set to NULL, will be detected automatically.
#' @param id_prefix character, prefix for autogenerated sample IDs, see Details
#' @param split_by character vector, in the case where all the modes are in the 
#' same Excel file, the column names of feature data used to separate the modes 
#' (usually Mode and Column)
#' @param name in the case where the Excel file only contains one mode, the 
#' name of the mode, such as "Hilic_neg"
#' @param mz_limits numeric vector of two, all m/z values should be in between 
#' these
#' @param rt_limits numeric vector of two, all retention time values should be 
#' in between these
#' @param skip_checks logical: skip checking data integrity. Not recommended, 
#' but sometimes useful when you just want to read the data in as is and fix 
#' errors later. NOTE: Sample_ID and QC columns will not be constructed.
#' The data integrity checks need to be passed when contstructing MetaboSet 
#' objects.
#'
#' @inherit construct_metabosets return examples
#'
#' @return A list of three data frames:
#' \itemize{
#'   \item exprs: the actual abundances, size features x samples
#'   \item pheno_data: sample information, size sample info variables x samples
#'   \item feature_data: information about the features, size features x 
#' feature info columns
#' }
#'
#' @details Only specify one of \code{split_by} and \code{name}. The feature 
#' data returned will contain a column named "Split", which is used to separate 
#' features from different modes. Unless a column named "Feature_ID"
#' is found in the file, a feature ID will be generated based on the value of 
#' "Split", mass and retention time. The function will try to find columns for 
#' mass and retention time by looking at a few common alternatives, and throw 
#' an error if no matching column is found. Sample information needs to contain 
#' a row called "Injection_order", and the values need to be unique. In 
#' addition, a possible sample identifier row needs to be named "Sample_ID",
#' or to be specified in \code{id_column}, and the values need to be unique, 
#' with an exception of QC samples: if there are any "QC" identifiers, they 
#' will be replaced with "QC_1", "QC_2" and so on.
#' If a "Sample_ID" row is not found, it will be created using the 
#' \code{id_prefix} and injection order.
#'
#' @export
read_from_excel <- function(file, sheet = 1, id_column = NULL, 
                            corner_row = NULL, corner_column = NULL,
                            id_prefix = "ID_", split_by = NULL, name = NULL,
                            mz_limits = c(10, 2000), rt_limits = c(0, 20),
                            skip_checks = FALSE) {
  if (is.null(split_by) && is.null(name)) {
    stop("Either name or split_by needs to be defined, see documentation.")
  } else if ((!is.null(split_by)) && (!is.null(name))) {
    stop("Only define split_by OR name, see documentation.")
  }
  extracted <- .extract_information(file, sheet, corner_row, 
                                    corner_column, name)
  pheno_data <- extracted$pheno_data
  feature_data <- .prepare_feature_data(name, split_by, extracted$feature_data)
  exprs_ <- extracted$exprs_
  # Skip checks
  if (!skip_checks) {
    pheno_data <- .check_pheno_data(x = pheno_data, id_prefix = id_prefix,
                                   id_column = id_column, log_messages = TRUE)
    exprs_ <- .check_exprs(exprs_, log_messages = TRUE)
    feature_data <- .check_feature_data(feature_data, 
                                       mz_limits = mz_limits, 
                                       rt_limits = rt_limits, 
                                       log_messages = TRUE)
  }

  rownames(exprs_) <- rownames(feature_data)
  colnames(exprs_) <- rownames(pheno_data)

  return(list(exprs = exprs_, 
              pheno_data = pheno_data, 
              feature_data = feature_data))
}

# Helper function to search for mass and retention time column names
.find_mz_rt_cols <- function(feature_data) {
  # Find mass and retention time columns
  mz_tags <- c("mass", "m.?z$", "molecular.?weight")
  rt_tags <- c("retention.?time", "^rt$", "(?=.*rt)(?=.*min)")
  mz_col <- NULL
  for (tag in mz_tags) {
    hits <- grepl(tag, tolower(colnames(feature_data)), perl = TRUE)
    if (any(hits)) {
      mz_col <- colnames(feature_data)[which(hits)[1]]
      break
    }
  }
  rt_col <- NULL
  for (tag in rt_tags) {
    hits <- grepl(tag, tolower(colnames(feature_data)), perl = TRUE)
    if (any(hits)) {
      rt_col <- colnames(feature_data)[which(hits)[1]]
      break
    }
  }
  # If not found, throw error
  if (is.null(mz_col)) {
    stop(paste0("No mass to charge ratio column found - should match one of:\n",
                paste(mz_tags, collapse = ", "), " (not case-sensitive)"))
  }
  if (is.null(rt_col)) {
    stop(paste0("No retention time column found - should match one of:\n",
                paste(rt_tags, collapse = ", "), " (not case-sensitive)"))
  }
  log_text(paste0("Identified m/z column ", mz_col,
                  " and retention time column ", rt_col))
                  
  return(list(mz_col = mz_col, rt_col = rt_col))
}

# Combines mode name, mass and retention time to create a Feature ID
.name_features <- function(feature_data) {
  cols <- .find_mz_rt_cols(feature_data)
  mz_col <- cols$mz_col
  rt_col <- cols$rt_col
  log_text(paste0("Identified m/z column ", mz_col, 
                  " and retention time column ", rt_col))
  log_text("Creating feature IDs from Split, m/z and retention time")

  round_mz <- as.numeric(feature_data[, mz_col]) %>%
    as.character() %>%
    gsub("[.]", "_", .)
  round_rt <- as.numeric(feature_data[, rt_col]) %>%
    as.character() %>%
    gsub("[.]", "_", .)
  feature_data$Feature_ID <- paste0(feature_data$Split, "_", 
                                    round_mz, "a", round_rt)
  if (anyDuplicated(feature_data$Feature_ID)) {
    duplicates <- paste0(
      feature_data$Feature_ID[duplicated(feature_data$Feature_ID)],
      collapse = ", ")
    stop("Could not create unique feature names from m/z and", 
         " retention time columns. Duplicated values: ", duplicates)
  }

  feature_data
}

#' An S4 class used to represent LC-MS datasets
#'
#' MetaboSet is the main class used to represent data in the notame package.
#' It is built upon the \code{\link[Biobase]{ExpressionSet}} class from the 
#' Biobase package. In addition to the slots inherited from 
#' \code{\link[Biobase]{ExpressionSet}}, \code{MetaboSet} has three slots of 
#' its own. The extra three slots hold special column names that are stored 
#' purely for convenience, as many functions use these as defaults.
#'
#' @slot group_col character, name of the column holding group information
#' @slot time_col character, name of the column holding time points
#' @slot subject_col character, name of the column holding subject identifiers
#'
#' @param object A MetaboSet object.
#'
#' @section Constructor:
#' See \code{\link{construct_metabosets}} for constructor function.
MetaboSet <- setClass("MetaboSet",
  slots = c(
    group_col = "character",
    time_col = "character",
    subject_col = "character"),
  contains = "ExpressionSet"
)

setValidity(
  "MetaboSet",
  function(object) {
    if (!is.na(group_col(object)) &
      !group_col(object) %in% colnames(pData(object))) {
        return(paste0("Column '", group_col(object),
                      "' not found in pheno data"))
    } else if (!is.na(time_col(object)) &
      !time_col(object) %in% colnames(pData(object))) {
        return(paste("Column", time_col(object), "not found in pheno data"))
    } else if (!is.na(subject_col(object)) &
      !subject_col(object) %in% colnames(pData(object))) {
        return(paste("Column", subject_col(object), "not found in pheno data"))
    } else if (!all(
      c("Injection_order", "Sample_ID", "QC") %in% colnames(pData(object)))) {
        return(paste0("Pheno data should contain columns", 
                      " Sample_ID, QC and Injection_order"))
    } else if (any(is.na(pData(object)[, "QC"]))) {
      return("QC column should not contain NAs")
    } else if (!"Flag" %in% colnames(fData(object))) {
      return("Flag column not found in fData")
    } else if (!identical(object$Sample_ID, rownames(pData(object)))) {
      return("Sample_ID does not match rownames in pheno data")
    } else if (!identical(fData(object)$Feature_ID, rownames(fData(object)))) {
      return("Feature_ID does not match rownames in feature data")
    } else {
      x <- .check_pheno_data(pData(object), id_prefix = "")
      x <- .check_exprs(exprs(object))
      x <- .check_feature_data(fData(object), check_limits = FALSE)
      TRUE
    }
  }
)

#' Construct MetaboSet objects
#'
#' Construct MetaboSet objects from input read by read_from_excel.
#' Returns a list of MetaboSet objects, one per mode. The modes are separated 
#' by the "Split" column in feature_data.
#'
#' @param exprs matrix, the feature abundances, size features x samples
#' @param pheno_data data frame, sample information, size sample info variables 
#' x samples
#' @param feature_data data frame, information about the features, size 
#' features x feature info columns
#' @param group_col character, the name of the column in pheno_data to use as 
#' the grouping variable
#' @param time_col character, the name of the column in pheno_data to use as 
#' the time variable
#' @param subject_col character, the name of the column in pheno_data to use as 
#' the subject ID variable
#' @param split_data logical, whether to split data by analytical mode recorded 
#' in the "Split" column of feature data. If TRUE (the default), will return a 
#' list of MetaboSet objects, one per analytical mode. If FALSE, will return a 
#' single MetaboSet object.
#'
#' @examples
#' data <- read_from_excel(
#'   file = system.file("extdata", "sample_data_whole.xlsx", 
#'   package = "notame"), sheet = 1, corner_row = 4, corner_column = "X",
#'   split_by = c("Column", "Ion mode"))
#' 
#' modes <- construct_metabosets(exprs = data$exprs, 
#'   pheno_data = data$pheno_data, feature_data = data$feature_data,
#'   group_col = "Group")
#'
#' @return A list of MetaboSet objects or a single MetaboSet object.
#'
#' @seealso \code{\link{read_from_excel}}
#'
#' @export
construct_metabosets <- function(exprs, pheno_data, feature_data,
                                 group_col = NA_character_, 
                                 time_col = NA_character_,
                                 subject_col = NA_character_,
                                 split_data = TRUE) {
  if (!"Flag" %in% colnames(feature_data)) {
    message("Initializing the object(s) with unflagged features.")
    feature_data$Flag <- NA
  }
  feature_data <- .check_feature_data(feature_data, check_limits = FALSE,
                                     log_messages = TRUE)
  exprs <- .check_exprs(exprs, log_messages = TRUE)
  log_text(paste0("Setting row and column names of exprs", 
                  " based on feature and pheno data"))
  rownames(exprs) <- rownames(feature_data)
  colnames(exprs) <- rownames(pheno_data)
  pheno_data <- Biobase::AnnotatedDataFrame(data = pheno_data)

  if (split_data) {
    # Split the data by the Split column of feature data
    parts <- unique(feature_data$Split)
    obj_list <- list()
    for (part in parts) {
      fd_tmp <- Biobase::AnnotatedDataFrame(
        data = feature_data[feature_data$Split == part, ])
      ad_tmp <- exprs[fd_tmp$Feature_ID, ]
      obj_list[[part]] <- MetaboSet(exprs = ad_tmp, phenoData = pheno_data,
                                    featureData = fd_tmp, group_col = group_col,
                                    time_col = time_col,
                                    subject_col = subject_col)
    }
    return(obj_list)
  } else {
    fd_tmp <- Biobase::AnnotatedDataFrame(feature_data)
    object <- MetaboSet(exprs = exprs, phenoData = pheno_data,
                        featureData = fd_tmp, group_col = group_col,
                        time_col = time_col, subject_col = subject_col)

    return(object)
  }
}


#' Write results to Excel file
#'
#' Writes all the data in a MetaboSet object to an Excel spreadsheet.
#' The format is similar to the one used to read data in, except for the fact 
#' that EVERYTHING NEEDS TO BE WRITTEN AS TEXT. To fix numeric values in Excel,
#' choose any cell with a number, press Ctrl + A, then go to the dropdown menu
#' in upper left corner and choose "Convert to Number". This will fix the file,
#' but can take quite a while.
#'
#' @param object a MetaboSet object
#' @param file path to the file to write
#' @param ... Additional parameters passed to \code{openxlsx::write.xlsx}
#'
#' @return None, the function is invoked for its side effect.
#'
#' @examples
#' \dontshow{.old_wd <- setwd(tempdir())}
#' write_to_excel(example_set, file = "example_set.xlsx")
#' \dontshow{setwd(.old_wd)}
#'
#' @export
write_to_excel <- function(object, file, ...) {
  # Bottom part consists of (from left to right):
  # - feature data with results
  # - abundance values
  bottom <- cbind(fData(object), exprs(object))
  # All columns must be characters to allow combination with the top block
  bottom <- bottom %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
    rbind(colnames(.), .)
  # Top block holds the sample information
  pd <- pData(object)
  datafile_cols <- colnames(pd)[grepl("Datafile", colnames(pd))]
  if (length(datafile_cols)) {
    last_datafile <- datafile_cols[length(datafile_cols)]
    pd <- pd[, c(setdiff(colnames(pd), last_datafile), last_datafile)]
    log_text(paste0("Moved ", last_datafile, " column to last",
                    " to get meaningful column names for abundances"))
  }
  top <- cbind(matrix(colnames(pd), ncol = 1), t(pd))
  # NA blocks to fill the empty space
  empty1 <- matrix(NA_character_,
                   nrow = nrow(top),
                   ncol = ncol(fData(object)) - 1)
  top <- cbind(empty1, top)
  colnames(top) <- colnames(bottom)
  # Replace exprs column names with the last column (now row) of sample info
  replace_idx <- (ncol(fData(object)) + 1):ncol(bottom)
  bottom[1, replace_idx] <- top[nrow(top), replace_idx]
  # All combined
  big <- rbind(top[seq_len(nrow(top) - 1), ], bottom)

  openxlsx::write.xlsx(big, file = file, colNames = FALSE, ...)
}


# ------------ Methods -----------------

# Print and show methods
.print_levels <- function(v) {
  t_groups <- table(v)
  if (is.factor(v)) {
    groups <- levels(v)
  } else {
    groups <- unique(v)
  }
  output <- vapply(groups, function(y) {
    obs <- t_groups[y]
    paste0(y, ": ", obs)
  }, character(1))
  output <- paste(output, collapse = ", ")
  message("  ", output)
}

setMethod("show", c(object = "MetaboSet"), 
  function(object) {
    cat(paste("MetaboSet object with", nrow(object), "features and",
              ncol(object), "samples.\n"))
    cat(paste(sum(object$QC == "QC"), "QC samples included\n"))
    cat(paste(sum(is.na(flag(object))), "non-flagged features,",
              sum(!is.na(flag(object))), "flagged features.\n\n"))
    if (!is.na(group_col(object))) {
      cat(paste0(group_col(object), ":\n"))
      .print_levels(pData(object)[, group_col(object)])
    }
    if (!is.na(time_col(object))) {
      cat(paste0(time_col(object)), ":\n")
      .print_levels(pData(object)[, time_col(object)])
    }
    if (!is.na(subject_col(object))) {
      cat(paste0(subject_col(object)), ":\n")
      subject <- as.character(pData(object)[, subject_col(object)])
      subject <- subject[!grepl("QC", subject)]
      cat(paste0("  ", length(unique(subject)), " distinct subjects\n  min:",
                 min(table(subject)), ", max:", max(table(subject)),
                 " observations per subject.\n"))
    }

    cat("\nThe object has the following parts (splits):\n")
    splits <- unique(fData(object)$Split)
    t_splits <- table(fData(object)$Split)
    for (split in splits) {
      cat(paste0("  ", split, ": ", t_splits[split], " features\n"))
    }
  }
)
#' Retrieve both sample information and features
#'
#' @param object a MetaboSet object
#'
#' @return A data frame with sample information plus all features as columns,
#' one row per sample.
#'
#' @examples
#' combined_data(example_set)
#'
#' @export
setGeneric("combined_data", signature = "object",
           function(object) standardGeneric("combined_data"))

#' @describeIn MetaboSet sample information and features combined to a single 
#' data frame, one row per sample.
#' @export
setMethod("combined_data", c(object = "MetaboSet"), 
  function(object) {
    cbind(pData(object), t(exprs(object)))
  }
)


#' Get and set name of the special column for group labels
#' @param object a MetaboSet object
#' 
#' @return Character, the name of the grouping variable.
#'
#' @examples
#' # Get name of grouping variable
#' group_col(example_set)
#'
#' @export
setGeneric("group_col", signature = "object",
           function(object) standardGeneric("group_col"))
#' @describeIn MetaboSet get name of column for group labels
#'
#' @export
setMethod("group_col", "MetaboSet", function(object) object@group_col)

#' @rdname group_col
#' @param object a MetaboSet object
#' @param value string, name of column to be designated for holding group labels
#' 
#' @return For the endomorphism, an object with the grouping variable set to 
#' the specified variable.
#'
#' @examples
#' # Set grouping variable
#' group_col(example_set) <- "Group"
#'
#' @export
setGeneric("group_col<-", signature = "object",
           function(object, value) standardGeneric("group_col<-"))

#' @describeIn MetaboSet set name of column for group labels
#' @param value string, name of column to be designated for holding group labels
#' @export
setMethod("group_col<-", "MetaboSet", 
  function(object, value) {
    object@group_col <- value
    if (validObject(object)) {
      return(object)
    }
  }
)
#' Get and set the name of the special column for time points
#' @param object a MetaboSet object
#'
#' @return Character, name of time variable.
#'
#' @examples
#' # Get name of time variable
#' time_col(example_set)
#'
#' @export
setGeneric("time_col", signature = "object",
           function(object) standardGeneric("time_col")
)

#' @describeIn MetaboSet get name of column for time points
#' @export
setMethod("time_col", "MetaboSet",
          function(object) object@time_col)

#' @rdname time_col
#' @param object a MetaboSet object
#' @param value string, name of column to be designated for holding time points
#'
#' @return For the endomorphism, an object with the time variable set to the 
#' specified variable.
#'
#' @examples
#' # Set time variable
#' time_col(example_set) <- "Time"
#'
#' @export
setGeneric("time_col<-", signature = "object",
           function(object, value) standardGeneric("time_col<-"))

#' @describeIn MetaboSet set name of column for time points
#' @export
setMethod("time_col<-", "MetaboSet",
  function(object, value) {
    object@time_col <- value
    if (validObject(object)) {
      return(object)
    }
  }
)


#' Get and set the name of special column for subject identifiers
#' @param object a MetaboSet object
#'
#' @return Character, the name of the subject variable.
#'
#' @examples
#' # Get name of subject variable
#' subject_col(example_set)
#' @export
setGeneric("subject_col", signature = "object",
           function(object) standardGeneric("subject_col"))

#' @describeIn MetaboSet get name of column for subject subject identifiers
#' @export
setMethod("subject_col", "MetaboSet",
          function(object) object@subject_col)

#' @rdname subject_col
#' @param object a MetaboSet object
#' @param value string, name of column to be designated for holding subject 
#' identifiers
#'
#' @return For the endomorphism, an object with the subject variable set to the 
#' specified variable.
#'
#' @examples
#' # Set subject variable
#' subject_col(example_set) <- "Subject_ID"
#' @export
setGeneric("subject_col<-", signature = "object",
           function(object, value) standardGeneric("subject_col<-"))

#' @describeIn MetaboSet set name of column for subject identifiers
#' @export
setMethod("subject_col<-", "MetaboSet",
  function(object, value) {
    object@subject_col <- value
    if (validObject(object)) {
      return(object)
    }
  }
)



#' Get and set the values in the flag column
#' @param object a MetaboSet object
#' 
#' @return Character vector of feature flags.
#'
#' @examples
#' # Get values in flag column of fData
#' flag(example_set)
#'
#' @export
setGeneric("flag", signature = "object",
           function(object) standardGeneric("flag"))

#' @describeIn MetaboSet get flags
#' @export
setMethod("flag", "MetaboSet",
          function(object) fData(object)$Flag)

#' @rdname flag
#' @param object a MetaboSet object
#' @param value character vector, values for flag column
#'
#' @return For the endomorphism, an object with a modified flag column.
#'
#' @examples
#' # Flag a suspicious feature manually
#' flag(example_set)[1] <- "Contaminant, known from experience"
#' @export
setGeneric("flag<-", signature = "object",
           function(object, value) standardGeneric("flag<-"))

#' @describeIn MetaboSet set flags
#' @export
setMethod("flag<-", "MetaboSet",
  function(object, value) {
    fData(object)$Flag <- value
    if (validObject(object)) {
      return(object)
    }
  }
)


#' Join new columns to feature data
#'
#' Join a new data frame of information to feature data of a MetaboSet object.
#' The data frame needs to have a column "Feature_ID".
#' This function is usually used internally by some of the functions in the 
#' package, but can be useful.
#'
#' @param object a MetaboSet object
#' @param dframe a data frame with the new information
#'
#' @examples
#' new_info <- data.frame(
#'   Feature_ID = featureNames(example_set),
#'   Feature_number = seq_len(nrow(example_set))
#' )
#' with_new_info <- join_fData(example_set, new_info)
#' colnames(fData(with_new_info))
#'
#' @return A MetaboSet object with the new information added to fData(object).
#'
#' @export
setGeneric("join_fData", signature = c("object", "dframe"),
           function(object, dframe) standardGeneric("join_fData"))

#' @describeIn MetaboSet join new information to fData
#' @param dframe a data frame with the new information
#' @export
setMethod("join_fData", c("MetaboSet", "data.frame"),
  function(object, dframe) {
    fData(object) <- dplyr::left_join(fData(object), dframe, by = "Feature_ID")
    rownames(fData(object)) <- fData(object)$Feature_ID
    if (validObject(object)) {
      return(object)
    }
  }
)

#' Join new columns to pheno data
#'
#' Join a new data frame of information to pheno data of a MetaboSet object.
#'
#' @param object a MetaboSet object
#' @param dframe a data frame with the new information
#'
#' @examples
#' new_info <- data.frame(
#'   Sample_ID = sampleNames(example_set),
#'   BMI = stats::runif(ncol(example_set), 22, 26)
#' )
#' with_new_info <- join_pData(example_set, new_info)
#' colnames(pData(with_new_info))
#'
#' @return A MetaboSet object with the new information added to pData(object).
#'
#' @export
setGeneric("join_pData", signature = c("object", "dframe"),
           function(object, dframe) standardGeneric("join_pData"))

#' @describeIn MetaboSet join new information to fData
#' @param dframe a data frame with the new information
#' @export
setMethod("join_pData", c("MetaboSet", "data.frame"),
  function(object, dframe) {
    pData(object) <- dplyr::left_join(pData(object), dframe)
    rownames(pData(object)) <- pData(object)$Sample_ID
    if (validObject(object)) {
      return(object)
    }
  }
)

# FeatureNames also changing Feature_ID column in featureData
setMethod("featureNames<-",
  signature = signature(object = "MetaboSet", value = "ANY"),
  function(object, value) {
    fd <- featureData(object)
    featureNames(fd) <- value
    ad <- assayData(object)
    featureNames(ad) <- value
    object@featureData <- fd
    object@assayData <- ad
    fData(object)$Feature_ID <- value
    if (validObject(object)) {
      return(object)
    }
  }
)


setMethod("sampleNames<-",
  signature = signature(object = "MetaboSet", value = "ANY"),
  function(object, value) {
    pd <- phenoData(object)
    sampleNames(pd) <- value
    ad <- assayData(object)
    sampleNames(ad) <- value
    prd <- protocolData(object)
    if (nrow(prd) == 0) {
      prd <- pd[, integer(0)]
    } else {
      sampleNames(prd) <- value
    }
    object@phenoData <- pd
    object@protocolData <- prd
    object@assayData <- ad
    pData(object)$Sample_ID <- value
    if (validObject(object)) {
      return(object)
    }
  }
)
