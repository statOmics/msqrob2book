createPairwiseContrasts <- function(formula, coldata, var, ridge = FALSE, nullHypothesis = " = 0") {
    params <- .getParamNames(formula, coldata, var, ridge)
    out <- combn(params, 2)
    out <- paste0(out[2, ], " - ", out[1, ])
    hasIntercept <- attr(terms(formula), "intercept")
    isFirst <- labels(terms(formula))[[1]] == var
    if (hasIntercept && isFirst)
        out <- c(out, params)
    paste0(out, nullHypothesis)
}

createOneVersusAllContrasts <- function(formula, coldata, var,
                                        ridge = FALSE,
                                        nullHypothesis = " = 0") {
    params <- .getParamNames(formula, coldata, var, ridge)
    averages <- sapply(seq_along(params), function(i) {
        paramSum <- paste(params[-i], collapse = " + ")
        paste0("(", paramSum, ") / ", length(params) - 1)
    })
    paste0(params, " - ", averages, nullHypothesis)
}

.getParamNames <- function(formula, coldata, var, ridge) {
    params <- lme4::nobars(formula) |>
        model.matrix(, data = coldata) |>
        colnames()
    params <- params[grepl(var, params)]
    if (length(params) == 0)
        stop(sQuote(var), " not found as fixed effect in formula.")
    if (any(grepl("[:]", params)))
        stop("You cannot automatically create pairwise contrasts for a variable involved in an interaction.")
    if (ridge) params <- paste0("ridge", params)
    unique(params)
}

.makeResultColumns <- function(contrast, resultsColumnNamePrefix) {
    if (is.null(colnames(contrast))) {
        if (resultsColumnNamePrefix == "")
            resultsColumnNamePrefix <- "msqrobResults"
        if (ncol(contrast) > 1)
            colnames(contrast) <- seq_len(ncol(contrast))
    }
    paste0(resultsColumnNamePrefix, colnames(contrast))
}

.getMsqrob2Results <- function(object, row_var) {
    if (!row_var %in% colnames(rowData(object)))
        stop(
            sQuote(row_var), "is not found in the rowData of the ",
            "SummarizedExperiment object. Make sure you set the ",
            "same 'contrast' and 'resultsColumnNamePrefix' as for ",
            "'hypothesisTest()'."
        )
    rowData(object)[[row_var]]
}

##' @param object A `SummarizedExperiment` object containing model
##'    inference performed by `hypothesisTest()`.
##' @param contrast A numeric `matrix` providing the contrast matrix
##'     used during model estimation.
##' @param resultsColumnNamePrefix A `character(1)` providing the
##'     prefix used when running `hypothesisTest()`.
##' @param combine A `logical(1)` indicating whether the result tables
##'     should be combined in a single table (default) or return as a
##'     list of tables. When combined, two new variables are created:
##'     1. `contrast`tells from which contrast (column of the `contrast`
##'     matrix) the results where obtained; 2. `feature` provides the
##'     name of the modelled feature (this information is taken from
##'     the rownames of the tables, but this are made unique upon
##'     combining).
msqrobCollect <- function(object, contrast, resultsColumnNamePrefix = "",
                           combine = TRUE) {
    resultCols <- .makeResultColumns(contrast, resultsColumnNamePrefix)
    out <- lapply(resultCols, .getMsqrob2Results, object = object)
    names(out) <- resultCols
    if (combine) {
        for (i in names(out)) {
            out[[i]]$contrast <- i
            out[[i]]$feature <- rownames(out[[i]])
        }
        out <- do.call(rbind, out)
    }
    out
}
