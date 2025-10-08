setGeneric("getStatModels", function(object, ...) standardGeneric("getStatModels"))

setMethod(
    "getStatModels",
    signature = "SummarizedExperiment",
    function(object, modelColumn) {
        rd <- rowData(object)
        models <- rd[[modelColumn]]
        if (is.null(models))
            stop(
                sQuote(modelColumn), " is not a rowData column in the ",
                "proviedd SummarizedExperiment object."
            )
        if (!inherits(models[[1]], "StatModel"))
            stop(
                "The ", sQuote(modelColumn), " rowData column does ",
                "not contain msqrob2 results."
            )
        models
    }
)

setMethod(
    "getStatModels",
    signature = "QFeatures",
    function(object, i, modelColumn) {
        getStatModels(object[[i]], modelColumn)
    }
)

setGeneric("getCoef", function(object, ...) standardGeneric("getCoef"))

#' @rdname statModelAccessors
setMethod("getCoef",
          signature = "StatModel",
          definition = function(object) object@params$coefficients
)


#' @rdname statModelAccessors
setMethod(
    "getCoef",
    signature = "SummarizedExperiment",
    definition = function(object, modelColumn, j) {
        models <- getStatModels(object, modelColumn)
        getCoef(models[[j]])
    }
)

#' @rdname statModelAccessors
setMethod(
    "getCoef",
    signature = "QFeatures",
    definition = function(object, i, modelColumn, j) {
        getCoef(object[[i]], modelColumn, j)
    }
)

setGeneric("getCoefNames", function(object, ...) standardGeneric("getCoefNames"))

setMethod(
    "getCoefNames",
    signature = "StatModel",
    definition = function(object, dropRandom = TRUE) {
        coefNames <- names(getCoef(object))
        if (dropRandom) {
            coefNames <- coefNames[!grepl("\\(Intercept\\).", coefNames)]
        }
        coefNames
    }
)

setMethod(
    "getCoefNames",
    signature = "SummarizedExperiment",
    definition = function(object, modelColumn, dropRandom = TRUE) {
        models <- getStatModels(object, modelColumn)
        coefNames <- sapply(models, getCoefNames, dropRandom)
        unique(unlist(coefNames))
    }
)

setMethod(
    "getCoefNames",
    signature = "QFeatures",
    definition = function(object, i, modelColumn, dropRandom = TRUE) {
        getCoefNames(object[[i]], modelColumn, dropRandom)
    }
)

detailPlot <- function(object,
                       featureName,
                       sample_id = "colname",
                       feature_id = "rowname",
                       type = "line",
                       i = NULL,
                       colourVar = NULL) {
    if (!featureName %in% unlist(rownames(object)))
        stop(
            "The featureName ", sQuote(featureName),
            " is not found in the data."
        )
    object <- object[featureName, , ]
    if (!is.null(i)) object <- object[, , i]
    df <- object |>
        longForm(colvars = colourVar) |>
        data.frame() |>
        mutate_if(is.character, as.factor) |>
        droplevels()
    dplot <- ggplot(df) +
        aes(x = colname,
            y = value,
            shape = rowname,
            group = rowname) +
        facet_grid(~ assay) +
        labs(
            title = paste0("Detail plot for ", featureName),
            x = "sample",
            y = "feature intensity"
        )
    if (type == "line") dplot <- dplot + geom_line()
    else if (type == "boxplot") dplot <- dplot + geom_boxplot()
    else stop(sQuote(type), " is not a valid type of detail plot.")
    dplot + geom_point(aes(colour = .data[[colourVar]]))
}

overdispersionPlot <- function(object, rowvar, modelColumn) {
    models <- getStatModels(se, modelColumn)
    data.frame(
        vars = sapply(models, getVar),
        x = rowData(se)[[rowvar]]
    ) |>
        ggplot() +
        aes(x = x,
            y = vars) +
        geom_point() +
        labs(x = rowvar,
             y = "Residual variance",
             title = paste0("Overdispersion plot trended by ",
                            sQuote(rowvar)))
}

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


setGeneric("lowCountIsNA", function(object, ...) standardGeneric("lowCountIsNA"))

setMethod(
    "lowCountIsNA",
    signature = "SummarizedExperiment",
    function(object, n = 2) {
        x <- assay(object)
        x[aggcounts(object) < n] <- NA
        assay(object) <- x
        x
    }
)

setMethod(
    "lowCountIsNA",
    signature = "QFeatures",
    function(object, i, n = 2) {
        object[[i]] <- lowCountIsNA(object[[i]], n)
        object
    }
)

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

library("iSEE")
OverdispersionPanel <- createCustomPlot(
    function(se, rows, cols, rowvar = ".n",
             modelColumn = "msqrobModels") {
        overdispersionPlot(se, rowvar, modelColumn) +
            theme_bw()
    }, fullName = "Overdispersion plot"
)

msqrob2iSEE <- function(object, i, resultsColumnName = NULL, ...) {
    i <- QFeatures:::.normIndex(object, i)
    se <- getWithColData(object, i)
    stopifnot(length(resultsColumnName) <= 1)
    args <- list(...)
    if (!is.null(resultsColumnName)) {
        require(iSEEu)
        results <- .getMsqrob2Results(se, resultsColumnName)
        rowData(se) <- cbind(rowData(se), results)
        args$extra <- c(args$extra, VolcanoPlot())
    }
    if (!is.null(rowData(se)$.n))
        args$extra <- c(args$extra, OverdispersionPanel())
    do.call(iSEE, c(list(se), args))
}

