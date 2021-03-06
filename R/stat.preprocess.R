#' @name preprocess
#' @aliases preprocess,NGSexpressionSet-method
#' @rdname preprocess-methods
#' @docType methods
#' @description create the count dataset and normalize it following DEseq standard
#' @description  function used by DEseq Do not use
#' @param x the NGSexpressionSet
#' @param condition the column in the samples table to use as grouing variable
#' @title description of function preprocess
#' @export
setGeneric('preprocess', ## Name
		function (x, condition ) { ## Argumente der generischen Funktion
			standardGeneric('preprocess') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('preprocess', signature = c ('NGSexpressionSet'),
		definition = function (x, condition) {
			if ( is.null( x@usedObj[['cds']]) ) {
				x@usedObj[['cds']] <- list()
			}
			if ( is.na( match( condition, names(x@usedObj[['cds']]) ) )) {
				#condition <- as.factor(x@samples[,condition])
				dat <- list()
				if ( nrow(x@raw) > 0 ) {
					t <- as.matrix(x@raw)
				}else {
					t <- as.matrix(x@data)
				}
				#colnames( t ) <- x@samples[,condition]
				colnames(t) <- NULL
				dat$cds <- DESeq2::DESeqDataSetFromMatrix(
					t,
					x@samples,
					formula( paste('~',condition ) )
				)
				#dat$cds <- DESeq2::estimateSizeFactors(dat$cds)
				#DESeq2::sizeFactors(x@cds)
				#dat$cds <- DESeq2::estimateDispersions(dat$cds)
				#dat$vsdFull = DESeq2::varianceStabilizingTransformation( dat$cds )
				x@usedObj[['cds']][[length(x@usedObj[['cds']])+1]] = dat$cds
				names(x@usedObj[['cds']])[length(x@usedObj[['cds']])] = condition
			}
			x
		}
)
