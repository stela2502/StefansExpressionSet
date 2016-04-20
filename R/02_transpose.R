#' @name t
#' @aliases t,02_transpose-method
#' @rdname t-methods
#' @docType methods
#' @description transpose the object DELETING all stats and usedObj!
#' @param x the StefansExpressionSet object
#' @title description of function t
#' @export 
setMethod('t', signature = c ('StefansExpressionSet'),
		definition = function (x) {
			if ( is.null(x@usedObj[['transposed']])){
				x@usedObj[['transposed']] = FALSE
			}
			x@usedObj[['transposed']]
			ret <- new( class(x)[1] , data = data.frame(t(x@data)), 
					samples = x@annotation, 
					annotation = x@samples, 
					outpath = x@outpath,
					sampleNamesCol = x@rownamescol,
					rownamescol = x@sampleNamesCol,
					snorm= x@snorm,
					zscored= x@zscored
			)
			ret@usedObj[['transposed']] = ! x@usedObj[['transposed']]
			ret
		} 
)
