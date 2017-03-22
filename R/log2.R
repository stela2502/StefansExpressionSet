#' @name log2
#' @aliases log2,StefansExpressionSet-method
#' @rdname log2-methods
#' @docType methods
#' @description Converts the x@data values to log2(x@data+1) values if not already done
#' @param x  the StefansExpressionSet
#' @title description of function log2
#' @export
setMethod('log2', signature = c ('StefansExpressionSet'),
	definition = function (x) {
		if ( is.null(x@usedObj$log2)){
			x@usedObj$log2 = 0
		}
	if ( ! x@usedObj$log2 ) {
		x@usedObj$log2 = 1
		x@data <- as.data.frame(apply(x@data,2,function(x) { log2(x+1) } ))
	}
	x
} )
