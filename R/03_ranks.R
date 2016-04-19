#' @name ranks
#' @aliases ranks,StefansExpressionSet-method
#' @rdname ranks-methods
#' @docType methods
#' @description  creates a new ranks vector in the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @title description of function ranks
#' @export 
setGeneric('ranks', ## Name
	function (x ) { ## Argumente der generischen Funktion
		standardGeneric('ranks') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('ranks', signature = c ('StefansExpressionSet') ,
	definition = function (x ) {
	if ( is.na(x@usedObj[['ranks']]) ){
		x@usedObj[['ranks']] <- apply( x@data,2,order)
		colnames( x@usedObj[['ranks']] ) <- colnames(x@data) 
		rownames( x@usedObj[['ranks']] ) <- rownames(x@data) 
	}
	x
})
