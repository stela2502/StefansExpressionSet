#' @name quality_of_fit
#' @aliases quality_of_fit,Clusters-method
#' @rdname quality_of_fit-methods
#' @docType methods
#' @description this calculated a quality of fit for the expression data
#' @param x the StefansExpressionSet object
#' @param what where to group the data on (cells or genes; default='cells')
#' @param col which samples/annoation column to use for the grouping (default='groupID')
#' @title description of function quality_of_fit
#' @export 
setGeneric('quality_of_fit', ## Name
	function ( x, what='cells', col='groupID' ) { ## Argumente der generischen Funktion
		standardGeneric('quality_of_fit') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('quality_of_fit', signature = c ('StefansExpressionSet'),
	definition = function ( x, what='cells', col='groupID' ) {
	
	test <- x@data
	test[which(test ==  0 ) ] = NA
	if ( what=='cells') {
		clusters <- x@samples[,col]
		ret <- list ( 'single' = apply(test,2, difference, as.character(clusters), max(clusters) ) )
		ret$sum = round(sum(ret$single))
	}
	else if ( what=='genes') {
		clusters <- x@annotation[,col]
		ret <- list ( 'single' = apply(test,1, difference, as.charcter(clusters), max(clusters) ) )
		ret$sum = round(sum(ret$single))
	}
	else {
		stop(paste( what,': what option is not supported!') )
	}
	ret
} )
