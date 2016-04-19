#' @name difference
#' @aliases difference,Clusters-method
#' @rdname difference-methods
#' @docType methods
#' @description this code has to be revised! NOT WORKING
#' @param x the StefansExpressionSet object
#' @param clusters the grouoing strings
#' @param groups.n the max number of groups
#' @title description of function difference
setGeneric('difference', ## Name
	function ( x, clusters, groups.n ) { ## Argumente der generischen Funktion
		standardGeneric('difference') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('difference', signature = c ('data.frame'),
	definition = function ( x, clusters, groups.n ) {
	ret = 0 
	for ( i in 1:groups.n  ) {
		a <- x[which( clusters == i)]
		a <- a[- (is.na(a))==F]
		if ( length(a) > 1 ) {  ret = ret + sum( (a- mean(a) )^2 ) }
	}
	ret
} )
