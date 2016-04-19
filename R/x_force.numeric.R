#' @name force.numeric
#' @aliases force.numeric,StefansExpressionSet-method
#' @rdname force.numeric-methods
#' @docType methods
#' @description  The moethod forces the values in the data matrix to be numbers.
#' @param dataObj the StefansExpressionSet object
#' @title description of function force.numeric
#' @export 
setGeneric('force.numeric', ## Name
	function (dataObj ) { ## Argumente der generischen Funktion
		standardGeneric('force.numeric') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('force.numeric', signature = c ('StefansExpressionSet') ,
	definition = function ( dataObj ) {
	for ( i in 1: ncol(dataObj@data) ) { 
		if ( !  paste(as.vector(dataObj@data[,i]), collapse=" ") == paste(as.vector(as.numeric(dataObj@data[,i])), collapse=" ") ) { 
			dataObj@data[,i] <- as.vector(as.numeric(dataObj@data[,i]))
		}
	}
	dataObj
})

