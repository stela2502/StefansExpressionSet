#' @name sd.filter
#' @aliases sd.filter,stat_sd.filter-method
#' @rdname sd.filter-methods
#' @docType methods
#' @description drops all genes that have an sd of 0 over the whole dataset
#' @param x  TEXT MISSING
#' @title description of function sd.filter
#' @export 
setGeneric('sd.filter', ## Name
		function ( x ) { ## Argumente der generischen Funktion
			standardGeneric('sd.filter') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('sd.filter', signature = c ('StefansExpressionSet'),
		definition = function ( x ) {
			reduce.Obj( x, names(which(apply(x@data,1,sd) != 0)), name=x@name)
		}
)

