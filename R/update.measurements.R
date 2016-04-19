#' @name update.measurements
#' @aliases update.measurements,SingleCellsNGS-method
#' @rdname update.measurements-methods
#' @docType methods
#' @description This method is inspired by PMID:24531970; it calculates the expression variance over all samples
#' @description and the fraction of cells expressing for each gene. This data is used by interesting_genes() to identify
#' @description putatively interesting genes in a SingleCellNGS object.
#' @param x the SingleCellsNGS object
#' @title description of function update.measurements
#' @export 
setGeneric('update.measurements', ## Name
		function ( x ) { ## Argumente der generischen Funktion
			standardGeneric('update.measurements') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('update.measurements', signature = c ('SingleCellsNGS'),
		definition = function ( x ) {
			x@annotation$mean_var <- apply( x@data , 1, function (x ) { 
							a<-x[which(x > 0)]
							if ( length(a) < 5 ){
								0
							}
							else {
								mean(a) / var(a) 
							}
						} )
			x@annotation$fraction_expressing <-apply( x@data , 1, function (x ) { length(which(x > 0)) / length(x) })
			x
		} 
)
