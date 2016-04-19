#' @name interesting_genes
#' @aliases interesting_genes,SingleCellsNGS-method
#' @rdname interesting_genes-methods
#' @docType methods
#' @description This method uses the values created by update.measurements() and selects the most interesting subset of genes.
#' @param x the SingleCellsNGS object
#' @param Min the minimal fraction of cells expresing the gene default=0.4
#' @param Max the maximal fraction of cells expresing the gene default=0.8 
#' @param genes how many genes you want to get back default=100
#' @title description of function interesting_genes
#' @export 
setGeneric('interesting_genes', ## Name
		function (x, Min=0.4, Max=0.8 , genes=100) { ## Argumente der generischen Funktion
			standardGeneric('interesting_genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('interesting_genes', signature = c ('SingleCellsNGS'),
		definition = function (x, Min=0.4, Max=0.8 , genes=100) {
			if ( is.null(x@annotation$mean_var) ) {
				x <- update.measurements ( x )
			}
			possible <- which( x@annotation$fraction_expressing > Min & x@annotation$fraction_expressing < Max )
			order.possible <- order(x@annotation$mean_var[possible])
			rownames(x@data)[possible[order.possible][1:genes]]
		} 
)




