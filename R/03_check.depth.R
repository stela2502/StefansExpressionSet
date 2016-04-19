#' @name check.depth
#' @aliases check.depth,NGSexpressionSet-method
#' @rdname check.depth-methods
#' @docType methods
#' @description this function checks how many reads are taken from the total dataset using the top 5 percent of genes only
#' @param x the NGSexpressionSet
#' @param gene an optional list of genes you want to select for
#' @param cutoff a cutoff that defines how many reads may be taken from the top 5 percent genes before a sample gets rejected
#' @title description of function check
#' @return a vector of sample names that did not pass the check
#' @export 
setGeneric('check.depth', ## Name
	function (x, genes=NULL, cutoff = 0.77 ) { ## Argumente der generischen Funktion
		standardGeneric('check.depth') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('check.depth', signature = c ('NGSexpressionSet'),
	definition = function (x, genes=NULL, cutoff = 0.77 ) {
	percent5 <-  reads.taken(x, 0.05, genes)
	names(which(percent5$reads.taken > cutoff )) ## obtained experimentally using Jennies HSC dataset
})

