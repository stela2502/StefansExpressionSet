#' @name getAnnotation4probesets
#' @aliases getAnnotation4probesets,StefansExpressionSet-method
#' @rdname getAnnotation4probesets-methods
#' @docType methods
#' @description  subsets the annotation table for a set of probesets
#' @param x the StefansExpressionSet object
#' @param probesets the vector of probesets to annotate (rownames(x@data))
#' @param colname a vector of colnames to annotate
#' @return a vector of annotation values
#' @title description of function getAnnotation4probesets
setGeneric('getAnnotation4probesets', ## Name
	function (x, probesets=c(), colname='Gene.Symbol' ) { ## Argumente der generischen Funktion
		standardGeneric('getAnnotation4probesets') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getAnnotation4probesets', signature = c ( 'StefansExpressionSet') ,
	definition = function (x, probesets=c(), colname='Gene.Symbol' ) {
	as.vector(x@annotation[match( probesets, rownames(x@data) ), colname ])
})

