#' @name addAnnotation
#' @aliases addAnnotation,StefansExpressionSet-method
#' @rdname addAnnotation-methods
#' @docType methods
#' @description  add aditional annotation to the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @param mart the annotation table (data.frame or mart object)
#' @param mart.col which column corresponds to the rownames(x@data)
#' @title description of function addAnnotation
#' @export 
setGeneric('addAnnotation', ## Name
	function (x ,mart, mart.col='refseq_mrna') { ## Argumente der generischen Funktion
		standardGeneric('addAnnotation') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('addAnnotation', signature = c ( 'StefansExpressionSet') ,
	definition = function (x ,mart, mart.col='refseq_mrna') {
	if ( ! class(mart) == 'data.frame' ){
		x@annotation <- cbind(x@annotation, mart[match(rownames(x@data),mart[,mart.col] ), ] )
	}
	else {
		x@annotation <- mart[is.na(match(mart[,mart.col], rownames(x@data)))==F,]
		rownames(x@annotation) <- rownames(x@data)
	}
	x
})

