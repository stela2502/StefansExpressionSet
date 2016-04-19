#' @name write.data
#' @aliases write.data,StefansExpressionSet-method
#' @rdname write.data-methods
#' @docType methods
#' @description  write the StefansExpressionSet data table to disk
#' @param x the StefansExpressionSet object
#' @param annotation a vector of annotation data column names to include in the written table (default=none)
#' @title description of function write.data
#' @export 
setGeneric('write.data', ## Name
	function ( x, annotation=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('write.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('write.data', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, annotation=NULL ) {
	if ( !is.null(annotation) ) {
		write.table( cbind( x@annotation[,annotation], x@data), file= paste( x@outpath,x@name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
	else {
		write.table( cbind( rownames(x@data), x@data), file= paste( x@outpath,x@name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
})

