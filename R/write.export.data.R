#' @name export.data
#' @aliases export.data,StefansExpressionSet-method
#' @rdname export.data-methods
#' @docType methods
#' @description  write the StefansExpressionSet data, annotation and samples tables to disk
#' @param x the StefansExpressionSet object
#' @title description of function write.data
#' @export 
setGeneric('export.data', ## Name
		function ( x ) { ## Argumente der generischen Funktion
			standardGeneric('export.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('export.data', signature = c ( 'StefansExpressionSet') ,
		definition = function ( x ) {
			write.table(cbind ( x@annotation, x@data ), file=paste(x@outpath,x@name,"_expressionValues.xls",sep=''),sep='\t', row.names=F,quote=F )
			write.table( x@samples, file=paste( x@outpath, x@name, "_sampleInformation.xls",sep='' ), sep='\t', row.names=F,quote=F )
		})

