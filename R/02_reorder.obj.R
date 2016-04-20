
#' @name reorder.samples
#' @aliases reorder.samples,StefansExpressionSet-method
#' @rdname reorder.samples-methods
#' @docType methods
#' @description this function reorderes the StefansExpressionSet object based on a column in the annotation table (e.g. for plotting)
#' @param dataObj the StefansExpressionSet object
#' @param column the annotation column to reorder on
#' @title description of function remove.genes
#' @export 
setGeneric('reorder.samples', ## Name
		function ( dataObj, column ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reorder.samples', signature = c ('StefansExpressionSet'),
		definition = function ( dataObj, column ) {
			dataObj@data <- dataObj@data[ , order( dataObj@samples[,column])]
			dataObj@samples <- dataObj@samples[order( dataObj@samples[,column]),]
			dataObj
		}
)

#' @name reorder.genes
#' @aliases reorder.genes,StefansExpressionSet-method
#' @rdname reorder.genes-methods
#' @docType methods
#' @description this function reorderes the StefansExpressionSet object based on a column in the samples table (e.g. for plotting)
#' @param dataObj the StefansExpressionSet object
#' @param column the samples column to reorder on
#' @title description of function remove.genes
#' @export 
setGeneric('reorder.genes', ## Name
		function ( dataObj, column ) { ## Argumente der generischen Funktion
			standardGeneric('reorder.genes') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('reorder.genes', signature = c ('StefansExpressionSet'),
		definition = function ( dataObj, column ) {
			dataObj@data <- dataObj@data[ order( dataObj@annotation[,column]),]
			dataObj@annotation <- dataObj@annotation[order( dataObj@annotation[,column]),]
			dataObj
		}
)

