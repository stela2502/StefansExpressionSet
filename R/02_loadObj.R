#' @name loadObj
#' @aliases loadObj,StefansExpressionSet-method
#' @rdname loadObj-methods
#' @docType methods
#' @description This function loads the a StefansExpressionSet file and returns the data R object int them.
#'  The function is mainly to keep the scripts simpler. There is no magic involved.
#' @param file the infile
#' @title description of function loadObj
#' @export 
setGeneric('loadObj', ## Name
		function ( file=NULL ){	## Argumente der generischen Funktion
			standardGeneric('loadObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('loadObj', signature = c ('StefansExpressionSet'),
		definition = function ( file=NULL ){
			if ( is.null(file)){
				stop( "Sorry I need to filename to load from" )
			}
			load( file.path(data@outpath, file) )
			data
		}
)


