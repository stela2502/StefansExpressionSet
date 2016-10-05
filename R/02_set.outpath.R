#' @name set.outpath
#' @aliases set.outpath,StefansExpressionSet-method
#' @rdname set.outpath-methods
#' @docType methods
#' @description Drop the samples, that have been selected!
#' @param x the StefansExpressionSet object
#' @param path the new outpath (defaults to working path)
#' @title description of function set.outpath
#' @export 
setGeneric('set.outpath', ## Name
		function ( x, path=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('set.outpath') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('set.outpath', signature = c ( 'StefansExpressionSet') ,
		definition = function ( x, path=NULL ) {
			if ( is.null(path) ){
				path = pwd()
			}
			x@outpath = path
			change_path <- function(name) {
				
			}
			if ( !is.null(x@usedObj$rfObj_row) ) {
				OPATH <- paste( x@outpath,"/",str_replace( x@name, '\\s', '_'), sep='')
				opath = paste( OPATH,"/RFclust.mp",sep='' )
				x@usedObj$rfObj_row <- lapply( x@usedObj$rfObj_row, function( RFobj ) { RFobj@tmp.path <- opath
							RFobj})
			}
			if ( !is.null(x@usedObj$rfObj) ) {
				OPATH <- paste( x@outpath,"/",str_replace( x@name, '\\s', '_'), sep='')
				opath = paste( OPATH,"/RFclust.mp",sep='' )
				x@usedObj$rfObj <- lapply( x@usedObj$rfObj, function( RFobj ) { RFobj@tmp.path <- opath 
							RFobj})
			}
			x
		})