#' @name colors_4
#' @aliases colors_4,StefansExpressionSet-method
#' @rdname colors_4-methods
#' @docType methods
#' @description Create the colour information for a samples or annotation column
#' @param x the StefansExpressionSet object
#' @param name the name of data column to colour
#' @param colFunc a colour function like default = function(x) {rainbow(x)}
#' @title description of function createRFgrouping_col
#' @export 
setGeneric('colors_4', ## Name
		function ( x, name,  colFunc = NULL  ) { ## Argumente der generischen Funktion
			standardGeneric('colors_4') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('colors_4', signature = c ('StefansExpressionSet'),
		definition = function ( x, name, colFunc = NULL ) {
			if ( is.null(x@usedObj[['colorRange']])){
				x@usedObj[['colorRange']] = list()
			}
			if ( is.null(colFunc) ){
				colFunc = function(x) {rainbow(x)}
			}
			
			mix <- function ( col, l ) {
				if ( l > 10 ){
					if ( l %% 2 == 0 ) {
						col <- c(rbind( col[1:(l/2)], col[(l/2+1):l] ))
					}else {
						col <- c(rbind( col[1:(l/2+1)], c(col[(l/2+1):l],'white') ))[1:l]
					}
				}
				col
			}
			if ( is.na( match( name, names(x@usedObj[['colorRange']])))) {
				if ( !is.na( match(name, colnames(x@samples)))){
					x@usedObj[['colorRange']][[name]] <- mix(
							colFunc( length(levels( x@samples[, name ]))) ,
							length(levels( x@samples[, name ])) 
					)
				}
				else if ( !is.na( match(name, colnames(x@annotation)))){
					x@usedObj[['colorRange']][[name]] <- mix( 
							colFunc( length(levels( x@annotation[, name ]))),
							length(levels( x@annotation[, name ]))
					)
				}
				else {
					stop( "Sorry this column is nether defined in the samples nor in the annotation table!" )
				}
			}
			x
		}
)