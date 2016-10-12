#' @name describe.probeset
#' @aliases describe.probeset,StefansExpressionSet-method
#' @rdname describe.probeset-methods
#' @docType methods
#' @description print all annotation values + all values + all stat results for one probeset
#' @param x the StefansExpressionSet object
#' @param probeset the probeset to describe
#' @title description of function describe.probeset
#' @export 
setGeneric('describe.probeset', ## Name
		function ( x, probeset ) { ## Argumente der generischen Funktion
			standardGeneric('describe.probeset') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('describe.probeset', signature = c ('StefansExpressionSet'),
		definition = function ( x, probeset ) {
			ret <- list()
			print ( paste ("Annoataion for probeset ",probeset ) )
			ret$annotation <- x@annotation[ probeset, ] 
			print ( ret$annotation )
			print ( "Values:" ) 
			ret$data <- x@data[probeset, ]
			print ( ret$data )
			ret$stats <- NULL
			#browser()
			if ( length(x@stats)>0 ) {
				for( i in 1:length(x@stats) ) {
					#		browser()
					ret$stats <- rbind(ret$stats, 
							cbind ( 
									rep(names(x@stats)[i],length(probeset) ), 
									x@stats[[i]][is.na(match( x@stats[[i]][,1], probeset))==F,] 
							) 
					) 
				}
				colnames(ret$stats) <- c('comparison', colnames(x@stats[[1]]) )
				print ( ret$stats )
			}
			invisible(ret)
		})


