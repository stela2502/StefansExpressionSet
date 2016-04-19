#' @name addGeneColGroup
#' @aliases addGeneColGroup,StefansExpressionSet-method
#' @rdname addGeneColGroup-methods
#' @docType methods
#' @description an internal function to add row level clusters to a merged data table
#' @param x the StefansExpressionSet
#' @param melted the melted table
#' @param colName the annotation column names to use for the coloring default=NULL
#' @title description of function addGeneColGroup
setGeneric('addGeneColGroup', ## Name
		function ( x, melted, colName=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('addGeneColGroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('addGeneColGroup', signature = c ('StefansExpressionSet'),
		definition = function ( x, melted, colName=NULL ) {
			if ( ! is.null(colName)) {
				grps <- NULL
				datarows = nrow(x@data)
				ret <- NULL
				for ( GNid in 1:length(n)){
					add  <- as.matrix(melted[1:datarows,])
					add[,3] <- as.vector(x@annotation[,colName[GNid]])
					add[,4] <-'GeneGroup'
					ret <- rbind( ret, add )
				}
			}else {
				ret <- melted
			}
			ret
		} )




