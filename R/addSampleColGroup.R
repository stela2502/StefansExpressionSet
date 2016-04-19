#' @name addSampleColGroup
#' @aliases addSampleColGroup,StefansExpressionSet-method
#' @rdname addSampleColGroup-methods
#' @docType methods
#' @description an internal method adding new color columns to the melted dataset
#' @param x  the StefansExpressionSet object
#' @param melted the already melted matrix
#' @param colName the colnames to add to the matrix (defult=NULL)
#' @title description of function addSampleColGroup
setGeneric('addSampleColGroup', ## Name
		function ( x, melted, colName=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('addSampleColGroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('addSampleColGroup', signature = c ('StefansExpressionSet'),
		definition = function ( x, melted, colName=NULL ) {
			if ( ! is.null(colName) ) {
				genes = nrow(x@data)
				samples = ncol(x@data)
				grps <- NULL
				for ( GNid in 1:length(colName)){
					le <- genes + GNid -1
					melted_new <- matrix(nrow=(nrow(melted)+GNid*samples), ncol=ncol(melted) )
					for ( sid in 0:(samples-1)) {
						for ( i in (1+(sid*genes)):((sid+1)*genes) ) {
							melted_new[ (i+GNid * sid),] <- as.vector(t(melted[ i,]))
						}
						l <- melted[((sid+1)*genes),]
						l[c(1,3)] <- c(colName[GNid], as.character(x@samples[sid+1,colName[GNid]]) )
						melted_new[ (genes+GNid)*(sid+1), ] <- as.vector(t(l))
					}
					melted <- data.frame(melted_new,row.names= 1:nrow(melted_new))
				}
			}
			melted
		} )
