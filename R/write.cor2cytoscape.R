#' @name cor2cytoscape
#' @aliases cor2cytoscape,StefansExpressionSet-method
#' @rdname cor2cytoscape-methods
#' @docType methods
#' @description  exports the correlation matrix in the format \url{http://www.cytoscape.org/} does
#' @description  support as network file.
#' @param M the correlation matrix obtained by a run of \code{\link{corMat}}
#' @title description of function cor2cytoscape
#' @export 
setGeneric('cor2cytoscape', ## Name
	function (M, file, cut=0.9 ) { ## Argumente der generischen Funktion
		standardGeneric('cor2cytoscape') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

###last

setMethod('cor2cytoscape', signature = c ( 'StefansExpressionSet') ,
	definition = function (M, file, cut=0.9 ) {
	edges <- NULL
	if ( class(M) == 'list' ) {
		edges <- list()
		for( i in 1:length(M) ) {
			edges[[i]] <- cor2cytoscape( M[[i]], file= paste(file,names(M)[i],'.txt',sep='_'), cut=cut) 		
		}
		names(edges) <- names(M)
	}
	else {
		edges <- NULL
		n <- as.vector(rownames(M))
		diag(M) <- 0
		for (i in 1:nrow(M)) {
			for (j in which(M[i,i:ncol(M)] >cut | M[i,i:ncol(M)] < -cut ) ) {
				edges <- rbind(edges, c(n[i], n[j], M[i,j] ))	
			}
		}
		#edges <- matrix( edges, ncol=3 )
		edges <- as.data.frame(edges)
		colnames(edges) <- c('Start', 'End', 'rho')
		edges$Type <- as.vector(rep( 1, nrow(edges) ))
		edges$Type[which(as.numeric(as.vector(edges$rho)) < 0)] <- -1
		write.table( edges, file , sep=" ",quote=F, row.names=F )
	}
	invisible(edges)
})
#
