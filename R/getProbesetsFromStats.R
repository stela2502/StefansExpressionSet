#' @name getProbesetsFromStats
#' @aliases getProbesetsFromStats,StefansExpressionSet-method
#' @rdname getProbesetsFromStats-methods
#' @docType methods
#' @description  getProbesetsFromStats returns a list of probesets (the rownames from the data matrix)
#' @description  for a restriction of a list of stat comparisons probes <- getProbesetsFromStats (
#' @description  x, v=1e-4, pos="adj.P.Val" )
#' @param v The cutoff value
#' @param pos The column in the stats tables to base the selection on
#' @param Comparisons A list of comparisons to check (all if left out)
#' @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
#' @return a list of probesets that shows an adjusted p value below 1e-4
#' @title description of function getProbesetsFromStats
#' @export 
setGeneric('getProbesetsFromStats', ## Name
		function ( x, v=1e-4, pos='padj', mode='less', Comparisons=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('getProbesetsFromStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('getProbesetsFromStats', signature = c ('StefansExpressionSet'),
		definition = function ( x, v=1e-4, pos='padj', mode='less', Comparisons=NULL ) {
			if ( is.null(Comparisons)){	Comparisons <- names(x@stats) }
			probesets <- NULL
			for ( i in match(Comparisons, names(x@stats) ) ) {
				switch( mode,
						'less' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] <= v),1] )),
						'more' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] > v),1] )), 
						'onlyless' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] < v),1] )),
						'equals' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] == v),1] ))
				)
			}
			unique(probesets)
		})

