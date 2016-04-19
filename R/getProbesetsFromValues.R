#' @name getProbesetsFromValues
#' @aliases getProbesetsFromValues,StefansExpressionSet-method
#' @rdname getProbesetsFromValues-methods
#' @docType methods
#' @description Using this function you can select probes based on there expression value in a set of samples
#' @param x the ExpreesionSet object
#' @param v the cut off- or matching value
#' @param sample sample names as used in data column names
#' @param mode one of less, more, onlyless or equals
#' @title description of function getProbesetsFromValues
#' @return a list of probesets that match the requirements
#' @export 
setGeneric('getProbesetsFromValues', ## Name
		function ( x, v='NULL', sample='NULL', mode='less' ) { ## Argumente der generischen Funktion
			standardGeneric('getProbesetsFromValues') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('getProbesetsFromValues', signature = c ('StefansExpressionSet'),
		definition = function ( x, v='NULL', sample='NULL', mode='less' ) {
			s <- FALSE
			if ( is.null(v) ){
				s<-TRUE
			}
			if ( is.null(sample) ){
				s<-TRUE
			}
			if ( s ) { stop ( "Please give me the required values for v and sample") }
			probesets <- NULL
			for ( s in sample ) {
				switch( mode,
						'less' = probesets <- c(probesets, as.vector(rownames(x@data)[which(x@data[,s] <= v)] ) ) ,
						'more' = probesets <- c(probesets, as.vector(rownames(x@data)[which(x@data[,s] > v)] ) ), 
						'onlyless' = probesets <- c(probesets,  as.vector(rownames(x@data)[which(x@data[,s] < v)] ) ),
						'equals' = probesets <- c(probesets, as.vector(rownames(x@data)[which(x@data[,s] == v)] ) )
				)
			}
			unique(probesets)
		})


