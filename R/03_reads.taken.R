#' @name reads.taken
#' @aliases reads.taken,NGSexpressionSet-method
#' @rdname reads.taken-methods
#' @docType methods
#' @description  this check is testing how many percent of the total reads are taken from the first
#' @description  x percent of the genes this function calls reads.taken.NGSexpressionSet internally.
#' @param x the NGSexpressionSet
#' @param genes see reads.taken 
#' @param cutoff = 0.77 a good expression dataset should not use more than 77\% of the reads in the top 5\%
#' @param tf if you supply a list of genes here the same values as for all genes are calulated for this list
#' @return a list of samples which have passed the test
#' @title description of function reads.taken
#' @return a complicated list with al measured values
#' @export 
setGeneric('reads.taken', ## Name
	function ( x, percentile= 0.05, tf=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('reads.taken') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reads.taken', signature = c ('NGSexpressionSet'),
	definition = function ( x, percentile= 0.05, tf=NULL ) {
	top.genes <- list()
	reads.taken <- vector( 'numeric', ncol(x@data))
	nTF <- vector('numeric',  ncol(x@data)+1)
	percentile= 1- percentile
	for ( i in 1:ncol(x@data) ){
		qua <- quantile(x@data[,i], percentile)
		reads.taken [i] <- sum (x@data[which(x@data[,i] > qua),i] ) / sum(x@data[,i])
		top.genes[[i]] <- rownames(x@data) [ which(x@data[,i] > qua) ]
		if ( ! is.null(tf) ){
		nTF[i] <- length( intersect( tf, str_replace(top.genes[[i]],  "\\.[0-9]*" ,'') ) )
		}
	}
	names( reads.taken) <- colnames(x@data)
	names(top.genes) <- colnames(x@data)
		
	inter <- intersect( top.genes[[1]], top.genes[[2]])
	for (i in 3:ncol(x@data) ) { inter <- intersect( inter, top.genes[[i]]) }
	if ( ! is.null(tf) ) {nTF[ncol(x@data)+1] <- length(intersect(str_replace(inter,  "\\.[0-9]*" ,''), tf ) )}
	reads.taken.intersect <- vector( 'numeric', ncol(x@data))
	for ( i in 1:ncol(x@data) ){
		reads.taken.intersect [i] <- sum ( x@data[inter ,i] ) / sum(x@data[,i])
	}
	names( reads.taken.intersect) <- colnames(x@data)
	list( reads.taken = reads.taken, top.genes = top.genes, intersect = inter,reads.taken.intersect = reads.taken.intersect, nTF = nTF )
})
	

