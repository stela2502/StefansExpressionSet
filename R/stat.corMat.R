#' @name corMat
#' @aliases corMat,StefansExpressionSet-method
#' @rdname corMat-methods
#' @docType methods
#' @description  calculate a correlation matrix for the data object
#' @param x the StefansExpressionSet varibale
#' @param sd_cut the cut off value for the sd check
#' @param method any method supported by \code{\link[stats]{cor.test}}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param groupCol the column name of the grouping variable in the samples table
#' @param name the name of the analysis
#' @return the correlation matrix
#' @title description of function corMat
#' @export 
setGeneric('corMat', ## Name
	function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) { ## Argumente der generischen Funktion
		standardGeneric('corMat') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('corMat', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	d <- reduce.Obj( x, rownames(x@data)[which( apply(x@data,1,sd) > sd_cut)], name = name )
	if ( ! is.null(groupCol) ){
		ret <- list()
		names <- unique(d$samples[,groupCol])
		for ( i in 1:length(names)) {
			a <- subset( d, column=groupCol, value=names[i], name= paste(d$name,names[i],sep='_'), mode='equals' )
			ret[[i]] = corMat( a, sd_cut= sd_cut,method=method, geneNameCol=geneNameCol )
		}
		names(ret) <- namestime.col
		ret
	}
	else {
		n = nrow(d$data)
		print ( paste(d$name,": I create a",n,'x', n,'matrix') )
		ret <- cor(t(d$data), method=method )
		colnames(ret) <- rownames(ret) <- forceAbsoluteUniqueSample( as.vector(d$annotation[,geneNameCol]) )

		ret
	}
})
