#' @name coexprees4gene
#' @aliases coexprees4gene,StefansExpressionSet-method
#' @rdname coexprees4gene-methods
#' @docType methods
#' @description  Calculate the coexpression for any given gene
#' @param x the StefansExpressionSet varibale
#' @param method any method supported by \link[=cor.test]{cor.test}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param padjMethod the method to calucate the FDR with \link[=p.adjust]{p.adjust}
#' @return a data.frame with the columns 'GeneSymbol', 'pval', 'cor', 'adj.p'
#' @title description of function coexprees4gene
#' @export 
setGeneric('coexprees4gene', ## Name
	function ( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) { ## Argumente der generischen Funktion
		standardGeneric('coexprees4gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('coexprees4gene', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) {
	ret <- NULL
	if ( ! is.null(gene) ){
		z <- as.vector( t(x@data[ gene[1] ,]) )
		pval <- vector( 'numeric', nrow(x@data))
		cor <- vector( 'numeric', nrow(x@data))
		for ( i in 1:nrow(x@data) ) {
			try( {	res <-  cor.test( z, as.vector(t(x@data[i,]), 'numeric') ,method=method)
			pval[i] <- res$p.value
			cor[i] <- res$estimate }, silent=T
			)
		}
		adj.p <- p.adjust(pval , method = padjMethod) #BH
		ret <- data.frame( pval = pval, cor = cor, adj.p = adj.p )
		ret <- cbind(x@annotation[,geneNameCol], ret )
		rownames(ret) <- rownames(x@data)
	}
	ret
})
 
