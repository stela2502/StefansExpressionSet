#' @name anayse_working_set
#' @aliases anayse_working_set,NGSexpressionSet-method
#' @rdname anayse_working_set-methods
#' @docType methods
#' @description  zscore the ngs dataset
#' @param m the NGSexpressionSet
#' @return The z scored NGSexpressionSet
#' @title description of function anayse_working_set
#' @export 
setGeneric('anayse_working_set', ## Name
	function ( dataOBJ, name,  p_values = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ),geneNameCol= "mgi_symbol", batchCorrect=NULL) { ## Argumente der generischen Funktion
		standardGeneric('anayse_working_set') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('anayse_working_set', signature = c ('NGSexpressionSet'),
	definition = function ( dataOBJ, name,  p_values = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ),geneNameCol= "mgi_symbol", batchCorrect=NULL) {
	dataOBJ <- preprocess.NGSexpressionSet ( dataOBJ )
	if ( ! is.null(batchCorrect) ){
		removeBatch.NGSexpressionSet( dataOBJ ,batchCorrect )
	}
	png(file=paste(dataOBJ@name,"_",version,'_interesting_samples_clean_up.png',sep=''), width=800, height=800) 
	plot(hclust(as.dist(1-cor(dataOBJ@data))),main=paste(dataOBJ@name,version, ' samples')) 
	dev.off()
	dataOBJ <- do_comparisons.NGSexpressionSet ( dataOBJ, geneNameCol=geneNameCol)
	#dataOBJ@expr <- exprs(dataOBJ@vsdFull)
	dataOBJ <- get_gene_list.NGSexpressionSet(dataOBJ,p_values, geneNameCol=geneNameCol)	
	dataOBJ
})


