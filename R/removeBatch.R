#' @name removeBatch
#' @aliases removeBatch,NGSexpressionSet-method
#' @rdname removeBatch-methods
#' @docType methods
#' @description  merge the groups in the dataset into one single column The sample description is
#' @description  restricted to the first entry in each group. Most variables in the sample description
#' @description  table might be useless after this collape
#' @param dataObj the NGSexpressionSet
#' @param what merge by 'median','mean','sd' or 'sum'
#' @return the collapsed NGSexpressionSet
#' @title description of function removeBatch
#' @export 
setGeneric('removeBatch', ## Name
	function ( x, phenotype ) { ## Argumente der generischen Funktion
		standardGeneric('removeBatch') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('removeBatch', signature = c ('NGSexpressionSet'),
	definition = function ( x, phenotype ) {
	if ( x@batchRemoved==1) {
		return (x)
	}
	browser()
	exprs =  dataOBJ$cds
	null <- which ( exprs == 0)
	exprs[null]<- 1
	log <- log(exprs)
	svseq <-  ComBat(dat=filtered.log , mod=mod1, batch=x@samples[,phenotype], par.prior=TRUE, prior.plots=FALSE )
	svseq <- exp(log)
	svseq[null] = 0
	x@cds <- svseq
	x@batchRemoved = 1
	x
})

