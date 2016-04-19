#' @name writeStatTables
#' @aliases writeStatTables,StefansExpressionSet-method
#' @rdname writeStatTables-methods
#' @docType methods
#' @description  export the statistic files
#' @param x the StefansExpressionSet
#' @title description of function writeStatTables
#' @export 
setGeneric('writeStatTables', ## Name
	function ( x, annotation=NULL, wData=F ) { ## Argumente der generischen Funktion
		standardGeneric('writeStatTables') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('writeStatTables', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, annotation=NULL, wData=F ) {
	opath = paste(x@outpath,'stats_Pval/', sep='')
	if ( wData ) {
		opath = paste( x@outpath,'stats_wData/',sep='')
	}
	write.table( x@samples, file= paste(opath,'stats_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	Comparisons <- make.names(names(x@stats))
	system ( paste('mkdir ',opath, sep='') )
	for ( i in 1:length(Comparisons )){
		fname <- paste(opath ,Comparisons[i], '.xls',sep='')
		if ( ! is.null(annotation) ) {
			if ( wData==F ) {
				write.table( cbind( x@annotaion[, annotation],  x@stats[[i]] ) , file= fname, row.names=F, sep="\t",quote=F )
			}else {
				write.table( cbind(x@annotation[,annotation], x@stats[[i]], x@data ), file= fname, row.names=F, sep="\t",quote=F )
			}
		}
		else {
			if ( wData==F ) {
				write.table( x@stats[[i]], file= fname, row.names=F, sep="\t",quote=F )
			}else {
				write.table( cbind(x@stats[[i]], x@data ), file= fname, row.names=F, sep="\t",quote=F )
			}
		}
		print ( paste ( "table ",fname," for cmp ",Comparisons[i],"written" ) )
	}
})

