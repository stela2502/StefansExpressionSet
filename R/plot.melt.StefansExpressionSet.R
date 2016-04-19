#' @name melt
#' @aliases melt,StefansExpressionSet-method
#' @rdname melt-methods
#' @docType methods
#' @description  met an StefansExpressionSet to be plotted using ggplot2 functions
#' @param data the StefansExpressionSet object
#' @param groupcol the column in the samples table to group the expression on
#' @param colCol samples column names(s) to color the grouping data
#' @param rowCol annotation column name(s) describing the gene groups
#' @param probeNames which probenames to use (ProbeSetID or Gene.Symbol ...)
#' @title description of function melt
#' @export 
setGeneric('melt.StefansExpressionSet', ## Name
		package = 'StefansExpressionSet',
	function ( dat, groupcol='GroupName', colCol=NULL, rowCol=NULL, probeNames=NULL,  na.rm = FALSE, value.name = "value") { ## Argumente der generischen Funktion
		standardGeneric('melt.StefansExpressionSet') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)


setMethod('melt.StefansExpressionSet',
	 signature= ('StefansExpressionSet' ),
	 function ( dat, groupcol='GroupName', colCol=NULL,  rowCol=NULL, probeNames=NULL, na.rm = FALSE, value.name = "value" ) {
	if ( is.null(probeNames)){
		probeNames <- dat@rownamescol
	}
	ma  <- dat@data[,order(dat@samples[,groupcol] )]
	#rownames(ma) <- forceAbsoluteUniqueSample(as.vector(dat@annotation[, probeNames]) )
	melted <- reshape2::melt( cbind(rownames(ma),ma) )
	dat@samples <- dat@samples[order(dat@samples[,groupcol]),]
	if ( length( which ( melted[,2] == '') ) > 0 ){
		melted <- melted[ - which ( melted[,2] == ''),]
	}
	melted[,3] <- as.numeric(as.character(melted[,3]))
	grps <- NULL
	for ( i in as.vector(dat@samples[,groupcol]) ){
		grps <- c( grps, rep( i, nrow(dat@data)))
	}
	cgrps <- NULL
	for ( i in as.vector(dat@samples[,colCol]) ){
		cgrps <- c( cgrps, rep( i, nrow(dat@data)))
	}
	colnames(melted) <- c('ProbeName', 'SampleName', 'Expression')
	melted$Group <- grps
	melted <- addSampleColGroup( dat, melted, colName= colCol)
	melted <- addGeneColGroup( dat, melted, colName= rowCol)
	melted
})

