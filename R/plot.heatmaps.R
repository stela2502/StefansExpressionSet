#' @name plot.heatmaps
#' @aliases plot.heatmaps,StefansExpressionSet-method
#' @rdname plot.heatmaps-methods
#' @docType methods
#' @description An rather ancient, but flexible heatmapping function that depends on heatmap.2 and allows for many data
#' @description selection/conversion options
#' @param dataOBJ the StefansExpressionSet object
#' @param groupCol the column in the samples table to group the samples on
#' @param gene.names the rownames(dataObj@data) level gene names
#' @param pvalue an optional cut off value to select genes from the statistical results tables
#' @param analysis_name the name for the outfiles
#' @param Subset an optional list of strings matching to the gene symbols used to select genes of interest
#' @param collaps collaps the sample groups into a single column of the heatmap using one of collaps=c('median', 'mean')
#' @title description of function plot.heatmaps
#' @export
setGeneric('plot.heatmaps', ## Name
	function ( dataOBJ, groupCol='GroupName', gene.names=NULL , pvalue=1, analysis_name =NULL, Subset=NULL, collaps=NULL, pdf=F,... ) { ## Argumente der generischen Funktion
		standardGeneric('plot.heatmaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.heatmaps', signature = c ( 'StefansExpressionSet') ,
	definition = function ( dataOBJ, groupCol='GroupName', gene.names=NULL , pvalue=1, analysis_name =NULL, Subset=NULL, collaps=NULL, pdf=F,... ) {
	dataOBJ <- normalize(dataOBJ)
	dataOBJ <- z.score( dataOBJ )
	dataOBJ = sd.filter(dataOBJ)
	
	if ( is.null(analysis_name)){
		analysis_name = dataOBJ@name
	}
	
	if ( ! is.null(gene.names)){
		dataOBJ <- reduce.Obj( dataOBJ, gene.names, name= dataOBJ@name )
	}
	
	if ( ! is.null(Subset) ) { ## 
		useful <- NULL
		for ( i in 1:length(Subset) ){
			useful <- c( useful, grep(Subset[i], rownames(dataOBJ@data) ) )
		}
		if ( length(useful) > 0 ){
			dataOBJ <- reduce.Obj(dataOBJ, rownames(dataOBJ@data)[useful], name= dataOBJ@name)
		}
	}
	
	if ( ! is.null(collaps)){
		dataOBJ <- collaps(dataOBJ, by=collaps, groupCol = groupCol )
	}else {
		collaps= ''
	}
	
	## estimate the plot dimensions
	s <- ceiling(20 * nrow(dataOBJ@data)/230 )
	if ( s < 5 ) {s <- 5}
	if ( pdf ) {
		pdf( file=file.path(dataOBJ@outpath,paste(dataOBJ@name,collaps,"_Heatmap_",pvalue,".pdf",sep='')) ,width=10,height=s)
	}else {
		devSVG( file=file.path(dataOBJ@outpath,paste(dataOBJ@name,collaps,"_Heatmap_",pvalue,".svg",sep='')) ,width=10,height=s)
	}
	colnames(dataOBJ@data) <- dataOBJ@samples[,groupCol]
	heatmap.2( as.matrix(dataOBJ@data), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", 
			distfun = function (x) {as.dist( 1- cor(t(x), method='pearson') ) } ,...)
	
	dev.off()
	
	#write.table( data.frame ( 'GeneSymbol' = rownames(geneNamesBased),geneNamesBased[,-1]),file= paste(fname,collaps,'_HeatmapID_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	#write.table( geneNamesBased, file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	#write.table( dataOBJ@samples, file= paste(fname,collaps,'_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	
})


