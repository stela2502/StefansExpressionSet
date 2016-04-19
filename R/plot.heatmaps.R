#' @name plot.heatmaps
#' @aliases plot.heatmaps,StefansExpressionSet-method
#' @rdname plot.heatmaps-methods
#' @docType methods
#' @description  A flexible heatmapping function that depends on heatmap.2 and allows for many data
#' @description  selection/conversion options
#' @param dataOBJ the StefansExpressionSet object
#' @param gene.names the rownames(dataObj@data) level gene names
#' @param pvalue an optional cut off value to select genes from the statistical results tables
#' @param analysis_name the name for the outfiles
#' @param gene_centered in case there are multiple probesets for each gene - sum the data or display each probset
#' @param Subset an optional list of strings matching to the gene symbols used to select genes of interest
#' @param collaps collaps the sample groups into a single column of the heatmap using one of collaps=c('median', 'mean')
#' @title description of function plot.heatmaps
setGeneric('plot.heatmaps', ## Name
	function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F,... ) { ## Argumente der generischen Funktion
		standardGeneric('plot.heatmaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.heatmaps', signature = c ( 'StefansExpressionSet') ,
	definition = function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F,... ) {
	dataOBJ <- normalize(dataOBJ)
	dataOBJ <- z.score( dataOBJ )
	
	if ( is.null(analysis_name)){
		analysis_name = dataOBJ@name
	}
	exprs = dataOBJ@data
	drop <- which( (apply(exprs,1,sd) == 0) == T)
	if ( length(drop) > 0 ){
		print ( paste("I need to drop",length(drop),"genes as the varianze is 0") )
		exprs = exprs[-drop,]
	}
	if ( ! is.null(gene.names)){
		geneNamesBased <- exprs[is.na(match(rownames(exprs),gene.names$all$id ))==F,]
		geneNamesBased <- data.frame ( 'GeneSymbol' = as.vector (dataOBJ@annotation[is.na(match ( dataOBJ@annotation$GeneID, gene.names$all$id ) ) == F, geneNameCol]),  geneNamesBased) 
	}
	else {
		geneNamesBased <- data.frame( 'GeneSymbol' = dataOBJ@annotation[rownames(exprs),geneNameCol], exprs )
	}
	if ( ! is.null(Subset) ) { ## 
		useful <- NULL
		for ( i in 1:length(Subset) ){
			useful <- c( useful, grep(Subset[i], geneNamesBased[,1] ) )
		}
		if ( length(useful) > 0 ){
			geneNamesBased <- geneNamesBased[ unique(useful), ]
		}
	}
	if ( gene_centered ){
		ok <- as.vector(which ( is.na(geneNamesBased[,1] )))
		grepped <- grep('_drop',forceAbsoluteUniqueSample(as.vector(geneNamesBased[,1]), '_drop_' ))
		if ( length(ok) > 1 ){
			geneNamesBased <- geneNamesBased[- grepped [ - match ( ok[-1] , grepped ) ], ]
		}
		else if ( length(grepped) > 0 ) {
			geneNamesBased <- geneNamesBased[-grepped, ]
		}
	}
	fname = paste(dataOBJ@outpath, 'GenesOnly_',dataOBJ@name,sep='')
	rn <- rownames(geneNamesBased)
	
	if ( ! is.null(collaps)){
		u <- unique(as.vector(dataOBJ@samples$GroupName))
		m <- length(u)
		mm <-  matrix ( rep(0,m * nrow(geneNamesBased)), ncol=m)
		colnames(mm) <- u
		rownames(mm) <- rownames(geneNamesBased)
		if ( collaps=="median") {
			for ( i in u ){
				## calculate the mean expression for each gene over all cells of the group
				mm[,i] <- apply( geneNamesBased[ , which(as.vector(dataOBJ@samples$GroupName) == i )+1],1,median)
			}
		}
		else {
			for ( i in u ){
				## calculate the mean expression for each gene over all cells of the group
				mm[,i] <- apply( geneNamesBased[ , which(as.vector(dataOBJ@samples$GroupName) == i )+1],1,mean)
			}
		}
		geneNamesBased <- data.frame( GeneSymbol = as.vector(geneNamesBased[,1]), mm )
	}else {
		collaps= ''
	}
	
	write.table( data.frame ( 'GeneSymbol' = rownames(geneNamesBased),geneNamesBased[,-1]),file= paste(fname,collaps,'_HeatmapID_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table(  geneNamesBased, file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table( dataOBJ@samples, file= paste(fname,collaps,'_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	
	if ( nrow(geneNamesBased) > 1 ){
		geneNamesBased <- data.frame(read.delim( file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''), header=T ))
		rownames(geneNamesBased) <- rn
		s <- ceiling(20 * nrow(geneNamesBased)/230 )
		if ( s < 5 ) {s <- 5}
		print (paste ("Height =",s))
		if ( pdf ) {
			pdf( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_IDs_",pvalue,".pdf",sep='') ,width=10,height=s)
		}else {
			devSVG( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_IDs_",pvalue,".svg",sep='') ,width=10,height=s)
		}
		print ( paste ("I create the figure file ",dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_IDs_",pvalue,".svg",sep=''))
		heatmap.2(as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", margin=c(10, 10),dendrogram="row",scale="row",
'distfun'= function (x) as.dist( 1- cor(t(x), method='pearson')),Colv=F, main=paste( analysis_name, pvalue ),...)
		dev.off()
		
		rownames(geneNamesBased) <- paste(geneNamesBased[,1] , rownames(geneNamesBased) )
		if ( pdf ) {
			pdf( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_GeneSymbols_",pvalue,".pdf",sep='') ,width=10,height=s)
		}else{
			devSVG( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_GeneSymbols_",pvalue,".svg",sep='') ,width=10,height=s)
		}
		
		print ( paste ("I create the figure file ",dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_GeneSymbols_",pvalue,".svg",sep=''))
		heatmap.2( as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", 
distfun = function (x) {as.dist( 1- cor(t(x), method='pearson') ) } ,...)
		dev.off()
	}
	else {
		print ( "Problems - P value cutoff to stringent - no selected genes!" );
	}
})


