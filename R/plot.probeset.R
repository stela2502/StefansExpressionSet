#' @name plot.probeset
#' @aliases plot.probeset,StefansExpressionSet-method
#' @rdname plot.probeset-methods
#' @docType methods
#' @description This function plots the expression data grouped by the GroupName using ggplot2.
#' @description Works only for one probeset.
#' @param x the StefansExpressionSet object
#' @param groupcol the grouping column in the samples data
#' @param colCol the coloring column in the sample data
#' @param probeNames the column in the annotation datacontaining the gene symbol
#' @param x the StefansExpressionSet object
#' @param probeset the probeset name (rownames(x@data))
#' @param boxplot (F or T) create a a dots- or box-plot
#' @param pdf save the file as pdf (default svg)
#' @param geneNameCol the column name for the gene symbol to use in the plots
#' @param sampleGroup the sample table column containing the grouping information
#' @title description of function plot.probeset
#' @examples 
#' p <- plot.probeset( red, rownames(red@data)[20], geneNameCol= 'GeneID', boxplot=T)
#' #produces the file '../../tmp/reducedSet_boxplot_Mybl1_expression.svg'
#' p <- plot.probeset( red, rownames(red@data)[20], geneNameCol= 'GeneID', boxplot=F)
#' #produces the file '../../tmp/reducedSet_points_Mybl1_expression.svg'
#' @export
setGeneric('plot.probeset', ## Name
	function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol", sampleGroup='GroupName' ) { ## Argumente der generischen Funktion
		standardGeneric('plot.probeset') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

####last

setMethod('plot.probeset', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol", sampleGroup='GroupName' ) {
	if ( sum(is.na(match(probeset, rownames(x@data)))==F) == 0 ){
		probeset <- rownames(x@data)[match( probeset, x@annotation[,geneNameCol] ) ]
	}
	if ( length(probeset) == 0 ) {
		stop( "I could not find the probeset in the data object!")
	}
	print ( probeset )
	x <- reduce.Obj(x, probeset )
	if ( nrow(x@data) == 1 ) {
		melted <- melt( x, probeNames=geneNameCol, groupcol = sampleGroup , colCol= sampleGroup )
		colnames(melted) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
		collaps <- '_points_'
		if (boxplot){
			collaps <-'_boxplot_'
		}
		if ( pdf ) {
			pdf( file=paste(x@outpath,x@name,collaps,x@annotation[,geneNameCol],"_expression.pdf",sep='') ,width=5,height=4)
		}else {
			devSVG( file=paste(x@outpath,x@name,collaps,x@annotation[,geneNameCol],"_expression.svg",sep='') ,width=5,height=4)
		}
		if ( boxplot ){
			p <- ggplot(melted ,aes(x=Group, y=Expression,color=Group)) + geom_boxplot(shape=19)+theme_bw()
		}else{
			#TODO: adjust jitter width based on sample group numbers
			p <-ggplot(melted ,aes(x=Group, y=Expression,color=Group)) +geom_jitter(shape=19, width=0.2, height=0)+theme_bw()
		}
		print ( p )
		dev.off()
		invisible(p)
	}
	else {
		stop("You should rather use a heatmap to display more than one probeset")
	}
})



