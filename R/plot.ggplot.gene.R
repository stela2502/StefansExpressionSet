#' @name ggplot.gene
#' @aliases ggplot.gene,StefansExpressionSet-method
#' @rdname ggplot.gene-methods
#' @docType methods
#' @description  Plot one gene in the StefansExpressionSet as boxplot or points plot (using ggplot2)
#' @param dat the StefansExpressionSet object
#' @param gene the gene of interest
#' @param colrs the grouping colors for the x axis (samples)
#' @param groupCol the samples clumn that contains the grouping information
#' @param colCol the sample column that contains the color information
#' @param log2 log2 transform data before plotting (default = F)
#' @title description of function ggplot.gene
#' @export 
setGeneric('ggplot.gene', ## Name
	function (dat,gene, colrs=NULL, groupCol='GroupID', colCol=NULL, boxplot=F, log2=F) { ## Argumente der generischen Funktion
		standardGeneric('ggplot.gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('ggplot.gene', signature = c ( 'StefansExpressionSet') ,
	definition = function (dat,gene, colrs=NULL, groupCol='GroupID', colCol=NULL, boxplot=F, log2=F) {
	not.in = 'NUKL'
	if ( is.null(colrs)){
		if ( is.null(colCol) ) {
			colrs = rainbow( length(unique(dat@samples[,groupCol])))
		}
		else {
			colrs = rainbow( length(unique(dat@samples[,colCol])))
		}	
	}
	dat <- reduce.Obj ( dat, gene, name=gene )
	if ( log2 ) {
		dat@data <- log2(dat@data+1)
	}
	g1 <- melt.StefansExpressionSet( dat, probeNames=dat@rownamescol, groupcol=groupCol,colCol=colCol)
	#g1 <- subset(dat, Gene.Symbol == gene)
	colnames(g1) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' ) [1:ncol(g1)]
	if ( log2 ) {
		exprsname <- 'Expression (log2)'
	}
	else{
		exprsname <- 'Expression'
	}
	g1$Group <- factor( g1$Group, levels=unique( as.character(g1$Group)))
	if ( length(g1) == 0){
		not.in = c( gene )
	}
	if ( boxplot ){
		list ( plot=ggplot2::ggplot(g1 ,ggplot2::aes(x=Group, y=Expression,color=Group)) 
						+ ggplot2::geom_boxplot(shape=19)+ggplot2::theme_bw() 
						+ ggplot2::scale_colour_manual( values= colrs, guide=FALSE )
						+ ggplot2::theme( axis.text.x= ggplot2::element_text(angle=90) )
						+ ggplot2::labs(title=gene, x='', y=exprsname),
				not.in = not.in
		)
		
	}else{
		list ( plot=ggplot2::ggplot(g1,ggplot2::aes(x=Group, y=Expression,color=Group)) 
						+ ggplot2::geom_jitter(shape=19)+ggplot2::theme_bw()
						+ ggplot2::scale_colour_manual( values= colrs, guide=FALSE )
						+ ggplot2::theme( axis.text.x= ggplot2::element_text(angle=90) )
						+ ggplot2::labs(title=gene, x='', y=exprsname),
				not.in = not.in
		)
	}
})

