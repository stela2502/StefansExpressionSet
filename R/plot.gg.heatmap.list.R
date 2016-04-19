#' @name gg.heatmap.list
#' @aliases gg.heatmap.list,StefansExpressionSet-method
#' @rdname gg.heatmap.list-methods
#' @docType methods
#' @description  uses ggplot2 to plot heatmaps
#' @param dat the StefansExpressionSet object
#' @param glist a list of probesets to plot (or all)
#' @param colrs a list of colors for the sample level boxes (or rainbow colors)
#' @param groupCol the column group in the samples table that contains the grouping strings
#' @param colCol the column group in the samples table that contains the color groups
#' @title description of function gg.heatmap.list
#' @export 
setGeneric('gg.heatmap.list', ## Name
		function (dat,glist=NULL, colrs=NULL, groupCol='GroupID', colCol=NULL) { ## Argumente der generischen Funktion
			standardGeneric('gg.heatmap.list') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('gg.heatmap.list', signature = c ( 'StefansExpressionSet') ,
		definition = function (dat,glist=NULL, colrs=NULL, groupCol='GroupID', colCol=NULL) {
			
			if ( ! is.null(glist) ) {
				isect <- reduce.Obj ( dat, glist)
			}else {
				isect <- dat
			}
			if ( is.null(colrs) ){
				colrs = rainbow( length(unique(isect@samples[,colCol])))
			}
			if ( ! isect@zscored ) {isect <- z.score(isect)}
			dat.ss <- melt.StefansExpressionSet ( isect, probeNames=isect@rownamescol, groupcol=groupCol,colCol=colCol)
			#dat.ss <- dat[which(is.na(match(dat$Gene.Symbol,isect))==F),]
			colnames(dat.ss) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 
					paste('ColorGroup', 1:10) )[1:ncol(dat.ss)]
			r <- defineHeatmapColors(dat, dat.ss )
			dat.ss <- r$melted
			ord.genes <- rownames(isect@data)[hclust(dist(isect@data),method="ward.D2")$order]
			if ( ! is.null(colCol) ){
				ord.genes <- c( ord.genes,colCol  )
			}
			dat.ss$Gene.Symbol <- with(dat.ss,factor(Gene.Symbol,levels =
									unique(as.character(ord.genes))))
			dat.ss$Sample <- with(dat.ss,factor(Sample,levels =
									unique(as.character(Sample))))
			dat.ss$Group <- with(dat.ss,factor(Group,levels =
									unique(as.character(Group))))
			dat.ss$colrss <- colrs[dat.ss$Group]
			ss <-dat.ss[which(dat.ss$Gene.Symbol==dat.ss$Gene.Symbol[1]),]
			p = ( ggplot2::ggplot(dat.ss, ggplot2::aes(x=Sample,y=Gene.Symbol))
						+ ggplot2::geom_tile(ggplot2::aes(fill=Expression)) 
						+ ggplot2::scale_fill_manual( values =  r$colors ) 
						+ ggplot2::theme(
								legend.position= 'bottom',
								axis.text.x=ggplot2::element_blank(),
								axis.ticks.length=ggplot2::unit(0.00,"cm")
						)+ ggplot2::labs( y='') )
			if ( ncol(dat.ss) == 6 ){
				p <- p + ggplot2::facet_grid( colrss ~ Group,scales="free", space='free')
			}else if ( ncol(dat.ss) == 5 ) {
				p <- p + ggplot2::facet_grid( . ~ Group,scales="free", space='free')
			}
			list ( plot = p, not.in = setdiff( glist, rownames(isect@data)) )
		}
)


