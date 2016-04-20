#' @name defineHeatmapColors
#' @aliases defineHeatmapColors,SingleCellsNGS-method
#' @rdname defineHeatmapColors-methods
#' @docType methods
#' @description  uses ggplot2 to plot heatmaps
#' @param merged the merged data object with the Expression column that should be colored
#' @param colrs and optional colors vector( gray + bluered for data and rainbow for samples)
#' @title description of function gg.heatmap.list
#' @return a list with the modified merged table and the colors vector
#' @export 
setGeneric('defineHeatmapColors', ## Name
		function (x, melted,..., colrs=NULL) { ## Argumente der generischen Funktion
			standardGeneric('defineHeatmapColors') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('defineHeatmapColors', signature = c ( 'SingleCellsNGS') ,
		definition = function (x, melted, colrs=NULL ){
			if ( is.factor( melted$Expression )) {
				## here might be some row grouping going on!
				d <- levels(melted$Expression)[melted$Expression]
				prob.id <- which(is.na(as.numeric(d))==T)
				treat.separate <- unique(d[prob.id])
				n <- as.numeric(d[-prob.id])
				m <- minValueExpr( x)
				brks= c( (m -.1), as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				d[-prob.id]  <- brks [cut( n, breaks= brks)]
				melted$Expression <- factor( d, levels= c(brks, treat.separate ) )
				colors= c(
						'gray', 
						gplots::bluered(length(brks) -2  ), ## the expression
						rainbow( length(treat.separate) ) ## the sample descriptions
				)
			}
			else {
				n <- as.numeric(melted$Expression )
				m <- minValueExpr( x)
				brks= c( (m-.1), as.vector(quantile(n[which(n != m)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				melted$Expression <- factor( brks [cut( n, breaks= brks)] , levels= c(brks) )
				
				colors= c(
								'gray', 
								gplots::bluered(length(brks) -2  ) ## the expression
				)
			}
			list (melted = melted, colors = colors)
		}
)


setMethod('defineHeatmapColors', signature = c ( 'StefansExpressionSet') ,
		definition = function (x, melted, colrs=NULL ){
			if ( is.factor( melted$Expression )) {
				## here might be some row grouping going on!
				d <- levels(melted$Expression)[melted$Expression]
				prob.id <- which(is.na(as.numeric(d))==T)
				treat.separate <- unique(d[prob.id])
				n <- as.numeric(d[-prob.id])
				brks= c( -20.1, as.vector(quantile(n[which(n != -20)],seq(0,1,by=0.1)) ))
				brks = unique(brks)
				d[-prob.id]  <- brks [cut( n, breaks= brks)]
				melted$Expression <- factor( d, levels= c(brks, treat.separate ) )
				colors= c(
						gplots::bluered(length(brks) -1  ), ## the expression
						rainbow( length(treat.separate) ) ## the sample descriptions
				)
			}
			else {
				n <- as.numeric(melted$Expression )
				brks= c(  as.vector(quantile(n,seq(0,1,by=0.1)) ))
				brks = unique(brks)
				melted$Expression <- factor( brks [cut( n, breaks= brks)] , levels= c(brks) )
				
				colors= c(
						'gray', 
						gplots::bluered(length(brks) -2  ) ## the expression
				)
			}
			list (melted = melted, colors = colors)
		}
)

#' @name minValueExpr
#' @aliases minValueExpr,SingleCellsNGS-method
#' @rdname minValueExpr-methods
#' @docType methods
#' @description just get me the min value for the object
#' @param x data object
#' @title description of function gg.heatmap.list
#' @return a list with the modified merged table and the colors vector
#' @export
setGeneric('minValueExpr', ## Name
		function (x) { ## Argumente der generischen Funktion
			standardGeneric('minValueExpr') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('minValueExpr', signature = c ( 'SingleCellsNGS') ,
		definition = function (x ){
			m <- 0
			if ( x@zscored ) {
				m <- -20
			}
			m
		}
)
setMethod('minValueExpr', signature = c ( 'StefansExpressionSet') ,
		definition = function (x ){
			m <- 0
			if ( x@zscored ) {
				m <- -324772345234 # do not use!
			}
			m
		}
)
