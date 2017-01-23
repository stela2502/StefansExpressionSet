#' @name plot.histogram
#' @aliases plot.histogram,StefansExpressionSet-method
#' @rdname plot.histogram-methods
#' @docType methods
#' @description This function plots one gene as histogram to check whether there ar4e clear expression differences in different plates.
#' @param dataObj the StefansExpressionSet object
#' @param probesetID the probeset id of the gene to plot
#' @param cuts the cuts are used for the 1D gene groups default=vector('list',1)
#' @param subpath the subpath to plot to ( default = preprocess)
#' @param colGroup the samples table column to color the plot
#' @param nameCol the gene name column to enhance the plot information
#' @param png create a png file (default =F)
#' @param breaks the amount of breaks in the hist default=15
#' @title description of function plot.histograms
#' @export 
setGeneric('plot.histogram', ## Name
		function ( dataObj, probesetID, cuts=vector('list',1), subpath='preprocess', colGroup='ArrayID', nameCol='gene_name', png=FALSE, breaks=15 ) { ## Argumente der generischen Funktion
			standardGeneric('plot.histograms') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.histogram', signature = c ('StefansExpressionSet'),
		definition = function ( dataObj, probesetID, cuts=vector('list',1), subpath='preprocess', colGroup='ArrayID', nameCol='gene_name', png=FALSE,breaks=15 ) {
			
			ma <- dataObj@data
			#if ( dataObj@wFACS ){
			#	ma <- rbind( ma,  dataObj@facs )
			#}
			n <- rownames(ma)
			names = names(table(dataObj@samples[,colGroup]))
			arrays <- length(names)
			dataObj <- colors_4(dataObj,colGroup)
			cols <- dataObj@usedObj$colorRange[[colGroup]]
			n.cuts <- names(cuts)
			if ( png ){
				opath = file.path(dataObj@outpath,subpath )
				dir.create(opath, showWarnings = FALSE)
			}
			
			plot.this <- function( i ) {
				temp <- vector('list',arrays)
				m <- NULL
				for (a in names) {
					temp[[a]] <- density(t(ma[i,which(dataObj@samples[,colGroup] == a )]))
					m <- c(m,max(temp[[a]]$y))
				}
				#h <- hist(ma[i,],main=n[i], xlab='expression values [raw]', freq=F, col=rgb(0, 1, 0, 0.5), cex.lab = 1.5, breaks = 15, ylim=c(0,max(m)) )
				h <- hist(t(ma[i,]), breaks = breaks,plot=F ) #, main= paste(dataObj@annotation[i,nameCol], i) )
				m <- c(m, max(h$density) )
				hist(t(ma[i,]), breaks = breaks, freq=F,
						main= paste(dataObj@annotation[i,nameCol], i), 
						col=rgb(0, 1, 0, 0.5), xlab="Ct", cex.lab = 1.5, 
						ylim=c(0,max(m))  
				)
				id = 1
				for (a in names ) {
					lines( temp[[a]] , col=cols[id], lwd=2)
					id = id +1
				}
				pos <- which( n.cuts == n[i] )
				if ( length(pos) > 0 ){
					for (c in 1:length(cuts[[pos]]) ) {
						abline( v= cuts[[pos]][c], col='black', lwd = 3, lty = 2 )
					}
				}
			}
			for ( i in probesetID ) {
				if ( png ) {
					png( file=file.path( opath, paste(i,'png',sep='.')),width=800, height=800 )
				}
				try(plot.this ( i ))
				if ( png ) {
					dev.off()
				}
			}
		} 
)

