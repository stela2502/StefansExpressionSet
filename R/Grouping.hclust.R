
#' @name group.hclust
#' @aliases group.hclust,StefansExpressionSet-method
#' @rdname group.hclust-methods
#' @docType methods
#' @description Create a groping based on a hclust run like hclust(dist(x@data),method="ward.D2")
#' @param x the StefansExpressionSet object
#' @param groups how many groups should I identify
#' @param name the new column name created as  paste( name, groups, 'groups' )
#' @param type a gene or sample grouping?
#' @title description of function group.hclust
#' @export 
setGeneric('group.hclust', ## Name
		function ( x, groups, name, type='genes', distfun= function(a) {dist(a)}, hclustMethod='ward.D2' ) { ## Argumente der generischen Funktion
			standardGeneric('group.hclust') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('group.hclust', signature = c ('StefansExpressionSet'),
		definition = function ( x, groups, name, type='gene', distfun= function(a) {dist(a)}, hclustMethod='ward.D2' ) {
			if ( type == 'gene' ){
				if ( is.null( x@usedObj[['hclust_gene']])){
					x@usedObj[['hclust_gene']] = list()
				}
				#hclust(dist(isect@data),method="ward.D2")$order
				cols = 2 
				if ( is.null(x@usedObj[['hclust_gene']][[ name ]]) ){
					x@usedObj[['hclust_gene']][[ name ]] <- hclust(distfun(x@data),method=hclustMethod)	
				#	x@usedObj[['hclust_gene']][[ name ]] <- hclust(as.dist( 1- cor(t(x@data), method='pearson') ),method="ward.D2")
					cols= 1:2
				}
				d = data.frame( 
						x@usedObj[['hclust_gene']][[ name ]]$order, 
						cutree(x@usedObj[['hclust_gene']][[ name ]],k=groups) 
				)
				colnames(d) <- c(paste( name, 'order' ), paste( name, groups, 'groups' ))
				x@annotation <- cbind(x@annotation, d[,cols] )
				x <- colors_4 ( x, paste( name, groups, 'groups' ), function( x ) { rainbow(x) } )
			}
			else { stop ( 'not implemented') }
			x
		}
)
