
#' @name group.hclust
#' @aliases group.hclust,StefansExpressionSet-method
#' @rdname group.hclust-methods
#' @docType methods
#' @description Create a groping based on a hclust run like hclust(dist(x@data),method="ward.D2")
#' @param x the StefansExpressionSet object
#' @param groups how many groups should I identify
#' @param name the new column name created as  paste( name, groups, 'groups' )
#' @param type a 'gene' or 'sample' grouping?
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
					cols= 1:2
				}
				d = data.frame( 
						x@usedObj[['hclust_gene']][[ name ]]$order, 
						factor(cutree(x@usedObj[['hclust_gene']][[ name ]],k=groups), labels=1:groups )
				)
				colnames(d) <- c(paste( name, 'order' ), paste( name, groups, 'groups' ))
				if(is.na( match( "d[, cols]", colnames(x@annotation)) )==F ){
					colnames(x@annotation)[match( "d[, cols]", colnames(x@annotation))] <- colnames(d)[2]
				}
				x@annotation <- cbind(x@annotation, d[,cols] )
				x <- colors_4 ( x, paste( name, groups, 'groups' ) )
			}
			else if ( type == 'sample' ){
				if ( is.null( x@usedObj[['hclust_sample']])){
					x@usedObj[['hclust_sample']] = list()
				}
				#hclust(dist(isect@data),method="ward.D2")$order
				cols = 2 
				if ( is.null(x@usedObj[['hclust_sample']][[ name ]]) ){
					x@usedObj[['hclust_sample']][[ name ]] <- hclust(distfun(t(x@data)),method=hclustMethod)	
					#	x@usedObj[['hclust_sample']][[ name ]] <- hclust(as.dist( 1- cor(t(x@data), method='pearson') ),method="ward.D2")
					cols= 1:2
				}
				d = data.frame( 
						x@usedObj[['hclust_sample']][[ name ]]$order, 
						factor(cutree(x@usedObj[['hclust_sample']][[ name ]],k=groups) , labels=1:groups )
				)
				colnames(d) <- c(paste( name, 'order' ), paste( name, groups, 'groups' ))
				x@samples <- cbind(x@samples, d[,cols] )
				if(is.na( match( "d[, cols]", colnames(x@samples)) )==F ){
					colnames(x@samples)[match( "d[, cols]", colnames(x@samples))] <- colnames(d)[2]
				}
				x <- colors_4 ( x, paste( name, groups, 'groups' ) )
			}
			else { stop ( 'not implemented') }
			x
		}
)
