#' @name clusters
#' @aliases clusters,StefansExpressionSet-method
#' @rdname clusters-methods
#' @docType methods
#' @description Culters the data either based on the raw data values or any MDS data type and adds the grouping into the samples table.
#' @param dataObj the StefansExpressionSet object
#' @param clusterby cluster on raw data or MDS clustered data default="raw"
#' @param useGrouping do nothing and simply use tis grouping default=NULL
#' @param groups.n how many groups should be detected default= 3
#' @param ctype cluster type - either 'hierarchical clust' or 'kmeans' default = 'hierarchical clust'
#' @param onwhat this option has been kept for the Fluidigm data as there FACS data can also be used default = 'Expression'
#' @param cmethod the method to used with the hclust clustering (default = 'ward.D2')
#' @param name the name for the new grouping (default = 'auto_clusters.1:n')
#' @title description of function clusters
#' @export 
setGeneric('clusters', ## Name
	function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3, 
			ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {## Argumente der generischen Funktion
		standardGeneric('clusters') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('clusters', signature = c ('StefansExpressionSet'),
	definition = function (dataObj,clusterby="raw", useGrouping=NULL, groups.n = 3,
				ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D2', name=NULL ) {
	
			clusters <- NULL
			hc <- NULL
			if(onwhat=="Expression"){
				tab <- dataObj@data
			}
			else {
				stop( paste("Sorry, the mds.type",mds.type,"is not supported") )
			}
			if ( ! is.null(useGrouping) ) {
				clusters <- dataObj@samples[,useGrouping]
				if ( is.factor( clusters)) {
					clusters = as.numeric(clusters)
				}
				dataObj <- colors_4 (dataObj, useGrouping )
			}
			else if ( clusterby=="raw"){#...do mds on tab
				if ( ctype=='hierarchical clust'){
					hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
					clusters <- cutree(hc,k=groups.n)
				}else if (  ctype=='kmeans' ) {
					hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
					clusters <- kmeans( dataObj@usedObj[['mds.proj']] ,centers=groups.n)$cluster
				}else if ( ctype =='mclust' ) {
					hc <- hc( as.dist( 1- cor(t(tab), method='pearson') ) )
					clusters <- hclass(hc, 12)
				}
				else { stop( paste('ctype',ctype, 'unknown!' ) )}
			}else { ## now the clusterby is a MDS algorithm name / MDS dataset name
				if ( is.null( dataObj@usedObj$MDS[[clusterby]] ) ) {
					dataObj <- mds( dataObj, onwhat="Expression", mds.type=clusterby)
				}
				if ( ctype=='hierarchical clust'){
					hc <- hclust(dist( dataObj@usedObj$MDS[[clusterby]] ),method = cmethod)
					clusters <- cutree(hc,k=groups.n)
				}else if (  ctype=='kmeans' ) {
					hc <- hclust(dist( dataObj@usedObj$MDS[[clusterby]] ),method = cmethod)
					clusters <- kmeans( dataObj@usedObj$MDS[[clusterby]] ,centers=groups.n)$cluster
				}else if ( ctype =='mclust' ) {
					hc <- hc( dataObj@usedObj$MDS[[clusterby]] )
					clusters <- hclass(hc, groups.n)
				}
				else { stop( paste('ctype',ctype, 'unknown!' ) )}
			}
			if ( is.null(useGrouping) ){
				## define the group name n and populate the samples table
				if ( is.null(name)){
					if(is.null(dataObj@usedObj[['auto_clusters']])){
					dataObj@usedObj[['auto_clusters']] = 0
				}
				dataObj@usedObj[['auto_clusters']] <- dataObj@usedObj[['auto_clusters']] +1
				name <- paste( 'auto_clusters', 
						dataObj@usedObj[['auto_clusters']] ,sep='.')
				}
				dataObj@samples <- cbind ( dataObj@samples, clusters )
				colnames(dataObj@samples)[ncol(dataObj@samples)] = name
				clusters <- dataObj@usedObj[['clusters']]
				dataObj@usedObj$usedGrouping <- name
				dataObj <- colors_4(dataObj, name )
				print ("used a new grouing")
			}else {
				print ( "reusing old grouping" )
				dataObj@usedObj$usedGrouping <- useGrouping
			}
			## now I want to create some gene clusters too based on hclust only
#			if ( is.null(dataObj@annotation$'hclust Order')){
#				hcG <- hclust(as.dist( 1- cor(dataObj@data, method='pearson') ),method = cmethod )
#				dataObj@annotation$'hclust Order' <- hcG$order
#				dataObj@annotation$'hclust 5 groups' <- factor(cutree(hcG,k=5) )
#				dataObj@annotation$'hclust 10 groups' <- factor(cutree(hcG,k=10) )
#				for ( i in c('hclust Order', 'hclust 5 groups', 'hclust 10 groups' )){
#					dataObj <- colors_4(dataObj, i )
#				}
#			}
			dataObj@usedObj[['clusters']] <- clusters
			dataObj@usedObj[['hc']] <- hc
			dataObj
} )
