# TODO: implement! Add Documentation when ready

#' Author: Stefan Lang
#' this implements clustering of genes for the ExpressionSet class
#' the clustering options are taken from the SCExV R lib
mds.and.clus <-function(dataObj,clusterby="raw",mds.type="PCA", groups.n, LLEK=2, cmethod='ward.D', ctype='hierarchical clust',onwhat="Expression",... ) {
	if(onwhat=="Expression"){
		tab <- dataObj$z$PCR
	} 
	else {
		print ( paste ( "I work on the FACS data!" ) )
		tab <- dataObj$FACS
	}
	mds.proj <- NULL
	pr <- NULL
	system ( 'rm  loadings.png' )
	if(mds.type == "PCA"){
		pr <- prcomp(tab)
		mds.proj <- pr$x[,1:3]
		png ( file='loadings.png', width=1000, height=1000 )
		plot (  pr$rotation[,1:2] , col='white' );
		text( pr$rotation[,1:2], labels= rownames(pr$rotation), cex=1.5 )
		dev.off()
		write.table( cbind( Genes = rownames(pr$rotation), pr$rotation[,1:2] ), file='gene_loadings.xls' , row.names=F, sep='\t',quote=F )
		#	mds.trans <- prcomp(t(tab))$x[,1:3]
	} else if ( mds.type == "LLE"){
		mds.proj <- LLE( tab, dim = 3, k = as.numeric(LLEK) )
		#	mds.trans <- LLE( t(tab), dim = 3, k = as.numeric(LLEK) )
	}else if ( mds.type == "ISOMAP"){
		mds.proj <- Isomap( tab, dim = 3, k = as.numeric(LLEK) )$dim3
		#	mds.trans <- Isomap( t(tab), dim = 3, k = as.numeric(LLEK) )$dim3
	}
	else {
		print( paste("Sory I can not work on the option",mds.type) )
	}
	
	geneC <- NULL
	if ( exists('geneGroups') ) {
		geneC <- geneGroups$groupID
	}
	dataObj$mds.coord <- mds.proj
	dataObj$geneC <- geneC
	dataObj <- clusters ( dataObj, onwhat=onwhat, clusterby=clusterby, mds.proj =  mds.proj, groups.n = groups.n, ctype = ctype, cmethod=cmethod )
	
	dataObj
}

#' cclusters calculates all possible clusters on the dataset
#' supported are hclust cclust and tclust with there respective options
#' @param dataObj the ExpressionSet object
#' @param clusterby (raw or MDS)
#' 
clusters <- function(dataObj,clusterby="raw", mds.proj=NULL,groups.n = 3, ctype='hierarchical clust',onwhat="Expression", cmethod='ward.D', tc_restr="eigen", tc_alpha=0.05, tc_nstart=50, tc_iter.max=20, tc_restr.fact=20 ){
	## custering	
	clusters <- NULL
	hc <- NULL
	if(onwhat=="Expression"){
		tab <- dataObj$z$PCR
	} 
	else {
		print ( paste ( "I work on the FACS data!" ) )
		tab <- dataObj$FACS
	}
	if ( exists('userGroups') ) {
		clusters <- userGroups$groupID
	}else if(clusterby=="MDS"){
		if ( ctype=='hierarchical clust'){
			hc <- hclust(dist( mds.proj ),method = cmethod)
			clusters <- cutree(hc,k=groups.n)
		}else if (  ctype=='kmeans' ) {
			hc <- hclust(dist( mds.proj ),method = cmethod)
			clusters <- kmeans( mds.proj ,centers=groups.n)$cluster
		}
	}
	else{#...do mds on tab
		if ( ctype=='hierarchical clust'){
			hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
			clusters <- cutree(hc,k=groups.n)
		}else if (  ctype=='kmeans' ) {
			hc <- hclust(as.dist( 1- cor(t(tab), method='pearson') ),method = cmethod)
			clusters <- kmeans( mds.proj ,centers=groups.n)$cluster
		}
	}
	if ( ! exists('userGroups') ){
		png ( file='hc_checkup_main_clustering_function.png', width=1600, height=800 )
		plot ( hc);
		dev.off()
	}
	dataObj$clusters <- clusters
	dataObj$hc <- hc
#	dataObj <- reorder_on_correlation ( dataObj )
#	print (cbind ( old= as.vector(dataObj$oldclusters), new=as.vector(dataObj$clusters)))
	dataObj
}