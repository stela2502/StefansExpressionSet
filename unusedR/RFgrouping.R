# TODO: Add comment
# 
# Author: stefan
###############################################################################

create_RFgroupingObj <- function( geneG='RandomForestdistRFobject.RData', sampleG='RandomForestdistRFobject_genes.RData' ) {
	x <- list()
	load(geneG)
	x$geneG <- distRF
	rm(distRF)
	load(sampleG)
	x$sampleG <- distRF
	rm(distRF)
	class( x ) <- append(  'RFgrouping', class(x))
	x
}


createGeneGroups <- function (dataObj, groups=10 ) {
	UseMethod('createGeneGroups', dataObj)
}	

createGeneGroups.RFgrouping <- function (dataObj, groups=10 ) {

	persistingGenes <- rownames( dataObj$data )
	
	res = pamNew(dataObj$geneG$cl1, groups )
	N <- names( res )
	## probably some cells have been kicked in the meantime - I need to kick them too
	N <- intersect( persistingGenes, N )
	geneGroups <- matrix(ncol=3, nrow=0)
	for ( a in 1:length(N) ){
		geneGroups <- rbind (geneGroups, c( N[a], as.numeric(res[[N[a]]]) ) )
	}
	colnames(geneGroups) <- c('cellName', 'groupID' )
	## write this information into a file that can be used as group
	geneGroups = data.frame( geneGroups)
	save ( geneGroups , file= paste("forest_Gene_group_n", groups,'.RData', sep=''))
	dataObj$geneGroups <- geneGroups
	dataObj$geneGroupsCounts <- max(geneGroups[,2])
	
	dataObj
}

createSampleGroups <- function (dataObj, groups=10 ) {
	UseMethod('createSampleGroups', dataObj)
}	

createSampleGroups.RFgrouping <- function (dataObj, groups=10 ) {
	
	if ( round(length(persistingSamples)/4) < groups ){
		groups <- round(length(persistingSamples)/4)
	}
	if (groups < 2 ){
		groups <- 2
	}
	res = pamNew(dataObj$sampleG$cl1, groups )
	N <- names( res )
	## probably some cells have been kicked in the meantime - I need to kick them too
	N <- intersect( persistingSamples , N )
	geneGroups <- matrix(ncol=3, nrow=0)
	for ( a in 1:length(N) ){
		geneGroups <- rbind (geneGroups, c( N[a], 'no info', as.numeric(res[[N[a]]]) ) )
	}
	colnames(geneGroups) <- c('geneName', 'userInput',  'groupID' )
	## write this information into a file that can be used as group
	geneGroups = data.frame( geneGroups)
	save ( geneGroups , file= paste("forest_gene_group_n", groups,'.RData', sep=''))
	
	fileConn<-file(paste("Gene_grouping.randomForest.txt", sep="") )
	writeLines(c(paste("load('forest_gene_group_n",groups,".RData')",sep=""),
					"geneGroups <- checkGrouping ( geneGroups[is.na(match(geneGroups$geneName, colnames(data.filtered$PCR) ))==F, ] )",
					"write.table( geneGroups[order(geneGroups[,3]),], file='GeneClusters.xls' , row.names=F, sep='\t',quote=F )"
			), fileConn)
	close(fileConn)
	
}
