library(affy)  # Affymetrix pre-processing
library(limma)
library(RSvgDevice)
library(gplots)
library(rgl)
library(stringr)
library(stats)


# Creaste a StefansStefansExpressionSet object (S3)
# This object is mainly used for subsetting of the data and plotting
# @param dat A LIMMA eset object
# @param Samples A sample description table
# @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
# @param name The name of this object is going to be used in the output of all plots and tables - make it specific
# @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
# @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
# @param outpath Where to store the output from the analysis
# @param annotation The annotation table from e.g. affymetrix csv data
# @param newOrder The samples column name for the new order (default 'Order')

ArrayStefansExpressionSet <- function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'Probe.Set.ID', usecol='Use', outpath = NULL, annotation=NULL, newOrder='Order' ) {
	S <- Samples
	if ( ! is.null(Analysis) ){
		if ( length(grep(usecol, colnames(Samples))) > 0 ){
			S <- Samples[which ( Samples$Analysis == Analysis & Samples[, usecol] == 1 ),]
		}else{
			S <- Samples[which ( Samples$Analysis == Analysis ),]
		}
	}
	ret <- exprs(dat)
	ret <- ret[, as.vector(S$filename) ]
	if ( exists(newOrder, where=S) ){
		ret <- ret[, order(S[,newOrder])  ]
		S <- S[order(S[,newOrder]), ]
	}
	colnames(ret) <- forceAbsoluteUniqueSample ( as.vector(S[, namecol]) )
	S$SampleName <- colnames(ret)
	annotation <- 	annotation[ match( rownames(ret), annotation[, namerow ] ), ]
	write.table (cbind(rownames(ret), ret ), file=paste(name, '_DataValues',".xls", sep=''), sep='\t',  row.names=F,quote=F )
	write.table (S, file=paste(name,'_Sample_Description', ".xls", sep=''), sep='\t',  row.names=F,quote=F )
	if ( is.null(outpath)){
		outpath = pwd()
	}else if ( ! file.exists(outpath)){
		dir.create( outpath )
	}
	data <- list ( eset=dat, 'data' = ret, samples = S, name= name, annotation = annotation, rownamescol = namerow, sampleNamesCol =  namecol, outpath = outpath )
	data$genes <- rownames(data$data)
	class(data) <- 'StefansStefansExpressionSet'
	data
}


addAnnotation <- function(x ,mart, mart.col='refseq_mrna' ) {
	UseMethod('addAnnotation', x)
}

addAnnotation.StefansStefansExpressionSet <- function(x ,mart, mart.col='refseq_mrna'){
	if ( ! class(mart) == 'data.frame' ){
		x$annotation <- cbind(x$annotation, mart[match(rownames(x$data),mart[,mart.col] ), ] )
	}
	else {
		x$annotation <- mart[is.na(match(mart[,mart.col], rownames(x$data)))==F,]
		rownames(x$annotation) <- rownames(x$data)
	}
	x
}
pwd <- function () {
	system( 'pwd > __pwd' )
	t <- read.delim( file = '__pwd', header=F)
	t <- as.vector(t[1,1])
	t <- paste(t,"/",sep='')
	unlink( '__pwd')
	t
}


# This function will use the LIMMA package to calculate statistics.
# @param x The StefansStefansExpressionSet object
# @param fcol The column in the samples table, that contains the group definition data
# @param perlscript The all_vs_all_model_matrix.pl perl script that creates the cont.matrix and Comparisons data
# @param files Whether to export the stat tables; default TRUE

#createStats <- function(x, fcol=NULL, perlscript= '/home/slang/workspace/BMC_NGS/bin/all_vs_all_model_matrix.pl', files=T) {
#	UseMethod('createStats', x)
#}

createStats.StefansStefansExpressionSet <- function (x, fcol=NULL, perlscript= '/home/slang/workspace/BMC_NGS/bin/all_vs_all_model_matrix.pl', files=T) {
	if ( length( grep(fcol, colnames(x$samples)) ) == 0 ){
		stop( paste("The sample description",fcol,"can not be found in the Samples table" ) )
	}
	Factors <- as.factor( x$samples[,fcol] )
	x$design <- model.matrix( ~ -1+Factors)
	colnames(x$design) <- str_replace_all( colnames(x$design), 'Factors', '' )
	if ( ! file.exists(perlscript) ){
		stop ("The required all_vs_all_model_matrix.pl perl script is missing - please install the BMC_NGS perl package.")
	}
	print ( paste ('perl ', perlscript, ' -where x -groups ', paste( colnames(x$design), collapse=' '), ' > ',x$outpath,'Groups.R',sep='' ) )
	system( paste ('perl ', perlscript, ' -where x -groups ', paste( colnames(x$design), collapse=' '), ' > ',x$outpath,'Groups.R',sep='' ) )
	print (paste( x$outpath,'Groups.R',sep='' ) )
	source (paste( x$outpath,'Groups.R',sep='' ) , local=TRUE )
	x$fit <- lm.series(x$eset, x$design)
	x$fit2 <- contrasts.fit(x$fit,x$cont.matrix)
	x <- createStatTable( x , files=files)
	x
}

remBatchEffect<- function(x, batch=c(), fcol=NULL , design=FALSE  ) {
	UseMethod('remBatchEffect', x)
}
remBatchEffect.StefansStefansExpressionSet <- function(x, batch=c(), fcol=NULL , design=FALSE ) {
	if ( exists ('BatchEffectRemodved', x) ){
		if ( x$BatchEffectRemodved) {stop ("Batch effect has already been removed!" ) }
	}
	if ( design ) {
	if ( ! exists( 'design', where=x) ){
		if ( length( grep(fcol, colnames(x$samples)) ) == 0 ){
			stop( paste("The sample description",fcol,"can not be found in the Samples table" ) )
		}
		Factors <- as.factor( x$samples[,fcol] )
		x$design <- model.matrix( ~ -1+Factors)
		colnames(x$design) <- str_replace_all( colnames(x$design), 'Factors', '' )
	}
	x$data <- removeBatchEffect( x$data, batch=batch, design=x$design )
	}else{
		x$data <- removeBatchEffect( x$data, batch=batch )
	}
	x$BatchEffectRemodved = TRUE
	x
}



createStats.StefansStefansExpressionSet <- function ( x, condition, files=T , A=NULL, B=NULL) {
	if ( length( grep ( condition, colnames(x$samples))) > 0 ) {
		condition = factor( x$samples[,condition] )
	}
	x$toptabs <- vector("list", length(x$Comparisons) )
	for ( i in 1:length(x$Comparisons) ) {
		x$toptabs[[i]] <- toptable(x$fit2, coef = i, adjust='fdr',number=1000000, genelist = x$annotation)
	}
	names(x$toptabs) <- x$Comparisons
	if ( files ){
		writeStatTables( x )
	}
	x
}

meanExpression4groups <- function ( x,  fcol=NULL ) {
	UseMethod('meanExpression4groups', x)
}
# meanExpression4groups calculates the mean expression for groups defined in the samples table
meanExpression4groups.StefansStefansExpressionSet <- function ( x, fcol=NULL ) {
	stop = FALSE
	if ( is.null (fcol) ){
		stop = TRUE
	}else if ( length(match( fcol, colnames(x$samples) ) ) == 0){
		stop = TRUE
	}
	if ( stop)	{stop('I need a grouping information column name in the samples table in order to built the groups!' )	}
	mdata <- x
	cols = unique(x$samples[,fcol])
	mdata$data <- matrix( 0, ncol=length(cols), nrow=nrow(mdata$data))
	colnames(mdata$data) <- cols
	mdata$oldSamples <- mdata$samples
	old_lines <- NULL
	pos1 <- vector('numeric', length(cols))
	for ( i in cols ) { 
		mdata$data[,i] <- as.vector(apply( x$data[,grep(i, x$samples[,fcol])],1,mean))
		old_lines[i] <- paste( grep(i, x$samples[,fcol]), collapse=" " )
		pos1[i] <- grep(i, x$samples[,fcol])[1]
	}
	mdata$samples <- mdata$samples[pos1,]
	mdata$samples$oldLines <- old_lines
	colnames(mdata$data) <- cols
	rownames(mdata$data) <- rownames(x$data)
	mdata
}

export4GEDI <- function( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL ) {
	UseMethod('export4GEDI', x)
}

export4GEDI.StefansStefansExpressionSet <- function( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) {
	if ( is.null(time.col)){
		stop ( paste( "choose time.col from:", paste( colnames(x$samples), collapse=", ") ) ) 
	}
	if ( is.null(tag.col)){
		stop ( paste( "choose tag.col from:", paste( colnames(x$samples), collapse=", ") ) ) 
	}
	groupnames <- vector('numeric', nrow(x$samples))
	for (i in 1:nrow(x$samples)) {
		groupnames[i] = paste( x$samples[i, tag.col], x$samples[i, time.col] , sep="_" )
	}
	if ( length(which(table(groupnames) > 1)) > 0 ){
		stop ( "Sorry, please calculate the mean expression for the data first ('meanExpression4groups')") 
	}
	

	treatments <- unique( as.vector(x$samples[ , tag.col] ))
	
	## now I need to find all possible days
	possible <- NULL
	for ( t in treatments ){
		possible <- c( possible, x$samples[ grep( t, x$samples[ , tag.col]) , time.col] )
	}
	required = names(which(table(possible) > minSample_PerTime) )
	passed <- NULL
	for ( t in treatments ){
		l <- rep ( -1, nrow(x$samples) )
		p <- grep( t, x$samples[ , tag.col])
		if ( length(p) >= length(required) ){
			l[p] <- x$samples[ p, time.col]
			l[ is.na(match(l, required))== T ] <- -1
			passed <- c(passed, paste( c(paste("}",t,sep=""),  l), collapse="\t") )
		}
	}
	
	fileConn<-file(fname)
	open( fileConn, "w" )
	writeLines( paste( "}Dynamic", length(passed), nrow(x$data), x$name, sep="\t"), con=fileConn)
	writeLines( paste(passed , collapse="\n") , con=fileConn )
	
	for ( i in 1:nrow(x$data) ){
		writeLines( paste( c(rownames(x$data)[i], x$data[i,]), collapse="\t" ),  con=fileConn )
	}
	
	close(fileConn)
	print( paste ("created file", fname))
}

getProbesetsFromStats <- function ( x, v=0.05, pos=6, mode='less', Comparisons=NULL ) {
	UseMethod('getProbesetsFromStats', x)
}


# getProbesetsFromStats returns a list of probesets (the rownames from the data matrix) for a restriction of a list of stat comparisons
# @param v The cutoff value
# @param pos The column in the stats tables to base the selection on
# @param Comparisons A list of comparisons to check (all if left out)
# @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
# @examples 
# probes <- getProbesetsFromStats ( x, v=1e-4, pos="adj.P.Val" ) 
# returns a list of probesets that shows an adjusted p value below 1e-4

getProbesetsFromStats.StefansStefansExpressionSet <- function ( x, v=0.05, pos=6, mode='less', Comparisons=NULL ) {
	if ( is.null(Comparisons)){	Comparisons <- names(x$toptabs) }
	probesets <- NULL
	for ( i in match(Comparisons, names(x$toptabs) ) ) {
		switch( mode,
			'less' = probesets <- c( probesets, as.vector(x$toptabs[[i]][which(x$toptabs[[i]][,pos] <= v),1] )),
			'more' = probesets <- c( probesets, as.vector(x$toptabs[[i]][which(x$toptabs[[i]][,pos] > v),1] )), 
			'onlyless' = probesets <- c( probesets, as.vector(x$toptabs[[i]][which(x$toptabs[[i]][,pos] < v),1] )),
			'equals' = probesets <- c( probesets, as.vector(x$toptabs[[i]][which(x$toptabs[[i]][,pos] == v),1] ))
		)
	}
	unique(probesets)
}

getProbesetsFromValues <- function ( x, v='NULL', sample='NULL', mode='less' ){
	UseMethod('getProbesetsFromValues', x)
}


# Select probesets, that show a certain level in expression for a single sample
# @param v The cutoff value
# @param sample The sample name
# @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
# @examples 
# probes <- getProbesetsFromStats ( x, v=10, sample="A" ) 
# returns a list of probesets that has a expression less than 10 in sample A

getProbesetsFromValues.StefansStefansExpressionSet <- function ( x, v='NULL', sample='NULL', mode='less' ){
	s <- FALSE
	if ( is.null(v) ){
		s<-TRUE
	}
	if ( is.null(sample) ){
		s<-TRUE
	}
	if ( s ) { stop ( "Please give me the required values for v and sample") }
	probesets <- NULL
	switch( mode,
			'less' = probesets <-  as.vector(rownames(x$data)[which(x$data[,sample] <= v)] ) ,
			'more' = probesets <- as.vector(rownames(x$data)[which(x$data[,sample] > v)] ), 
			'onlyless' = probesets <- as.vector(rownames(x$data)[which(x$data[,sample] < v)] ),
			'equals' = probesets <- as.vector(rownames(x$data)[which(x$data[,sample] == v)] )
	)
	probesets
}
reduce.Obj <- function  ( x, probeSets=c(), name="reducedSet" ) {
	UseMethod('reduce.Obj', x)
}
reduce.Obj.StefansStefansExpressionSet <- function ( x, probeSets=c(), name="reducedSet" ) {
	retObj <- x
	retObj$data <- matrix(x$data[match(probeSets, rownames(x$data)),],ncol=ncol(x$data) )
	rownames(retObj$data) <- probeSets
	colnames(retObj$data) <- colnames(x$data)
	retObj$genes <- probeSets
	retObj$name = name
	retObj$annotation <- retObj$annotation[match(probeSets, retObj$annotation[,retObj$rownamescol]),]
	if ( exists( 'toptabs', where=retObj) ){
		for ( i in 1:length(names(retObj$toptabs))){
			retObj$toptabs[[i]]= retObj$toptabs[[i]][ match(probeSets ,retObj$toptabs[[i]][,1] ),]
		}
	}
	retObj[ match( c( 'design', 'cont.matrix', 'eset','fit', 'fit2' ), names(retObj) )] <- NULL
	retObj
}

pca.rows <- function (x, groups=2, rad=0.01, info = 'no_more_info', nmax=3000, dname='genes', method='ward.D2', group_on_pca=FALSE )  {
	UseMethod('pca.rows', x)
}
pca.rows.StefansStefansExpressionSet <- function (x, groups=2, rad=0.01, info = 'no_more_info', nmax=3000, dname='genes', method='ward.D2', group_on_pca=FALSE )  {
	## plot the pca for the rows (genes)
	## the groups will be created from a hclust call
	ma <- x$data
	cols <- rainbow( groups )
	if ( ncol(ma) <= nmax ) {
		pc <- prcomp(ma)
		hc <- NULL
		n <- NULL
		if ( group_on_pca ){
			hc <- hclust(dist( pc$x ),method = method)
			n <- cutree(hc, groups)
		}else{
			hc <- hclust( as.dist( 1- cor(t(ma))), method= paste(method) )
			n <- cutree( hc, groups)
		}
		genes <- vector('list',groups)
		names(genes) <- c( paste( 'Group ', 1:groups) )
		for( i in 1:groups ) {
			rgl.spheres( pc$x[grep(i, n),1:3], col=cols[i],  radius=rad )
			genes[[i]] <- as.vector( rownames(ma)[grep(i, n)] )
			#print( paste ("Group ",i," in ",cols[i]," rep. ",length(genes[[i]]), ' ',dname, sep ='')) 
		}
		axes3d(edges = "bbox", col='black' )
		writeWebGL( width=800, height=800, dir = paste(outpath,"webGL",info, sep='') )
		
		lst <- list ( dist = hc, groups =  genes )
		ofile = paste(outpath,'PCA_analysis',info,sep='')
		genesTable <- matrix( rep ( 0 ,groups*3), ncol=3 )
		colnames(genesTable) <- c( "group", 'probesets', 'gene symbols' )
		groupTags <- vector(mode="character", length= groups)
		for ( a in 1:groups) {
			devSVG( file= paste(ofile, 'Group',a,'.svg',sep=''), width=6, height=6)
			groupTags[a] = paste('group ',a," (n=",length(lst$groups[[a]]),")", sep=' ')
			genesTable[a,] = c( 
					groupTags[a], 
					paste( lst$groups[[a]], collapse=" " ), 
					paste( x$names[match (lst$groups[[a]],x$genes)], collapse=" ") 
			)
			simpleLinePlot ( reduce.Obj( x, lst$groups[[a]]), mar= c(5,4,1,1), main=groupTags[a] )
			dev.off()
		}
		## plot the 3D legend in 2D
		png( file=paste(outpath,"webGL",info,"/legend.png", sep=''), width=400, height=800)
		plot( 1, type='n',xaxt='n', xlab='',ylab="", xlim=c(1,12), ylim=c(1,5) ,main= "Group coloring legend")
		legend( 'topleft', groupTags, fill=cols )
		dev.off()
		
		## print the table of genes used
		write.table(genesTable, file=  paste(outpath,'PCA_analysis_',info,'_gene_groups.xls',sep=''), sep='\t',  row.names=F,quote=F )
		w <- round ( nrow(ma) * 40 ) 
		if ( w > 2000 ) {w <- 2000}
		png ( file= paste(outpath,'Dendrogram_PCA_analysis',info,'.png',sep=''), height = 600, width = w )
		plot( lst$dist )
		dev.off()
		multiplot.4.n.figures( groups )
		for ( a in 1:groups) {
			simpleLinePlot ( reduce.Obj( x, lst$groups[[a]]), main=paste('group ',a," (n=",length(lst$groups[[a]]),")", sep=' ') )
		}
		##return the list containing the dist = hc and the groups = genes
		x$groups <- genes
		x$hc <- hc
		x$pc <- pc
		x
	}
	else {
		print( paste('Sorry too many selected genes',ncol(ma)))
		x
	}
}
simpleLinePlot <- function ( x, main="", mar=c(5.1, 4.1, 4.1, 2.1) ){
	UseMethod('simpleLinePlot', x)
}
simpleLinePlot.StefansStefansExpressionSet <- function ( x, main="", mar=c(5.1, 4.1, 4.1, 2.1) ){
	if ( length(x$genes) > 1 ){
		x <- z.score(x)
	}
	par(mar=mar)
	plot( 1, type='n',xaxt='n', xlab='', ylab="relative expression", xlim=c(1,ncol(x$data)), ylim=c(min(x$data),max(x$data)), main=main )
	xvals <- 1:ncol(x$data)
	axis(1, at=xvals, labels=colnames(x$data), cex.axis=0.4, las=2 )
	for ( i in 1:nrow(x$data) ) {
		lines ( xvals,x$data[i,], type='l' )
	}
	par(mar=c(5.1,4.1,4.1,2.1))
}

z.score <- function ( m ){
	UseMethod('z.score', m)
}
z.score.StefansStefansExpressionSet <- function(m) {
	rn <- rownames( m$data )
	me <- apply( m$data, 1, mean )
	sd <- apply( m$data, 1, sd )
	sd[which(sd==0)] <- 1e-8
	m$data <- (m$data - me) /sd
	rownames(m$data) <- rn
	m
}


multiplot.4.n.figures <- function ( n ) {
	ret <- NULL
	if ( n > 1) {
		x <- sqrt(n)
		if ( as.integer(x) == x ){
			par ( mfrow=c (x,x) )
			ret <- c( x,x)
		}
		else {
			x= as.integer(x)
			if ( as.integer(sqrt(n)+0.5) != x ){
				par ( mfrow=c (x+1,x+1) )
				ret <- c( x+1,x+1)
			}
			else {
				par ( mfrow=c (x+1,x) )
				ret <- c( x+1,x)
			}
		}
	}
	ret
}

plot.StefansStefansExpressionSet <- function ( x, pvalue=c( 0.1,1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), Subset=NULL , Subset.name= NULL, comp=NULL, gene_centered=F, collaps=NULL,geneNameCol= "mgi_symbol") {
	if ( !is.null(comp) ){
		print ("At the moment it is not possible to reduce the plotting to one comparison only" )
		return (x)
	}
	add = ''
	orig.name = x$name
	if ( gene_centered ) {
		add = 'GenesOnly'
	}
	#browser()
	if ( ! is.null(collaps)){
		if ( nchar(add) > 1 ){
			add = paste( 'mean', add, sep='_')
		}else{
			add = 'mean'
		}
	}
	if ( ! is.null(Subset) ){
		if ( is.null(Subset.name)){
			Subset.name = 'subset_name_unset'
		}
		if ( nchar(add) > 1 ){
			add = paste( add, Subset.name,sep='_')
		}else{
			add = Subset.name
		}
	}
	if ( nchar(add) > 0 ){
		x$name = paste( add,x$name, sep='_')
	}
	print ( x$name )
	for (i in match( names(x$sig_genes), pvalue ) ){
		try( plot.heatmaps( 
						x, 
						x$sig_genes[[i]],
						names(x$sig_genes)[i],
						analysis_name = paste (x$name, version), 
						gene_centered = gene_centered,
						Subset = Subset,
						collaps = collaps,
						geneNameCol= geneNameCol
				), silent=F)
	}
	x$name = orig.name
}


plot.heatmaps <- function ( dataOBJ, gene.names , pvalue=1, analysis_name ='Unnamed', gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol" ) {
	
	UseMethod('plot.heatmaps', dataOBJ)
}


plot.heatmaps.StefansStefansExpressionSet <- function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol" ) {	
	print ( "The function is untested - use with care!" )
	dataOBJ <- z.score( dataOBJ )
	
	if ( is.null(analysis_name)){
		analysis_name = dataOBJ$name
	}
	exprs = dataOBJ$data
	drop <- which( (apply(exprs,1,sd) == 0) == T)
	if ( length(drop) > 0 ){
		print ( paste("I need to drop",length(drop),"genes as the varianze is 0") )
		exprs = exprs[-drop,]
	}
	if ( ! is.null(gene.names)){
		geneNamesBased <- exprs[is.na(match(rownames(exprs),gene.names$all$id ))==F,]
		geneNamesBased <- data.frame ( 'GeneSymbol' = as.vector (dataOBJ$annotation[is.na(match ( dataOBJ$annotation$GeneID, gene.names$all$id ) ) == F, geneNameCol]),  geneNamesBased) 
	}
	else {
		geneNamesBased <- data.frame( 'GeneSymbol' = rownames(exprs), exprs )
	}
	if ( ! is.null(Subset) ) { ## subset is a list of srings matching to the gene symbols
		useful <- NULL
		for ( i in 1:length(Subset) ){
			useful <- c( useful, grep(Subset[i], geneNamesBased[,1] ) )
		}
		if ( length(useful) > 0 ){
			geneNamesBased <- geneNamesBased[ unique(useful), ]
		}
	}
	if ( gene_centered ){
		ok <- as.vector(which ( is.na(geneNamesBased[,1] )))
		grepped <- grep('_drop',forceAbsoluteUniqueSample(as.vector(geneNamesBased[,1]), '_drop_' ))
		if ( length(ok) > 1 ){
			geneNamesBased <- geneNamesBased[- grepped [ - match ( ok[-1] , grepped ) ], ]
		}
		else if ( length(grepped) > 0 ) {
			geneNamesBased <- geneNamesBased[-grepped, ]
		}
	}
	fname = paste(dataOBJ$outpath, 'GenesOnly_',dataOBJ$name,sep='')
	rn <- rownames(geneNamesBased)
	
	if ( ! is.null(collaps)){
		u <- unique(as.vector(dataOBJ$samples$GroupName))
		m <- length(u)
		mm <-  matrix ( rep(0,m * nrow(geneNamesBased)), ncol=m)
		colnames(mm) <- u
		rownames(mm) <- rownames(geneNamesBased)
		if ( collaps=="median") {
			for ( i in u ){
				## calculate the mean expression for each gene over all cells of the group
				mm[,i] <- apply( geneNamesBased[ , which(as.vector(dataOBJ$samples$GroupName) == i )+1],1,median)
			}
		}
		else {
			for ( i in u ){
				## calculate the mean expression for each gene over all cells of the group
				mm[,i] <- apply( geneNamesBased[ , which(as.vector(dataOBJ$samples$GroupName) == i )+1],1,mean)
			}
		}
		geneNamesBased <- data.frame( GeneSymbol = as.vector(geneNamesBased[,1]), mm )
	}
	write.table( data.frame ( 'GeneSymbol' = rownames(geneNamesBased),geneNamesBased[,-1]),file= paste(fname,'_HeatmapID_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table(  geneNamesBased, file= paste(fname,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	if ( nrow(geneNamesBased) > 1 ){
		geneNamesBased <- data.frame(read.delim( file= paste(fname,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''), header=T ))
		rownames(geneNamesBased) <- rn
		s <- ceiling(20 * nrow(geneNamesBased)/230 )
		if ( s < 5 ) {s <- 5}
		print (paste ("Height =",s))
		devSVG( file=paste(dataOBJ$outpath, dataOBJ$name,"_Heatmap_IDs_",pvalue,".svg",sep='') ,width=10,height=s)
		print ( paste ("I create the figure file ",fname,"_Heatmap_IDs_",pvalue,".svg",sep=''))
		heatmap.2(as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", margin=c(10, 10),dendrogram="row",scale="row",
				distfun=function (x) as.dist( 1- cor(t(x), method='pearson')),Colv=F, main=paste( analysis_name, pvalue ))	
		dev.off()
		rownames(geneNamesBased) <- paste(geneNamesBased[,1] , rownames(geneNamesBased) )
		devSVG( file=paste(ataOBJ$outpath, dataOBJ$name,"_Heatmap_GeneSymbols_",pvalue,".svg",sep='') ,width=10,height=s)
		print ( paste ("I create the figure file ",fname,"_Heatmap_GeneSymbols_",pvalue,".svg",sep=''))
		heatmap.2( as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", 
				margin=c(10, 10),dendrogram="row",scale="row",distfun=function (x) as.dist( 1- cor(t(x), method='pearson')),
				Colv=F, main=paste( analysis_name, pvalue ))
		dev.off()
	}
	else {
		print ( "Problems - P value cutoff to stringent - no gselected genes!" );
	}
}