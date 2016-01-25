
# this file contains all generic fnction for data export and ploting


# Create an ExpressionSet object (S3)
# This object is mainly used for subsetting of the data and plotting
# @param dat data frame or matrix containing all expression data
# @param Samples A sample description table
# @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
# @param name The name of this object is going to be used in the output of all plots and tables - make it specific
# @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
# @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
# @param outpath Where to store the output from the analysis
# @param annotation The annotation table from e.g. affymetrix csv data
# @param newOrder The samples column name for the new order (default 'Order')

createWorkingSet <- function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL ) {
	S <- Samples
	if ( ! is.null(Analysis) ){
		S <- Samples[which ( Samples$Analysis == Analysis ), ]
	}
	if ( ! is.null(usecol) ) {
		S <- Samples[which ( Samples[, usecol] == 1 ),]
	}
	if ( exists('filename',S) ) {
		ret <- dat[, as.vector(S$filename) ]
		annotation <- dat[, is.na(match( colnames(dat), as.vector(S$filename) ))==T ]
	}else{
		ret <- dat[, as.vector(S[,namecol]) ]
		annotation <- dat[, is.na(match( colnames(dat), as.vector(S[,namecol]) ))==T ]
	}
	
	if ( exists( 'Order', S)){
		ret <- ret[, order(S$Order)  ]
		S <- S[order(S$Order), ]
	}
	
	colnames(ret) <- make.names(forceAbsoluteUniqueSample ( as.vector(S[, namecol]) ))
	S$SampleName <- colnames(ret)
	if ( is.null(dim(annotation))){
		## this xcould be a problem... hope we at least have a vector
		if ( length(annotation) == nrow(ret)) {
			rownames(ret) <- annotation
		}
		else {
			stop ( "Sorry, please cbind the rownames to the expression values before creating this object!")
		}
	}else{
		rownames(ret) <- annotation[,namerow]
	}
	write.table (cbind(rownames(ret), ret ), file=paste(name, '_DataValues',".xls", sep=''), sep='\t',  row.names=F,quote=F )
	write.table (S, file=paste(name,'_Sample_Description', ".xls", sep=''), sep='\t',  row.names=F,quote=F )
	data <- list ( 'data' = data.frame(ret), samples = S, name= name, annotation = annotation, rownamescol = namerow )
	data$sampleNamesCol <- namecol
	class(data) <- 'ExpressionSet'
	data$batchRemoved=0
	
	if ( is.null(outpath)){
		outpath = pwd()
	}else if ( ! file.exists(outpath)){
		dir.create( outpath )
	}
	data$outpath <- outpath
	data
}

coexprees4gene <- function( x, gene=NULL, method='spearman', padjMethod='BH'  ){
	UseMethod('coexprees4gene', x)
}

coexprees4gene.ExpressionSet <- function( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) {
	ret <- NULL
	if ( ! is.null(gene) ){
		z <- as.vector( t(x$data[ gene[1] ,]) )
		pval <- vector( 'numeric', nrow(x$data))
		cor <- vector( 'numeric', nrow(x$data))
		for ( i in 1:nrow(x$data) ) {
			try( {	res <-  cor.test( z, as.vector(t(x$data[i,]), 'numeric') ,method=method)
			pval[i] <- res$p.value
			cor[i] <- res$estimate }, silent=T
			)
		}
		adj.p <- p.adjust(pval , method = padjMethod) #BH
		ret <- data.frame('GeneSymbol' = x$annotation[,geneNameCol], pval = pval, cor = cor, adj.p = adj.p )
		rownames(ret) <- rownames(x$data)
	}
	ret
}

corMat.Pvalues <-function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL , name='tmp' ) {
	UseMethod('corMat.Pvalues', x)
}


# 80.505 for a 355 x 355 matrix 65 deep

corMat.Pvalues.ExpressionSet <-function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	d <- reduce.Obj( x, rownames(x$data)[which( apply(x$data,1,sd) > sd_cut)], name =name )
	if ( ! is.null(groupCol) ){
		ret <- list()
		names <- unique(d$samples[,groupCol])
		for ( i in 1:length(names)) {
			a <- subset( d, column=groupCol, value=names[i], name= paste(d$name,names[i],sep='_'), mode='equals' )
			
			ret[[i]] = corMat( a, sd_cut= sd_cut,method=method, geneNameCol=geneNameCol )
		}
		names(ret) <- names
		ret
	}
	else {
		n = nrow(d$data)
		print ( paste(d$name,": I create a",n,'x', n,'matrix') )
		ret <- cor(t(d$data), method=method )
		colnames(ret) <- rownames(ret) <- forceAbsoluteUniqueSample( as.vector(d$annotation[,geneNameCol]) )
		
		ret
	}
}

corMat <-function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	UseMethod('corMat', x)
}

corMat.ExpressionSet <-function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	d <- reduce.Obj( x, rownames(x$data)[which( apply(x$data,1,sd) > sd_cut)], name = name )
	if ( ! is.null(groupCol) ){
		ret <- list()
		names <- unique(d$samples[,groupCol])
		for ( i in 1:length(names)) {
			a <- subset( d, column=groupCol, value=names[i], name= paste(d$name,names[i],sep='_'), mode='equals' )
			ret[[i]] = corMat( a, sd_cut= sd_cut,method=method, geneNameCol=geneNameCol )
		}
		names(ret) <- names
		ret
	}
	else {
		n = nrow(d$data)
		print ( paste(d$name,": I create a",n,'x', n,'matrix') )
		ret <- cor(t(d$data), method=method )
		colnames(ret) <- rownames(ret) <- forceAbsoluteUniqueSample( as.vector(d$annotation[,geneNameCol]) )

		ret
	}
}

cor2cytoscape <- function (M, file, cut=0.9 ){
	edges <- NULL
	if ( class(M) == 'list' ) {
		edges <- list()
		for( i in 1:length(M) ) {
			edges[[i]] <- cor2cytoscape( M[[i]], file= paste(file,names(M)[i],'.txt',sep='_'), cut=cut) 		
		}
		names(edges) <- names(M)
	}
	else {
		edges <- NULL
		n <- as.vector(rownames(M))
		diag(M) <- 0
		for (i in 1:nrow(M)) {
			for (j in which(M[i,i:ncol(M)] >cut | M[i,i:ncol(M)] < -cut ) ) {
				edges <- rbind(edges, c(n[i], n[j], M[i,j] ))	
			}
		}
		#edges <- matrix( edges, ncol=3 )
		edges <- as.data.frame(edges)
		colnames(edges) <- c('Start', 'End', 'rho')
		edges$Type <- as.vector(rep( 1, nrow(edges) ))
		edges$Type[which(as.numeric(as.vector(edges$rho)) < 0)] <- -1
		write.table( edges, file , sep=" ",quote=F, row.names=F )
	}
	invisible(edges)
}


melt.ExpressionSet <- function( x, groupcol='GroupName', colCol='GroupName', probeNames="Gene.Symbol" ) {
	ma  <- x$data[,order(x$samples[,groupcol] )]
	rownames(ma) <- forceAbsoluteUniqueSample(x$annotation[, probeNames] )
	melted <- melt( cbind(rownames(ma),ma) )
	x$samples <- x$samples[order(x$samples[,groupcol]),]
	if ( length( which ( melted[,2] == '') ) > 0 ){
		melted <- melted[ - which ( melted[,2] == ''),]
	}
	melted[,3] <- as.numeric(as.character(melted[,3]))
	grps <- NULL
	for ( i in as.vector(x$samples[,groupcol]) ){
		grps <- c( grps, rep( i, nrow(x$data)))
	}
	cgrps <- NULL
	for ( i in as.vector(x$samples[,colCol]) ){
		cgrps <- c( cgrps, rep( i, nrow(x$data)))
	}
	colnames(melted) <- c('ProbeName', 'SampleName', 'Expression')
	melted$Group <- grps
	melted$ColorGroup <- cgrps
	melted
}


## plot grouped probesets creates ONE plot for ONE group of probesets
## If you want a multi group plot create it qourself from the single ones.



plot.probeset <- function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol" ) {
	UseMethod('plot.probeset', x)
}

plot.probeset.ExpressionSet <- function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol" ) {
	if ( sum(is.na(match(probeset, rownames(x$data)))==F) == 0 ){
		probeset <- rownames(x$data)[match( probeset, x$annotation[,geneNameCol] ) ]
	}
	if ( length(probeset) == 0 ) {
		stop( "I could not find the probeset in the data object!")
	}
	print ( probeset )
	x <- reduce.Obj(x, probeset )
	if ( nrow(x$data) == 1 ) {
		melted <- melt( x )
		collaps <- '_points_'
		if (boxplot){
			collaps <-'_boxplot_'
		}
		if ( pdf ) {
			pdf( file=paste(x$outpath,x$name,collaps,x$annotation[,geneNameCol],"_expression.pdf",sep='') ,width=5,height=4)
		}else {
			devSVG( file=paste(x$outpath,x$name,collaps,x$annotation[,geneNameCol],"_expression.svg",sep='') ,width=5,height=4)
		}
		if ( boxplot ){
			p <- ggplot(subset(dat, Gene.Symbol == gene) ,aes(x=Group, y=Expression,color=Group)) + geom_boxplot(shape=19)+theme_bw()
		}else{
			p <-ggplot(g1,aes(x=Group, y=Expression,color=Group)) +geom_jitter(shape=19)+theme_bw()
		}
		print ( p )
		dev.off()
		invisible(p)
	}
	else {
		stop("You should rather use a heatmap to display more than one probeset")
	}
}

drop.samples <- function ( x, samplenames=NULL, name='dopped_samples' ){
	UseMethod('drop.samples', x)
}

drop.samples.ExpressionSet <- function ( x, samplenames=NULL, name='dopped_samples' ){
	if ( ! is.null(samplenames)){
		S <- NULL
		S <- x$samples[ is.na(match(x$samples$SampleName, samplenames  ) ) == T ,]
		print ( paste( "Dropping", length(samplenames), "samples (", paste( samplenames, collapse=", "),")") )
		x$data <- x$data[, as.vector(S[, 'SampleName' ])]
		x$samples <- S
		x$name <- name
		if ( exists(where=x, 'vsdFull')){
			x$vsdFull <- NULL
			x$cds <- NULL
		}
		if ( exists(where=x, 'stats')){
			x$stats <- NULL
		}
		if ( exists(where=x, 'sig_genes')){
			x$sig_genes <- NULL
		}
		colnames(x$data) <- forceAbsoluteUniqueSample ( as.vector(S[, x$sampleNamesCol ]) )
		x$samples[,'SampleName'] <- colnames(x$data)
	}
	x
}

subset <-  function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ){
	UseMethod('subset', x)
}

subset.ExpressionSet <- function( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ){
	
	S <- NULL
	if ( is.null(value)) {
		stop( "the value must not be NULL!")
	}
	switch( mode,
			'less' = S <- x$samples[which ( x$samples[,column] <=  value), ], 
			'more' = S <- x$samples[which ( x$samples[,column] > value ), ], 
			'onlyless' = S <- x$samples[which ( x$samples[,column]  < value ), ],
			'equals' = S <- x$samples[which ( x$samples[,column] ==  value), ]
	)
	
	x$data <- x$data[, as.vector(S[, 'SampleName' ])]
	x$samples <- S
	x$name <- name
	if ( exists(where=x, 'vsdFull')){
		x$vsdFull <- NULL
		x$cds <- NULL
	}
	if ( exists(where=x, 'stats')){
		x$stats <- NULL
	}
	if ( exists(where=x, 'sig_genes')){
		x$sig_genes <- NULL
	}
	colnames(x$data) <- forceAbsoluteUniqueSample ( as.vector(S[, x$sampleNamesCol ]) )
	x$samples[,'SampleName'] <- colnames(x$data)
	write.table (S, file=paste(name,'_Sample_Description', ".xls", sep=''), sep='\t',  row.names=F,quote=F )
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

forceAbsoluteUniqueSample <- function ( x ,separator='_') {
	last = ''
	ret <- vector(length=length(x))
	for ( i in 1:length(x) ){
		if ( is.null(ret) ){
			last = x[i]
			ret[i] <- last
		}
		else{
			last = x[i]
			if ( ! is.na(match( last, ret )) ){
				last <- paste(last,separator,sum( ! is.na(match( x[1:i], last )))-1, sep = '')
			}
			ret[i] <- last
		}
	}
	ret
}

addAnnotation <- function(x ,mart, mart.col='refseq_mrna' ) {
	UseMethod('addAnnotation', x)
}

addAnnotation.NGSexpressionSet <- function(x ,mart, mart.col='refseq_mrna'){
	if ( ! class(mart) == 'data.frame' ){
		x$annotation <- cbind(x$annotation, mart[match(rownames(x$data),mart[,mart.col] ), ] )
	}
	else {
		x$annotation <- mart[is.na(match(mart[,mart.col], rownames(x$data)))==F,]
		rownames(x$annotation) <- rownames(x$data)
	}
	x
}

## get a vator from the annotation dataset
getAnnotation4probesets <- function (x, probesets=c(), colname='Gene.Symbol' ) {
	UseMethod('getAnnotation4probesets', x)
}
getAnnotation4probesets <- function (x, probesets=c(), colname='Gene.Symbol' ) {
	as.vector(x$annotation[match( probesets, rownames(x$data) ), colname ])
}

rank <- function( x ){
	UseMethod('rank', x)
}
rank.ExpressionSet <- function(x ) {
	if ( ! exists ( 'ranks', where =x ) ){
		x$ranks <- apply( x$data,2,order)
		colnames( x$ranks ) <- colnames(x$data) 
		rownames( x$ranks ) <- rownames(x$data) 
	}
	x
}

reduce.Obj <- function  ( x, probeSets=c(), name="reducedSet" ) {
	UseMethod('reduce.Obj', x)
}
reduce.Obj.ExpressionSet <- function ( x, probeSets=c(), name="reducedSet" ) {
	retObj <- x
	useOnly <- match(probeSets, rownames(x$data))
	not.matched <- probeSets[is.na(useOnly)]
	if ( length(not.matched) > 0 ){
		print (paste('Problematic genes:', paste(not.matched,sep=', ')))
		probeSets <- probeSets[ ! is.na(useOnly)]
		useOnly <- useOnly[ ! is.na(useOnly) ]
	}
	retObj$data <- data.frame( x$data[ useOnly ,] )
	rownames(retObj$data) <- probeSets
	colnames(retObj$data) <- colnames(x$data)
	retObj$name = name
	retObj$annotation <- x$annotation[useOnly,]
	if ( exists( 'stats', where=retObj) ){
		for ( i in 1:length(names(retObj$stats))){
			retObj$stats[[i]]= x$stats[[i]][ match(probeSets ,x$stats[[i]][,1] ),]
		}
	}
	if ( exists ( 'ranks', where=x)){
		retObj$ranks = x$ranks[useOnly,]
	}
	if ( exists ( 'raw', where=x)){
		retObj$raw = x$raw[useOnly,]
	}
	retObj[ match( c( 'design', 'cont.matrix', 'eset','fit', 'fit2' ), names(retObj) )] <- NULL
	retObj
}

### plotting data

ggplot.gene <- function (dat,gene, colrs, groupCol='GroupID', colCol='GroupID', boxplot=F) {
	UseMethod('ggplot.gene', dat)
}
ggplot.gene.ExpressionSet <- function(dat,gene, colrs, groupCol='GroupID', colCol='GroupID', boxplot=F){
	not.in = 'NUKL'
	g1 <- melt(reduce.Obj ( dat, gene, name=gene ), probeNames=isect$rownamescol, groupcol=groupCol,colCol=colCol)
	#g1 <- subset(dat, Gene.Symbol == gene)
	colnames(g1) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
	g1$Group <- factor( g1$Group, levels=unique( as.character(g1$Group)))
	if ( length(g1) == 0){
		not.in = c( gene )
	}
	
	if ( boxplot ){
		list ( plot=ggplot(g1 ,aes(x=Group, y=Expression,color=Group)) 
						+ geom_boxplot(shape=19)+theme_bw() 
						+ scale_colour_manual( values= colrs, guide=FALSE )
						+ theme( axis.text.x= element_text(angle=90) )
						+ labs(title=gene, x=''),
				not.in = not.in
		)
		
	}else{
		list ( plot=ggplot(g1,aes(x=Group, y=Expression,color=Group)) 
						+ geom_jitter(shape=19)+theme_bw()
						+ scale_colour_manual( values= colrs, guide=FALSE )
						+ theme( axis.text.x= element_text(angle=90) )
						+ labs(title=gene, x=''),
				not.in = not.in
		)
	}
}

gg.heatmap.list <- function (dat,glist, colrs, groupCol='GroupID', colCol='GroupID') {
	UseMethod('gg.heatmap.list', dat)
}
gg.heatmap.list.ExpressionSet <- function(dat,glist, colrs, groupCol='GroupID', colCol='GroupID'){
	
	isect <- reduce.Obj ( dat, glist)
	#browser()
	dat.ss <- melt ( isect, probeNames=isect$rownamescol, groupcol=groupCol,colCol=colCol)
	#dat.ss <- dat[which(is.na(match(dat$Gene.Symbol,isect))==F),]
	colnames(dat.ss) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
	dat.ss$z <- ave(dat.ss$Expression, dat.ss$Gene.Symbol, 
			FUN= function (x) { 
				n <- which(x==0)
				if ( length(x) - length(n) > 1 ){
					x[-n] <- scale(as.vector(t(x[-n])))
				}
				else {
					x[] = -20
				}
				x[n] <- -20
				x
			}
	)
	
	#dat.ss$z[which(dat.ss$z < -5)] <- -5
	#dat.ss$z[which(dat.ss$z > 5)] <- 5
	samp.cast <- dcast(dat.ss,Gene.Symbol~Sample,mean,value.var="z")
	samp.mat <- as.matrix(samp.cast[,2:ncol(samp.cast)])
	ord.genes <-
			as.vector(samp.cast[hclust(dist(samp.mat),method="ward.D")$order,1])
	dat.ss$Gene.Symbol <- with(dat.ss,factor(Gene.Symbol,levels =
							unique(as.character(ord.genes))))
	dat.ss$Sample <- with(dat.ss,factor(Sample,levels =
							unique(as.character(Sample))))
	dat.ss$Group <- with(dat.ss,factor(Group,levels =
							unique(as.character(Group))))
	dat.ss$colrss <- colrs[dat.ss$Group]
	ss <-dat.ss[which(dat.ss$Gene.Symbol==dat.ss$Gene.Symbol[1]),]
	brks= c( -20.1, quantile(dat.ss$z[which(dat.ss$z != -20)],seq(0,1,by=0.1)) )
	brks[length(brks)] = brks[length(brks)] + 0.1
	dat.ss$z <- cut( dat.ss$z, breaks= brks)
	
	list ( plot = ggplot(dat.ss, aes(x=Sample,y=Gene.Symbol))
					+ geom_tile(aes(fill=z))
					+ scale_fill_manual( values = c( 'gray', bluered(10)) ) 
					+ theme(
							legend.position= 'bottom',
							axis.text.x=element_blank(),
							axis.ticks.x=element_line(color=ss$colrss),
							axis.ticks.length=unit(0.6,"cm")
					)
					+ labs( y=''),
			not.in = setdiff( glist, rownames(isect$data)) )
}

plot.heatmaps <- function ( dataOBJ, gene.names , pvalue=1, analysis_name ='Unnamed', gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F ,... ) {
	
	UseMethod('plot.heatmaps', dataOBJ)
}


plot.heatmaps.ExpressionSet <- function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F,... ) {
	dataOBJ <- normalize(dataOBJ)
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
		geneNamesBased <- data.frame( 'GeneSymbol' = dataOBJ$annotation[rownames(exprs),geneNameCol], exprs )
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
	}else {
		collaps= ''
	}
	
	write.table( data.frame ( 'GeneSymbol' = rownames(geneNamesBased),geneNamesBased[,-1]),file= paste(fname,collaps,'_HeatmapID_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table(  geneNamesBased, file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table( dataOBJ$samples, file= paste(fname,collaps,'_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	
	if ( nrow(geneNamesBased) > 1 ){
		geneNamesBased <- data.frame(read.delim( file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''), header=T ))
		rownames(geneNamesBased) <- rn
		s <- ceiling(20 * nrow(geneNamesBased)/230 )
		if ( s < 5 ) {s <- 5}
		print (paste ("Height =",s))
		if ( pdf ) {
			pdf( file=paste(dataOBJ$outpath,dataOBJ$name,collaps,"_Heatmap_IDs_",pvalue,".pdf",sep='') ,width=10,height=s)
		}else {
			devSVG( file=paste(dataOBJ$outpath,dataOBJ$name,collaps,"_Heatmap_IDs_",pvalue,".svg",sep='') ,width=10,height=s)
		}
		print ( paste ("I create the figure file ",dataOBJ$outpath,dataOBJ$name,collaps,"_Heatmap_IDs_",pvalue,".svg",sep=''))
		heatmap.2(as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", margin=c(10, 10),dendrogram="row",scale="row",
				distfun=function (x) as.dist( 1- cor(t(x), method='pearson')),Colv=F, main=paste( analysis_name, pvalue ),...)	
		dev.off()
		
		rownames(geneNamesBased) <- paste(geneNamesBased[,1] , rownames(geneNamesBased) )
		if ( pdf ) {
			pdf( file=paste(dataOBJ$outpath,dataOBJ$name,collaps,"_Heatmap_GeneSymbols_",pvalue,".pdf",sep='') ,width=10,height=s)
		}else{
			devSVG( file=paste(dataOBJ$outpath,dataOBJ$name,collaps,"_Heatmap_GeneSymbols_",pvalue,".svg",sep='') ,width=10,height=s)
		}
		
		print ( paste ("I create the figure file ",dataOBJ$outpath,dataOBJ$name,collaps,"_Heatmap_GeneSymbols_",pvalue,".svg",sep=''))
		heatmap.2( as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", 
				margin=c(10, 12),dendrogram="row",scale="row",distfun=function (x) as.dist( 1- cor(t(x), method='pearson'),...),
				Colv=F, main=paste( analysis_name, pvalue ))
		dev.off()
	}
	else {
		print ( "Problems - P value cutoff to stringent - no selected genes!" );
	}
}

groups.boxplot <- function( x, SampleCol='GroupName', clusters=NULL, svg=F, fname='group_', width=800, height=800, mar=NULL, Collapse=NULL, ... ) {
	UseMethod('groups.boxplot',x)
}

#, las=2, cex.axis=2.5, cex.main=3, lwd=2, cex.lab=2, mar=c( 6,4.1,4.1,2.1 )
groups.boxplot.ExpressionSet <- function( x, SampleCol='GroupName', clusters, svg=F, fname='group_', width=800, height=800,mar=NULL, Collapse=NULL, ...) {
	maxG <- max( clusters )
	ret <- list()
	r=1
	gnames <- unique(as.vector(x$samples[,SampleCol]))
	if ( ! is.null( Collapse) ){
		x <- collaps( x, what = Collapse )
		fname = paste( fname, Collapse,"_",sep='')
	}
	x <- z.score(x)
	fnames <- vector('character', maxG )
	ma = -100
	mi = +100
	for ( i in 1:maxG ){
		if ( svg ) {
			fnames[i] = paste(x$outpath,fname,i,"_boxplot.C.svg",sep='')
			devSVG ( file= fnames[i], width=width/130, height=height/130 )
		}else{
			fnames[i] =paste(x$outpath,fname,i,"_boxplot.png",sep='')
			png( file=fnames[i],width=width,height=height )
		}
		
		robj <- reduce.Obj( x, names(clusters)[which(clusters==i)], name=paste("group_",i,sep='') )
		
		ret[[r]] <- rownames(robj$data)
		r = r+1
		names(ret)[length(ret)] = robj$name
		a= 1;
		d <- list()
		for ( n in as.vector(gnames) ){
			d[[a]] = as.vector(unlist(robj$data[,which(robj$samples[,SampleCol] == n)]))
			a = a+1
		}
		names(d) <- gnames
		A <- NULL
		if ( ! is.null(mar) ){
			A <- boxplot(d, main=paste('Cluster', i, ' - ',nrow(robj$data)," genes", sep='' ),outline=FALSE, par(mar=mar), ...  )
		}else {
			A <- boxplot(d, main=paste('Cluster', i, ' - ',nrow(robj$data)," genes", sep='' ),outline=FALSE, ...  )
		}
		mi <- min(c(mi,A$stats[1,]))
		ma <- max(c(ma, A$stats[5,]))
		dev.off()
	}
	print ( paste(  "min",mi, 'max', ma) )
	#print (paste('montage', paste(fnames, collapse= " "), "-geometry", paste("'", width, "x", height,"+0+0'",sep=''),paste(x$outpath,fname,"montage.png",sep=''), sep=' ' ))
	try( file.remove(  paste(x$outpath,fname,"montage.png",sep='') ) , silent=T )
	system ( paste('montage', fnames, "-geometry", paste("'", width, "x", height,"+0+0'",sep=''),paste(x$outpath,fname,"montage.png",sep='')," 2>/dev/null", collapse=' ' ) )
	ret
}

## obtained from https://stackoverflow.com/questions/8261590/write-list-to-a-text-file-preserving-names-r
write.list <- function( x, fname, sep=' ') {
	z <- deparse(substitute(x))
	cat(z, "\n", file=fname)
	nams=names(x) 
	for (i in seq_along(x) ){ cat(nams[i], sep,  x[[i]], "\n", 
				file=fname, append=TRUE) 
	}
}

simpleAnova <- function(x, samples.col='GroupName', padjMethod='BH' ) {
	UseMethod ('simpleAnova', x)
}

simpleAnova.ExpressionSet <- function ( x, samples.col='GroupName', padjMethod='BH' ) {
	x <- normalize(x)
	significants <- apply ( x$data ,1, function(x) { anova( lm (x ~ Samples[, samples.col]))$"Pr(>F)"[1] } )
	adj.p <- p.adjust( significants, method = padjMethod)
	res <- cbind(significants,adj.p )
	res <- data.frame(cbind( rownames(res), res ))
	colnames(res) <- c('genes', 'pvalue', paste('padj',padjMethod) )
	res <- list ( 'simpleAnova' = res )
	if ( exists( 'stats', where=x )) {
		x$stats <- c( x$stats, res)
	}else {
		x$stats <- res
	}
	x
}

## constructors that have to be implemented in the classes

createStats <- function ( x, condition, files=F, A=NULL, B=NULL) {
	UseMethod('createStats', x)
}
createStats.ExpressionSet <- function ( x, condition, files=F, A=NULL, B=NULL ) {
	stop( 'Not implemented' )
}

normalize <- function ( x ){
	UseMethod('normalize',x )
}

normalize.ExpressionSet <- function ( x ) {
	x
}

force.numeric <- function (dataObj ) {
	UseMethod('force.numeric', dataObj)
}
force.numeric.ExpressionSet <- function ( dataObj ) {
	for ( i in 1: ncol(dataObj$data) ) { 
		if ( !  paste(as.vector(dataObj$data[,i]), collapse=" ") == paste(as.vector(as.numeric(dataObj$data[,i])), collapse=" ") ) { 
			dataObj$data[,i] <- as.vector(as.numeric(dataObj$data[,i]))
		}
	}
	dataObj
}

### exporting data

print.ExpressionSet <- function (x) {
	cat (paste("An object of class", class(x)),"\n" )
	cat("named ",x$name,"\n")
	cat (paste( 'with',nrow(x$data),'genes and', ncol(x$data),' samples.'),"\n")
	cat (paste("Annotation datasets (",paste(dim(x$annotation),collapse=','),"): '",paste( colnames(x$annotation ), collapse="', '"),"'  ",sep='' ),"\n")
	cat (paste("Sample annotation (",paste(dim(x$samples),collapse=','),"): '",paste( colnames(x$samples ), collapse="', '"),"'  ",sep='' ),"\n")
	if ( exists(where=x, 'vsdFull')){
		cat ( "expression values are preprocessed\n" )
	}
	if ( exists(where=x, 'stats')){
		cat ( "P values were calculated for ", length(names(x$stats)) -1, " condition(s)\n")
	}
	if ( exists(where=x, 'sig_genes')){
		cat ( "significant genes were accessed for ", length(names(x$sig_genes)), " p value(s)\n")
		cat ( paste( names(x$sig_genes)),"\n")
	}
}


# write.data
# @param x the NGSexpressionSet

write.data <- function ( x, annotation=NULL ) {
	UseMethod('write.data', x)
}
write.data.ExpressionSet <- function ( x, annotation=NULL ) {
	if ( !is.null(annotation) ) {
		write.table( cbind( x$annotation[,annotation], x$data), file= paste( x$outpath,x$name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
	else {
		write.table( cbind( rownames(x$data), x$data), file= paste( x$outpath,x$name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
}

# export the statistic files
# @param x the NGSexpressionSet

writeStatTables <- function ( x, annotation=NULL, wData=F ) {
	UseMethod('writeStatTables', x)
}
writeStatTables.ExpressionSet <- function( x, annotation=NULL, wData=F ) {
	opath = paste(x$outpath,'stats_Pval/', sep='')
	if ( wData ) {
		opath = paste( x$outpath,'stats_wData/',sep='')
	}
	write.table( x$samples, file= paste(opath,'stats_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	Comparisons <- make.names(names(x$stats))
	system ( paste('mkdir ',opath, sep='') )
	for ( i in 1:length(Comparisons )){
		fname <- paste(opath ,Comparisons[i], '.xls',sep='')
		if ( ! is.null(annotation) ) {
			if ( wData==F ) {
				write.table( cbind( x$annotaion[, annotation],  x$stats[[i]] ) , file= fname, row.names=F, sep="\t",quote=F )
			}else {
				write.table( cbind(x$annotation[,annotation], x$stats[[i]], x$data ), file= fname, row.names=F, sep="\t",quote=F )
			}
		}
		else {
			if ( wData==F ) {
				write.table( x$stats[[i]], file= fname, row.names=F, sep="\t",quote=F )
			}else {
				write.table( cbind(x$stats[[i]], x$data ), file= fname, row.names=F, sep="\t",quote=F )
			}
		}
		print ( paste ( "table ",fname," for cmp ",Comparisons[i],"written" ) )
	}
}

export4GEDI <- function( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL ) {
	UseMethod('export4GEDI', x)
}

export4GEDI.ExpressionSet <- function( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) {
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
		stop ( "Sorry, please calculate the mean expression for the data first ('collaps()')") 
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
