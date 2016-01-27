setClass( 
		Class='ExpressionSet', 
		representation = representation ( 
			data='data.frame',
			samples='data.frame',
			annotation='data.frame',
			outpath='character',
			name='character',
			rownamescol='character',
			sampleNamesCol='character',
			stats = 'list'
		)
)
#' this file contains all generic fnction for data export and ploting
#' Create an ExpressionSet object (S3)
#' This object is mainly used for subsetting of the data and plotting
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
#' @param outpath Where to store the output from the analysis
#' @param annotation The annotation table from e.g. affymetrix csv data
#' @param newOrder The samples column name for the new order (default 'Order')
#' @export 
setGeneric("ExpressionSet", ## Name
		function( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL){ ## Argumente der generischen Funktion
			standardGeneric("ExpressionSet") ## der Aufruf von standardGeneric sorgt für das Dispatching
		})

setMethod("ExpressionSet", signature = c ('ExpressionSet'), 
	definition = function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL ) {
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
	setClass ( Class = 'ExpressionSet', representation = representation ( data) )
	data
})


#' Calculate the coexpression for any given gene
#' 
#' @param x the ExpressionSet varibale
#' @param method any method supported by \link[=cor.test]{cor.test}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param padjMethod the method to calucate the FDR with \link[=p.adjust]{p.adjust}
#' 
#' @return a data.frame with the columns 'GeneSymbol', 'pval', 'cor', 'adj.p'
#' @export 
setGeneric('coexprees4gene', ## Name
	function ( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) { ## Argumente der generischen Funktion
		standardGeneric('coexprees4gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('coexprees4gene', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) {
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
})


#' calculate a p_value matrix for the data object
#' 
#' @param x the ExpressionSet varibale
#' @param sd_cut the cut off value for the sd check
#' @param method any method supported by \link[=cor.test]{cor.test}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param groupCol the column name of the grouping variable in the samples table
#' @param name the name of the analysis
#' 
#' THIS FUNCTION IS NOT DOING THE RIGTH THING
#' BROKEN
#' 
#' @export 
setGeneric('corMat.Pvalues', ## Name
	function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) { ## Argumente der generischen Funktion
		standardGeneric('corMat.Pvalues') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('corMat.Pvalues', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	# TODO: implement the p value calculation!
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
})



#' calculate a correlation matrix for the data object
#' 
#' @param x the ExpressionSet varibale
#' @param sd_cut the cut off value for the sd check
#' @param method any method supported by \link[=cor.test]{cor.test}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param groupCol the column name of the grouping variable in the samples table
#' @param name the name of the analysis
#' 
#' @return the correlation matrix
#' @export 
setGeneric('corMat', ## Name
	function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) { ## Argumente der generischen Funktion
		standardGeneric('corMat') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('corMat', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
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
})
#' exports the correlation matrix in the format \href{http://www.cytoscape.org/}{CytoScape} 
#' does support as network file.
#' 
#' @param M the correlation matrix obtained by a run of \code{\link{corMat}}
#' @param file the outfile
#' @param cut the cut of value for the correlation rho value (0.9)
#' @export 
setGeneric('cor2cytoscape', ## Name
	function (M, file, cut=0.9 ) { ## Argumente der generischen Funktion
		standardGeneric('cor2cytoscape') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('cor2cytoscape', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (M, file, cut=0.9 ) {
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
})

#' melts the object using  \code{\link[reshape2]{melt}}
#' 
#' @param x the ExpressionSet object
#' @param groupcol the grouping column in the samples data
#' @param colCol the coloring column in the sample data
#' @param probeNames the column in the annotation datacontaining the gene symbol
#' @export 
setGeneric('melt', ## Name
	function ( x, groupcol='GroupName', colCol='GroupName', probeNames="Gene.Symbol" ) { ## Argumente der generischen Funktion
		standardGeneric('melt') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('melt', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, groupcol='GroupName', colCol='GroupName', probeNames="Gene.Symbol" ) {
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
})


#' plot grouped probesets creates ONE plot for ONE group of probesets
#' If you want a multi group plot create it yourself from the single ones.
#' @param x the ExpressionSet object
#' @param probeset the probeset name (rownames(x$data))
#' @param boxplot (F or T) create a a dots- or box-plot
#' @param pdf save the file as pdf (default svg)
#' @param geneNameCol the column name for the gene symbol to use in the plots
#' @export 
setGeneric('plot.probeset', ## Name
	function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol" ) { ## Argumente der generischen Funktion
		standardGeneric('plot.probeset') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.probeset', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol" ) {
	if ( sum(is.na(match(probeset, rownames(x$data)))==F) == 0 ){
		probeset <- rownames(x$data)[match( probeset, x$annotation[,geneNameCol] ) ]
	}
	if ( length(probeset) == 0 ) {time.col
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
})



#' drops samples from the ExpressionSet
#' @param x the ExpressionSet object
#' @param samplenames which samples to drop (samples like colnames(x$data))
#' @param name the name of the new ExpressionSet object
#' 
#' @export 
setGeneric('drop.samples', ## Name
	function ( x, samplenames=NULL, name='dopped_samples' ) { ## Argumente der generischen Funktion
		standardGeneric('drop.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('drop.samples', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, samplenames=NULL, name='dopped_samples' ) {
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
})

#' restrictSamples the ExpressionSet based on a variable in the samples table
#' @param x the ExpressionSet object
#' @param name the name of the new ExpressionSet
#' @param column which column to analyze in the samples table
#' @param value which value to take as 'cut off'
#' @param mode one of 'less', 'more', 'onlyless', 'equals'
#' @export 
setGeneric('restrictSamples', ## Name
	function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ) { ## Argumente der generischen Funktion
		standardGeneric('restrictSamples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('restrictSamples', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ) {
	
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
})

#' return the working directory using the linux pwd command
#' @export 
setGeneric('pwd', ## Name
	function () { ## Argumente der generischen Funktion
		standardGeneric('pwd') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('pwd', signature = c ('ExpressionSet') ,
	definition = function () {
	system( 'pwd > __pwd' )
	t <- read.delim( file = '__pwd', header=F)
	t <- as.vector(t[1,1])
	t <- paste(t,"/",sep='')
	unlink( '__pwd')
	t
})

#' force absolute unique names in a vector by adding _<amount of repeats> to each value if
#' there is more than one repeat opf the value
#' @param x the ExpressionSet object
#' @param separator '_' or anything you want to set the separator to
setGeneric('forceAbsoluteUniqueSample', ## Name
	function ( x ,separator='_') { ## Argumente der generischen Funktion
		standardGeneric('forceAbsoluteUniqueSample') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('forceAbsoluteUniqueSample', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x ,separator='_') {
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
})


#' add aditional annotation to the ExpressionSet
#' @param x the ExpressionSet object
#' @param mart the annotation table (data.frame or mart object)
#' @param mart.col which column corresponds to the rownames(x$data)
#' @export 
setGeneric('addAnnotation', ## Name
	function (x ,mart, mart.col='refseq_mrna') { ## Argumente der generischen Funktion
		standardGeneric('addAnnotation') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('addAnnotation', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (x ,mart, mart.col='refseq_mrna') {
	if ( ! class(mart) == 'data.frame' ){
		x$annotation <- cbind(x$annotation, mart[match(rownames(x$data),mart[,mart.col] ), ] )
	}
	else {
		x$annotation <- mart[is.na(match(mart[,mart.col], rownames(x$data)))==F,]
		rownames(x$annotation) <- rownames(x$data)
	}
	x
})

#' subsets the annotation table for a set of probesets
#' @param x the ExpressionSet object
#' @param probesets the vector of probesets to annotate (rownames(x$data))
#' @param colname a vector of colnames to annotate
#' @return a vector of annotation values
setGeneric('getAnnotation4probesets', ## Name
	function (x, probesets=c(), colname='Gene.Symbol' ) { ## Argumente der generischen Funktion
		standardGeneric('getAnnotation4probesets') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getAnnotation4probesets', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (x, probesets=c(), colname='Gene.Symbol' ) {
	as.vector(x$annotation[match( probesets, rownames(x$data) ), colname ])
})


#' creates a new ranks vector in the ExpressionSet
#' @param  x the ExpressionSet object
#' @export 
setGeneric('rank', ## Name
	function (x ) { ## Argumente der generischen Funktion
		standardGeneric('rank') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('rank', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (x ) {
	if ( ! exists ( 'ranks', where =x ) ){
		x$ranks <- apply( x$data,2,order)
		colnames( x$ranks ) <- colnames(x$data) 
		rownames( x$ranks ) <- rownames(x$data) 
	}
	x
})
#' reduces the dataset based on genes e.g. dropps genes from the ExpressionSet
#' @param x the ExpressionSet object
#' @param probeSets a list of probesets to reduce the data to
#' @param name the new ExpressionSet name
#' @export 
setGeneric('reduce.Obj', ## Name
	function ( x, probeSets=c(), name="reducedSet" ) { ## Argumente der generischen Funktion
		standardGeneric('reduce.Obj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reduce.Obj', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, probeSets=c(), name="reducedSet" ) {
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
})

#' plot one gene using either a points plot or a boxplot version of ggplot2
#' @param dat the ExpressionSet object
#' gene which gene to plot
#' @param colrs wich colors to use for the sample groups
#' @param groupCol which grouping variable to use for the grouping
#' @param colCol which grouping variable to use for the coloring (not used here)
#' @parma boxplot whther to plot a boxplot or points plot (default points)
#' @export
setGeneric('ggplot.gene', ## Name
	function (dat,gene, colrs, groupCol='GroupID', colCol='GroupID', boxplot=F) { ## Argumente der generischen Funktion
		standardGeneric('ggplot.gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('ggplot.gene', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (dat,gene, colrs, groupCol='GroupID', colCol='GroupID', boxplot=F) {
	not.in = NULL
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
})

#' draw a ggplot heatmap from a subset of genes in the Expression set
#' @param dat a ExpressionSet object
#' @param glist a list of genes to select using reduce.Obj
#' @param groupCol which grouping variable to use for the grouping
#' @param colCol which grouping variable to use for the coloring (not used here)
#' @export
setGeneric('gg.heatmap.list', ## Name
	function (dat,glist, colrs, groupCol='GroupID', colCol='GroupID') { ## Argumente der generischen Funktion
		standardGeneric('gg.heatmap.list') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('gg.heatmap.list', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (dat,glist, colrs, groupCol='GroupID', colCol='GroupID') {
	
	isect <- reduce.Obj ( dat, glist)
	#browser()
	dat.ss <- melt ( isect, probeNames=isect@rownamescol, groupcol=groupCol,colCol=colCol)
	#dat.ss <- dat[which(is.na(match(dat$Gene.Symbol,isect))==F),]
	colnames(dat.ss) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
	dat.ss$z <- ave(dat.ss$Expression, dat.ss$Gene.Symbol, FUN = function (x)  {
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
			not.in = setdiff( glist, rownames(isect@data)) )
})

#' draw a heatmap of a list of genes using the heatmap.2 function
#' @param dataOBJ the expressionSet object
#' @param gene.names the names of the genes to plot or nothing (all)
#' @param pvalue just an addition to the exported files (unusable)
#' @param analysis_name the name for this analysis (will be used for all files)
#' @param gene_centered collapse all data for one gene to one entry in the heatmap
#' @param Subset example: you want to look at J segments in a dataset, but not all separate but together in one 
#' row - use 'J_segment' for this option (not sure it works)
#' @param collaps merge all samples in one groups into one tile of the heatmap applying this grouing function (mean, median, max, sum)
#' @export 
setGeneric('plot.heatmaps', ## Name
	function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F,... ) { ## Argumente der generischen Funktion
		standardGeneric('plot.heatmaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.heatmaps', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F,... ) {
	dataOBJ <- normalize(dataOBJ)
	dataOBJ <- z.score( dataOBJ )
	
	if ( is.null(analysis_name)){
		analysis_name = dataOBJ@name
	}
	exprs = dataOBJ@data
	drop <- which( (apply(exprs,1,sd) == 0) == T)
	if ( length(drop) > 0 ){
		print ( paste("I need to drop",length(drop),"genes as the varianze is 0") )
		exprs = exprs[-drop,]
	}
	if ( ! is.null(gene.names)){
		geneNamesBased <- exprs[is.na(match(rownames(exprs),gene.names$all$id ))==F,]
		geneNamesBased <- data.frame ( 'GeneSymbol' = as.vector (dataOBJ@annotation[is.na(match ( dataOBJ@annotation$GeneID, gene.names$all$id ) ) == F, geneNameCol]),  geneNamesBased) 
	}
	else {
		geneNamesBased <- data.frame( 'GeneSymbol' = dataOBJ@annotation[rownames(exprs),geneNameCol], exprs )
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
	fname = paste(dataOBJ@outpath, 'GenesOnly_',dataOBJ@name,sep='')
	rn <- rownames(geneNamesBased)
	
	if ( ! is.null(collaps)){
		u <- unique(as.vector(dataOBJ@samples$GroupName))
		m <- length(u)
		mm <-  matrix ( rep(0,m * nrow(geneNamesBased)), ncol=m)
		colnames(mm) <- u
		rownames(mm) <- rownames(geneNamesBased)
		if ( collaps=="median") {
			for ( i in u ){
				## calculate the mean expression for each gene over all cells of the group
				mm[,i] <- apply( geneNamesBased[ , which(as.vector(dataOBJ@samples$GroupName) == i )+1],1,median)
			}
		}
		else {
			for ( i in u ){
				## calculate the mean expression for each gene over all cells of the group
				mm[,i] <- apply( geneNamesBased[ , which(as.vector(dataOBJ@samples$GroupName) == i )+1],1,mean)
			}
		}
		geneNamesBased <- data.frame( GeneSymbol = as.vector(geneNamesBased[,1]), mm )
	}else {
		collaps= ''
	}
	
	write.table( data.frame ( 'GeneSymbol' = rownames(geneNamesBased),geneNamesBased[,-1]),file= paste(fname,collaps,'_HeatmapID_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table(  geneNamesBased, file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''),sep='\t', row.names=F,quote=F  )
	write.table( dataOBJ@samples, file= paste(fname,collaps,'_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	
	if ( nrow(geneNamesBased) > 1 ){
		geneNamesBased <- data.frame(read.delim( file= paste(fname,collaps,'_Heatmap_GeneSymbols_',pvalue,'_data4Genesis.txt', sep=''), header=T ))
		rownames(geneNamesBased) <- rn
		s <- ceiling(20 * nrow(geneNamesBased)/230 )
		if ( s < 5 ) {s <- 5}
		print (paste ("Height =",s))
		if ( pdf ) {
			pdf( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_IDs_",pvalue,".pdf",sep='') ,width=10,height=s)
		}else {
			devSVG( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_IDs_",pvalue,".svg",sep='') ,width=10,height=s)
		}
		print ( paste ("I create the figure file ",dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_IDs_",pvalue,".svg",sep=''))
		heatmap.2(as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", margin=c(10, 10),dendrogram="row",scale="row",
'distfun'= function (x) as.dist( 1- cor(t(x), method='pearson')),Colv=F, main=paste( analysis_name, pvalue ),...)
		dev.off()
		
		rownames(geneNamesBased) <- paste(geneNamesBased[,1] , rownames(geneNamesBased) )
		if ( pdf ) {
			pdf( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_GeneSymbols_",pvalue,".pdf",sep='') ,width=10,height=s)
		}else{
			devSVG( file=paste(dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_GeneSymbols_",pvalue,".svg",sep='') ,width=10,height=s)
		}
		
		print ( paste ("I create the figure file ",dataOBJ@outpath,dataOBJ@name,collaps,"_Heatmap_GeneSymbols_",pvalue,".svg",sep=''))
		heatmap.2( as.matrix(geneNamesBased[,-1]), lwid = c( 1,6), lhei=c(1,5), cexRow=0.4,cexCol=0.7,col = greenred(30), trace="none", 
distfun = function (x) {as.dist( 1- cor(t(x), method='pearson') ) } ,...)
		dev.off()
	}
	else {
		print ( "Problems - P value cutoff to stringent - no selected genes!" );
	}
})

#' groups.boxplot is old code, that I am likely removing soon
setGeneric('groups.boxplot', ## Name
	function ( x, SampleCol='GroupName', clusters, svg=F, fname='group_', width=800, height=800,mar=NULL, Collapse=NULL, ...) { ## Argumente der generischen Funktion
		standardGeneric('groups.boxplot') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('groups.boxplot', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, SampleCol='GroupName', clusters, svg=F, fname='group_', width=800, height=800,mar=NULL, Collapse=NULL, ...) {
	maxG <- max( clusters )
	ret <- list()
	r=1
	gnames <- unique(as.vector(x@samples[,SampleCol]))
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
			fnames[i] = paste(x@outpath,fname,i,"_boxplot.C.svg",sep='')
			devSVG ( file= fnames[i], width=width/130, height=height/130 )
		}else{
			fnames[i] =paste(x@outpath,fname,i,"_boxplot.png",sep='')
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
	#print (paste('montage', paste(fnames, collapse= " "), "-geometry", paste("'", width, "x", height,"+0+0'",sep=''),paste(x@outpath,fname,"montage.png",sep=''), sep=' ' ))
	try( file.remove(  paste(x@outpath,fname,"montage.png",sep='') ) , silent=T )
	system ( paste('montage', fnames, "-geometry", paste("'", width, "x", height,"+0+0'",sep=''),paste(x@outpath,fname,"montage.png",sep='')," 2>/dev/null", collapse=' ' ) )
	ret
})

#' obtained from https://stackoverflow.com/questions/8261590/write-list-to-a-text-file-preserving-names-r
#' writed a list to file preserving the names of the list.
#' @param x a list you want to save (human readable)
#' @param fname the outfile
#' @param sep default ' '
#' @export 
setGeneric('write.list', ## Name
	function ( x, fname, sep=' ') { ## Argumente der generischen Funktion
		standardGeneric('write.list') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('write.list', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, fname, sep=' ') {
	z <- deparse(substitute(x))
	cat(z, "\n", file=fname)
	nams=names(x) 
	for (i in seq_along(x) ){ cat(nams[i], sep,  x[[i]], "\n", 
				file=fname, append=TRUE) 
	}
})

#' Calculates a simple anova model on the data combining all informatrion into one single model.
#' @param x the ExpressionSet object
#' @param samples.col the column in the samples table that contains the grouping string (e.g. GroupName)
#' @param padjMethod anything the \code{\link[stats]{p.adjust}} method does support ('BH')
#' @export 
setGeneric('simpleAnova', ## Name
	function ( x, samples.col='GroupName', padjMethod='BH' ) { ## Argumente der generischen Funktion
		standardGeneric('simpleAnova') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('simpleAnova', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, samples.col='GroupName', padjMethod='BH' ) {
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
})

#' constructors that have to be implemented in the classes
#' the createStats constructor that has to be implemented in the data specific packages
setGeneric('createStats', ## Name
	function ( x, condition, files=F, A=NULL, B=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('createStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('createStats', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, condition, files=F, A=NULL, B=NULL ) {
	stop( 'Not implemented' )
})

#' constructors that have to be implemented in the classes
#' the normalize constructor that has to be implemented in the data specific packages
setGeneric('normalize', ## Name
	function ( x ) { ## Argumente der generischen Funktion
		standardGeneric('normalize') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('normalize', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x ) {
	x
})

#' probably uselecc function that makes sure all data values are numeric
#' @param dataObj the ExpressionSet object
setGeneric('force.numeric', ## Name
	function (dataObj ) { ## Argumente der generischen Funktion
		standardGeneric('force.numeric') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('force.numeric', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( dataObj ) {
	for ( i in 1: ncol(dataObj@data) ) { 
		if ( !  paste(as.vector(dataObj@data[,i]), collapse=" ") == paste(as.vector(as.numeric(dataObj@data[,i])), collapse=" ") ) { 
			dataObj@data[,i] <- as.vector(as.numeric(dataObj@data[,i]))
		}
	}
	dataObj
})

#' print the ExpressionSet 
#' @param x the ExpressionSet object
#' @return nothing
#' @example 
#' 
#' @export 
setGeneric('show', ## Name
	function (x) { ## Argumente der generischen Funktion
		standardGeneric('show') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('show', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function (x) {
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
})



#' write the ExpressionSet data table to disk
#' @param x the ExpressionSet object
#' @param annotation a vector of annotation data column names to include in the written table (default=none)
#' @export 
setGeneric('write.data', ## Name
	function ( x, annotation=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('write.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('write.data', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, annotation=NULL ) {
	if ( !is.null(annotation) ) {
		write.table( cbind( x$annotation[,annotation], x$data), file= paste( x$outpath,x$name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
	else {
		write.table( cbind( rownames(x$data), x$data), file= paste( x$outpath,x$name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
})

#' returns a list of probesets (the rownames from the data matrix) for a restriction of a list of stat comparisons
#' 
#' @param v The cutoff value
#' @param pos The column in the stats tables to base the selection on
#' @param Comparisons A list of comparisons to check (all if left out)
#' @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
#' @examples 
#' probes <- getProbesetsFromStats ( x, v=1e-4, pos="adj.P.Val" ) 
#' @return a list of probesets that shows an adjusted p value below v (1e-4 in the example)
#' @export 
setGeneric('getProbesetsFromStats', ## Name
		function ( x, v=0.05, pos=6, mode='less', Comparisons=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('getProbesetsFromStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('getProbesetsFromStats', signature = c ('ExpressionSet','ExpressionSet') ,
definition = function ( x, v=0.05, pos=6, mode='less', Comparisons=NULL ) {
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
})

#' export the statistic files
#' @param x the NGSexpressionSet
#' @export 
setGeneric('writeStatTables', ## Name
	function ( x, annotation=NULL, wData=F ) { ## Argumente der generischen Funktion
		standardGeneric('writeStatTables') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('writeStatTables', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, annotation=NULL, wData=F ) {
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
})

#' exports the ExpressionSet in the format \href{https://apps.childrenshospital.org/clinical/research/ingber/GEDI/gedihome.htm}[GEDI] can import.
#' @param x the ExpressionSet object
#' @param fname the filename to export the data to
#' @param tag.col the sample name column in the samples table
#' @param time.col the time column in the samples table
#' @param minSample_PerTime drop a timepoint if not at least (1) sample is in the group
#' @export 
setGeneric('export4GEDI', ## Name
	function ( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) { ## Argumente der generischen Funktion
		standardGeneric('export4GEDI') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('export4GEDI', signature = c ('ExpressionSet','ExpressionSet') ,
	definition = function ( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) {
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
})
