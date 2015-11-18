

read.bams <- function ( bamFiles, annotation, GTF.featureType='exon', GTF.attrType = "gene_id", isPairedEnd = FALSE, nthreads = 2){
	if (file.exists(bamFiles)){
		bamFiles <- readLines(bamFiles)
	}
	counts <- featureCounts(files =bamFiles,annot.ext = annotation ,isGTFAnnotationFile = TRUE,GTF.featureType = GTF.featureType,
		GTF.attrType = GTF.attrType,allowMultiOverlap=T, isPairedEnd =isPairedEnd , nthreads = nthreads)
	counts.tab <- cbind(counts$annotation,counts$counts)  # combine the annotation and counts
	counts.tab
}


NGSexpressionSet <- function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) {
	x <- createWorkingSet ( dat, Samples,  Analysis = Analysis, name=name, namecol=namecol, namerow= namerow, usecol= usecol, outpath =  outpath)
	class( x ) <- append(  'NGSexpressionSet', class(x))
	x
}

getProbesetsFromStats <- function ( x, v=0.05, pos=6, mode='less', Comparisons=NULL ) {
	UseMethod('getProbesetsFromStats', x)
}

# normalizeReads.NGSexpressionSet <- function ( x, readCounts  )
# @param x The NGSexpressionSet
# @param readCounts The number of reads from each bam file

normalize.NGSexpressionSet <- function (  x, readCounts=NULL, to_gene_length=FALSE, geneLengthCol='transcriptLength' ) {
	if ( is.null( readCounts ) ) {
		readCounts <- as.vector( estimateSizeFactorsForMatrix ( x$data) )
	}
	x$samples$SizeFactor <- readCounts
	if (! exists('normalized', where=x ) ) {
		x$normalized = 0
	}
	if ( x$normalized == 0 ) {
		x$raw <- x$data
		x$data =  t(apply(x$data,1, function(a) { a / readCounts } ) )
		x$normalized = 1
		
		if (to_gene_length){
			for ( i in 1:nrow(x$data)) {
				x$data[i,] <- x$data[i,]/ (x$annotation[i,geneLengthCol ] / 1000)
			}
		}
	}
	x
}

# getProbesetsFromStats returns a list of probesets (the rownames from the data matrix) for a restriction of a list of stat comparisons
# @param v The cutoff value
# @param pos The column in the stats tables to base the selection on
# @param Comparisons A list of comparisons to check (all if left out)
# @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
# @examples 
# probes <- getProbesetsFromStats ( x, v=1e-4, pos="adj.P.Val" ) 
# returns a list of probesets that shows an adjusted p value below 1e-4

getProbesetsFromStats.NGSexpressionSet <- function ( x, v=0.05, pos=6, mode='less', Comparisons=NULL ) {
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

getProbesetsFromValues.NGSexpressionSet <- function ( x, v='NULL', sample='NULL', mode='less' ){
	s <- FALSE
	if ( is.null(v) ){
		s<-TRUE
	}
	if ( is.null(sample) ){
		s<-TRUE
	}
	if ( s ) { stop ( "Please give me the required values for v and sample") }
	probesets <- NULL
	for ( s in sample ) {
	switch( mode,
			'less' = probesets <- c(probesets, as.vector(rownames(x$data)[which(x$data[,s] <= v)] ) ) ,
			'more' = probesets <- c(probesets, as.vector(rownames(x$data)[which(x$data[,s] > v)] ) ), 
			'onlyless' = probesets <- c(probesets,  as.vector(rownames(x$data)[which(x$data[,s] < v)] ) ),
			'equals' = probesets <- c(probesets, as.vector(rownames(x$data)[which(x$data[,s] == v)] ) )
	)
	}
	unique(probesets)
}


name_4_IDs.NGSexpressionSet <- function ( x, ids=NULL, geneNameCol='mgi_symbol' ) {
	if ( is.null(ids) ) {
		ids <- as.vector(colnames(x$data) )
	}
	as.vector(x$annotation[match( ids,x$annotation[, x$rownamescol]),geneNameCol])
}

get_gene_list.NGSexpressionSet <- function (x, p_value = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), geneNameCol='mgi_symbol' ) {
	if ( exists(where=x, 'stats')){
		cmps <- names(x$stats)
		p_value <- as.numeric( p_value )
		x$sig_genes <- vector ('list', length(p_value))
		names( x$sig_genes ) <- as.vector(p_value)
		for ( p in 1:length(p_value)){
			x$sig_genes[[p]] <- vector ('list', length(cmps)-1)
			for ( i in 2:length(cmps) ){
				sig.genes <- x$stats[[i]][which(x$stats[[i]]$padj < p_value[p] ), ] 
				sig.names <- name_4_IDs.NGSexpressionSet( x, sig.genes[,1], geneNameCol)
				sig.tab <- cbind(sig.names,sig.genes ) 
				if ( ncol(sig.tab) > 0 ){
					write.table(sig.tab,paste(x$outpath,x$name,'_',cmps[i],p_value[p] ,".xls",sep=''),col.names=T,row.names=F,sep="\t",quote=F)
				}
				x$sig_genes[[p]][[i-1]] = list (id = sig.genes[,1], names=sig.names )
			}
			x$sig_genes[[p]]$all <- list( id = unique(unlist(lapply (x$sig_genes[[p]] , function(a) { a$id } ))), names= unique(unlist(lapply ( x$sig_genes[[p]], function(a) { a$names } ))) )
		}
	}
	x
}

plot.NGSexpressionSet <- function ( x, pvalue=c( 0.1,1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), Subset=NULL , Subset.name= NULL, comp=NULL, gene_centered=F, collaps=NULL,geneNameCol= "mgi_symbol") {
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


collaps <- function ( dataObj, what=c('median','mean','sd','sum' ) ) {
	UseMethod('collaps', dataObj)
}
collaps.NGSexpressionSet<- function(dataObj, what=c('median','mean','sd','sum' ) ) {
	u <- unique(as.vector(dataObj$samples$GroupName))
	m <- length(u)
	mm <-  matrix ( rep(0,m * nrow(dataObj$data)), ncol=m)
	colnames(mm) <- u
	rownames(mm) <- rownames(dataObj$data)
	f <- NULL
	switch( what,
			median = f<- function (x ) { median(x) },
			mean = f <- function(x) { mean(x) },
			sd = f <- function(x) { sd(x) },
			sum = f <-function(x) { sum(x)}
	);
	if ( is.null(f) ) {
		stop("Please set what to one of 'median','mean','sd','sum'" )
	}
	new_samples <- NULL
	for ( i in u ){
		all <- which(as.vector(dataObj$samples$GroupName) == i )
		new_samples <- rbind ( new_samples, dataObj$samples[all[1],] )
		mm[,i] <- apply( dataObj$data[ , all],1,f)
	}
	dataObj$name = paste( dataObj$name, what, sep='_')
	dataObj$data <- mm
	dataObj$samples <- new_samples
	dataObj
}

removeBatch.NGSexpressionSet <- function( x, phenotype ){
	if ( x$batchRemoved==1) {
		return (x)
	}
	browser()
	exprs =  dataOBJ$cds
	null <- which ( exprs == 0)
	exprs[null]<- 1
	log <- log(exprs)
	svseq <-  ComBat(dat=filtered.log , mod=mod1, batch=x$samples[,phenotype], par.prior=TRUE, prior.plots=FALSE )
	svseq <- exp(log)
	svseq[null] = 0
	x$cds <- svseq
	x$batchRemoved = 1
	x
}

check <- function ( x , genes=NULL, cutoff=0.77 ) {
	UseMethod('check', x)
}

check.NGSexpressionSet <- function (x, genes=NULL, cutoff = 0.77 ) {
	percent5 <-  reads.taken(x, 0.05, genes)
	names(which(percent5$reads.taken > cutoff )) ## obtained experimentally using Jennies HSC dataset
}

reads.taken <- function ( x , percentile= 0.05, tf=NULL) {
	UseMethod('reads.taken', x)
}
reads.taken.NGSexpressionSet <- function ( x, percentile= 0.05, tf=NULL ) {
	top.genes <- list()
	reads.taken <- vector( 'numeric', ncol(x$data))
	nTF <- vector('numeric',  ncol(x$data)+1)
	percentile= 1- percentile
	for ( i in 1:ncol(x$data) ){
		qua <- quantile(x$data[,i], percentile)
		reads.taken [i] <- sum (x$data[which(x$data[,i] > qua),i] ) / sum(x$data[,i])
		top.genes[[i]] <- rownames(x$data) [ which(x$data[,i] > qua) ]
		if ( ! is.null(tf) ){
		nTF[i] <- length( intersect( tf, str_replace(top.genes[[i]],  "\\.[0-9]*" ,'') ) )
		}
	}
	names( reads.taken) <- colnames(x$data)
	names(top.genes) <- colnames(x$data)
		
	inter <- intersect( top.genes[[1]], top.genes[[2]])
	for (i in 3:ncol(x$data) ) { inter <- intersect( inter, top.genes[[i]]) }
	if ( ! is.null(tf) ) {nTF[ncol(x$data)+1] <- length(intersect(str_replace(inter,  "\\.[0-9]*" ,''), tf ) )}
	reads.taken.intersect <- vector( 'numeric', ncol(x$data))
	for ( i in 1:ncol(x$data) ){
		reads.taken.intersect [i] <- sum ( x$data[inter ,i] ) / sum(x$data[,i])
	}
	names( reads.taken.intersect) <- colnames(x$data)
	list( reads.taken = reads.taken, top.genes = top.genes, intersect = inter,reads.taken.intersect = reads.taken.intersect, nTF = nTF )
}
	
	
preprocess.NGSexpressionSet <- function (x) {
	if ( ! exists(where=x, 'vsdFull')){
		condition <- as.factor(x$samples$GroupName)
		print ( condition )
		x$cds <- newCountDataSet(x$data, condition)
		x$cds <- estimateSizeFactors(x$cds) 
		sizeFactors(x$cds)
		x$cds <- estimateDispersions(x$cds) 
		x$vsdFull = varianceStabilizingTransformation( x$cds )
	}
	x
}

z.score <- function ( m ){
	UseMethod('z.score', m)
}
z.score.matrix <- function (m ) {
	rn <- rownames( m )
	me <- apply( m, 1, mean )
	sd <- apply( m, 1, sd )
	sd[which(sd==0)] <- 1e-8
	m <- (m - me) /sd
	rownames(m) <- rn
	m
}
z.score.NGSexpressionSet <- function(m) {
	#m$data <- z.score( as.matrix( m$data ))
	rn <- rownames( m$data )
	me <- apply( m$data, 1, mean )
	sd <- apply( m$data, 1, sd )
	sd[which(sd==0)] <- 1e-8
	m$data <- (m$data - me) /sd
	rownames(m$data) <- rn
	m
}


anayse_working_set <- function ( dataOBJ, name,  p_values = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ),geneNameCol= "mgi_symbol", batchCorrect=NULL)  {
	dataOBJ <- preprocess.NGSexpressionSet ( dataOBJ )
	if ( ! is.null(batchCorrect) ){
		removeBatch.NGSexpressionSet( dataOBJ ,batchCorrect )
	}
	png(file=paste(dataOBJ$name,"_",version,'_interesting_samples_clean_up.png',sep=''), width=800, height=800) 
	plot(hclust(as.dist(1-cor(dataOBJ$data))),main=paste(dataOBJ$name,version, ' samples')) 
	dev.off()
	dataOBJ <- do_comparisons.NGSexpressionSet ( dataOBJ, geneNameCol=geneNameCol)
	#dataOBJ$expr <- exprs(dataOBJ$vsdFull)
	dataOBJ <- get_gene_list.NGSexpressionSet(dataOBJ,p_values, geneNameCol=geneNameCol)	
	dataOBJ
}
simpleAnova <- function(x, samples.col='GroupName', padjMethod='BH' ) {
	UseMethod ('simpleAnova', x)
}

simpleAnova.NGSexpressionSet <- function ( x, samples.col='GroupName', padjMethod='BH' ) {
	x <- normalize(x)
	Samples <- x$samples
	significants <- apply ( x$data ,1, function(x) { anova( lm (x ~ Samples[, samples.col]))$"Pr(>F)"[1] } )
	adj.p <- p.adjust( significants, method ='BH' )
	res <- cbind(significants,adj.p )
	res <- data.frame(cbind( rownames(res), res ))
	colnames(res) <- c('genes', 'pvalue', paste('padj',padjMethod) )
	res[,2] <- as.numeric(as.vector(res[,2]))
	res[,3] <- as.numeric(as.vector(res[,3]))
	res <- list ( 'simpleAnova' = res )
	if ( exists( 'stats', where=x )) {
		x$stats <- c( x$stats, res)
	}else {
		x$stats <- res
	}
	x
}



createStats.NGSexpressionSet <- function (x, condition, files=F, A=NULL, B=NULL) {
	if ( nrow(x$data) < 2e+4 ) {
	    stop ( "Please calculate the statistics only for the whole dataset!" )
	}
	if ( length( grep ( condition, colnames(x$samples))) > 0 ) {
		condition = factor( x$samples[,condition] )
	}
	if ( exists( 'raw', where=x) ){
		cds <- newCountDataSet(x$raw, condition)
	}else {
		cds <- newCountDataSet(x$data, condition)
	}
	cds <- estimateSizeFactors(cds)
	# sizeFactors(cds)
	cds <- estimateDispersions(cds)
	vsdFull = varianceStabilizingTransformation( cds )
	res <- list()
	ln= 0
	na <- c()
	conditions <- as.vector(unique(condition))
	if ( ! is.null(A) && ! is.null(B)) {
		ln = ln +1
		na[ln] = paste( A, B ,sep=' vs. ')
		res[[ln]] <- nbinomTest(cds, A, B )
	}
	else {
		for ( i in 1:(length(conditions)-1) ){
			for ( a in (i+1):length(conditions) ){
				ln = ln +1
				print ( paste( conditions[i], conditions[a], sep='_vs_') )
				na[ln] = paste( conditions[i], conditions[a],sep=' vs. ')
				res[[ln]] <- nbinomTest(cds, conditions[i], conditions[a])
				#res[[ln]] <- res[[res[[1]]+1]][which(res[[res[[1]]+1]]$padj < 0.2),]
			}
		}
	}
	names(res) <- na
	if ( exists( 'stats', where=x )) {
		x$stats <- c( x$stats, res)
	}else {
		x$stats <- res
	}
	if ( files ) {
		writeStatTables( x )
	}
	x
}

describe.probeset <-  function ( x, probeset ) {
	        UseMethod('describe.probeset', x)
}
describe.probeset.NGSexpressionSet <-  function ( x, probeset ) {
	ret <- list()
	print ( paste ("Annoataion for probeset ",probeset ) )
	ret$annotation <- x$annotation[ probeset, ] 
	print ( ret$annotation )
	print ( "Values:" ) 
	ret$data <- x$data[probeset, ]
	print ( ret$data )
	ret$stats <- NULL
	#browser()
	if ( exists( 'stats', where=x) ) {
		for( i in 1:length(x$stats) ) {
	#		browser()
			ret$stats <- rbind(ret$stats, 
					   cbind ( 
						  rep(names(x$stats)[i],length(probeset) ), 
						  x$stats[[i]][is.na(match( x$stats[[i]][,1], probeset))==F,] 
						  ) 
					   ) 
		}
		colnames(ret$stats) <- c('comparison', colnames(x$stats[[1]]) )
		print ( ret$stats )
	}
	invisible(ret)
}


# getProbesetsFromStats returns a list of probesets (the rownames from the data matrix) for a restriction of a list of stat comparisons
# @param v The cutoff value
# @param pos The column in the stats tables to base the selection on
# @param Comparisons A list of comparisons to check (all if left out)
# @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
# @examples 
# probes <- getProbesetsFromStats ( x, v=1e-4, pos="adj.P.Val" ) 
# returns a list of probesets that shows an adjusted p value below 1e-4

getProbesetsFromStats.NGSexpressionSet <- function ( x, v=1e-4, pos='padj', mode='less', Comparisons=NULL ) {
	if ( is.null(Comparisons)){	Comparisons <- names(x$stats) }
	probesets <- NULL
	for ( i in match(Comparisons, names(x$stats) ) ) {
		switch( mode,
				'less' = probesets <- c( probesets, as.vector(x$stats[[i]][which(x$stats[[i]][,pos] <= v),1] )),
				'more' = probesets <- c( probesets, as.vector(x$stats[[i]][which(x$stats[[i]][,pos] > v),1] )), 
				'onlyless' = probesets <- c( probesets, as.vector(x$stats[[i]][which(x$stats[[i]][,pos] < v),1] )),
				'equals' = probesets <- c( probesets, as.vector(x$stats[[i]][which(x$stats[[i]][,pos] == v),1] ))
		)
	}
	unique(probesets)
}

getProbesetsFromValues <- function ( x, v=NULL, sample=NULL, mode='less' ){
	UseMethod('getProbesetsFromValues', x)
}

# Select probesets, that show a certain level in expression for a single sample
# @param v The cutoff value
# @param sample The sample name
# @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
# @examples 
# probes <- getProbesetsFromStats ( x, v=10, sample="A" ) 
# returns a list of probesets that has a expression less than 10 in sample A

getProbesetsFromValues.NGSexpressionSet <- function ( x, v=NULL, sample=NULL, mode='less' ){
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
