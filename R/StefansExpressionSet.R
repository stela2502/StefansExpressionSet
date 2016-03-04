##perl -e ' foreach (@ARGV) { print "library($_)\nrequire($_)\n" }' 'reshape2' 'gplots' 'stringr' 'RSvgDevice' 'rgl' 'ggplot2' 
library(reshape2)
#require(reshape2)
require(gplots)
require(stringr)
require(RSvgDevice)
require(rgl)
require(ggplot2)

#' @name StefansExpressionSet
#' @title StefansExpressionSet
#' @description  An S4 class to visualize Expression data.
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot zscored genes are normalized?
#' @slot snorm samples normalized?
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot stats the stats list for all stats created in the object
#' @slot usedObj here a set of used and probably lateron important objects can be saved. Be very carful using any of them!
#' @exportClass StefansExpressionSet
setClass(
		Class='StefansExpressionSet', 
		representation = representation ( 
			data='data.frame',
			samples='data.frame',
			ranks='numeric',
			raw="data.frame",
			annotation='data.frame',
			snorm='logical',
			outpath='character',
			name='character',
			rownamescol='character',
			sampleNamesCol='character',
			stats = 'list',
			sig_genes = 'list',
			zscored = 'logical',
			simple = 'character',
			usedObj = 'list'
		),
		prototype(outpath ='', name = 'StefansExpressionSet',
				sampleNamesCol=NA_character_, 
				stats=list(),
				snorm=F,
				zscored=F,
				sig_genes=list(),
				usedObj= list(),
				simple= c( 'outpath', 'rownamescol', 'sampleNamesCol', 'simple', 'snorm', 'zscored') )
)


#' @name PMID25158935exp
#' @title Read counts for the expression data described in PMID25158935
#' @description The data was re-mapped against mouse mm10 using HISAT
#' @description and quantified using the R subreads package.
#' @docType data
#' @usage PMID25158935exp
#' @format data.frame
'PMID25158935exp'

#' @name PMID25158935samples
#' @title Read counts for the sample data for the expression information described in PMID25158935
#' @description The data was collected from the NCBI SRA archive
#' @docType data
#' @usage PMID25158935samples
#' @format data.frame
'PMID25158935samples'

#' @name red
#' @title reduced PMID25158935 dataset to a 100x15 StefansExpressionSet
#' @description Reduced StefansExpressionSet from the PMID25158935exp + PMID25158935samples dataset
#' @docType data
#' @usage red
#' @format StefansExpressionSet
'red'

#' @name StefansExpressionSet
#' @aliases StefansExpressionSet,data.frame-method
#' @rdname StefansExpressionSet-methods
#' @docType methods
#' @description  this file contains all generic fnction for data export and ploting Create an StefansExpressionSet
#' @description  object (S3) This object is mainly used for subsetting of the data and plotting @export
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
#' @param outpath Where to store the output from the analysis
#' @param annotation The annotation table from e.g. affymetrix csv data
#' @param newOrder The samples column name for the new order (default 'Order')
#' @title description of function StefansExpressionSet
#' @exportMethod  StefansExpressionSet
setGeneric("StefansExpressionSet", ## Name
		function( dat, Samples, class='StefansExpressionSet',  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = ''){ ## Argumente der generischen Funktion
			standardGeneric("StefansExpressionSet") ## der Aufruf von standardGeneric sorgt für das Dispatching
		})

setMethod("StefansExpressionSet", signature = c ('data.frame'), 
	definition = function ( dat, Samples, class='StefansExpressionSet',  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = '' ) {
	S <- Samples
	if ( ! is.null(Analysis) ){
		S <- Samples[which ( Samples$Analysis == Analysis ), ]
	}
	if ( ! is.null(usecol) ) {
		S <- Samples[which ( Samples[, usecol] == 1 ),]
	}
	if ( exists('filename',S) ) {
		n <- make.names(as.vector(S$filename))
		ret <- dat[, n ]
		annotation <- dat[, is.na(match( colnames(dat), n ))==T ]
	}else{
		n <- make.names(as.vector(S[,namecol]))
		ret <- dat[, n ]
		annotation <- dat[, is.na(match( colnames(dat), n ))==T ]
	}
	
	if ( exists( 'Order', S)){
		ret <- ret[, order(S$Order)  ]
		S <- S[order(S$Order), ]
	}
	
	if ( outpath == '' ){
		outpath = pwd()
	}
	if ( ! file.exists(outpath)){
		dir.create( outpath )
	}
	data <- list ( 'data' = data.frame(ret), samples = S, name= name, annotation = annotation, rownamescol = namerow )
	data$outpath <- outpath
	
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
	data$sampleNamesCol <- namecol
	class(data) <- 'StefansExpressionSet'
	data$batchRemoved=0
	
	r <- new ( class, data = data$data, samples = data$samples, name = data$name, annotation = data$annotation, rownamescol= data$rownamescol,sampleNamesCol = data$sampleNamesCol , outpath= data$outpath )
	force.numeric(r)
})


#' @name coexprees4gene
#' @aliases coexprees4gene,StefansExpressionSet-method
#' @rdname coexprees4gene-methods
#' @docType methods
#' @description  Calculate the coexpression for any given gene
#' @param x the StefansExpressionSet varibale
#' @param method any method supported by \link[=cor.test]{cor.test}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param padjMethod the method to calucate the FDR with \link[=p.adjust]{p.adjust}
#' @return a data.frame with the columns 'GeneSymbol', 'pval', 'cor', 'adj.p'
#' @title description of function coexprees4gene
#' @export 
setGeneric('coexprees4gene', ## Name
	function ( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) { ## Argumente der generischen Funktion
		standardGeneric('coexprees4gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('coexprees4gene', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, gene=NULL, method='spearman', geneNameCol='gene_name', padjMethod='BH' ) {
	ret <- NULL
	if ( ! is.null(gene) ){
		z <- as.vector( t(x@data[ gene[1] ,]) )
		pval <- vector( 'numeric', nrow(x@data))
		cor <- vector( 'numeric', nrow(x@data))
		for ( i in 1:nrow(x@data) ) {
			try( {	res <-  cor.test( z, as.vector(t(x@data[i,]), 'numeric') ,method=method)
			pval[i] <- res$p.value
			cor[i] <- res$estimate }, silent=T
			)
		}
		adj.p <- p.adjust(pval , method = padjMethod) #BH
		ret <- data.frame('GeneSymbol' = x@annotation[,geneNameCol], pval = pval, cor = cor, adj.p = adj.p )
		rownames(ret) <- rownames(x@data)
	}
	ret
})
 
#' @name pcorMat
#' @aliases pcorMat,StefansExpressionSet-method
#' @rdname pcorMat-methods
#' @docType methods
#' @description  calculate a p_value matrix for the data object THIS FUNCTION IS NOT DOING THE RIGTH
#' @description  THING BROKEN
#' @param x the StefansExpressionSet varibale
#' @param sd_cut the cut off value for the sd check
#' @param method any method supported by \code{\link[stats]{cor.test}}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param groupCol the column name of the grouping variable in the samples table
#' @param name the name of the analysis
#' @title description of function corMat.Pvalues
#' @export 
setGeneric('pcorMat', ## Name
	function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) { ## Argumente der generischen Funktion
		standardGeneric('pcorMat') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('pcorMat', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	# TODO: implement the p value calculation!
	d <- reduce.Obj( x, rownames(x@data)[which( apply(x@data,1,sd) > sd_cut)], name =name )
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



#' @name corMat
#' @aliases corMat,StefansExpressionSet-method
#' @rdname corMat-methods
#' @docType methods
#' @description  calculate a correlation matrix for the data object
#' @param x the StefansExpressionSet varibale
#' @param sd_cut the cut off value for the sd check
#' @param method any method supported by \code{\link[stats]{cor.test}}
#' @param geneNameCol the name of the gene column (gene_name)
#' @param groupCol the column name of the grouping variable in the samples table
#' @param name the name of the analysis
#' @return the correlation matrix
#' @title description of function corMat
#' @export 
setGeneric('corMat', ## Name
	function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) { ## Argumente der generischen Funktion
		standardGeneric('corMat') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('corMat', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, sd_cut=1, method='spearman', geneNameCol='gene_name', groupCol=NULL, name='tmp' ) {
	d <- reduce.Obj( x, rownames(x@data)[which( apply(x@data,1,sd) > sd_cut)], name = name )
	if ( ! is.null(groupCol) ){
		ret <- list()
		names <- unique(d$samples[,groupCol])
		for ( i in 1:length(names)) {
			a <- subset( d, column=groupCol, value=names[i], name= paste(d$name,names[i],sep='_'), mode='equals' )
			ret[[i]] = corMat( a, sd_cut= sd_cut,method=method, geneNameCol=geneNameCol )
		}
		names(ret) <- namestime.col
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
#' @name cor2cytoscape
#' @aliases cor2cytoscape,StefansExpressionSet-method
#' @rdname cor2cytoscape-methods
#' @docType methods
#' @description  exports the correlation matrix in the format \url{http://www.cytoscape.org/} does
#' @description  support as network file.
#' @param M the correlation matrix obtained by a run of \code{\link{corMat}}
#' @title description of function cor2cytoscape
#' @export 
setGeneric('cor2cytoscape', ## Name
	function (M, file, cut=0.9 ) { ## Argumente der generischen Funktion
		standardGeneric('cor2cytoscape') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

###last

setMethod('cor2cytoscape', signature = c ( 'StefansExpressionSet') ,
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
#
#' @name melt
#' @aliases melt,StefansExpressionSet-method
#' @rdname melt-methods
#' @docType methods
#' @description  met an StefansExpressionSet to be plotted using ggplot2 functions
#' @param data the StefansExpressionSet object
#' @param groupcol the column in the samples table to group the expression on
#' @param colCol samples column names(s) to color the grouping data
#' @param rowCol annotation column name(s) describing the gene groups
#' @param probeNames which probenames to use (ProbeSetID or Gene.Symbol ...)
#' @title description of function melt
#' @export 
setGeneric('melt.StefansExpressionSet', ## Name
		package = 'StefansExpressionSet',
	function ( dat, groupcol='GroupName', colCol=NULL, rowCol=NULL, probeNames=NULL,  na.rm = FALSE, value.name = "value") { ## Argumente der generischen Funktion
		standardGeneric('melt.StefansExpressionSet') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)


setMethod('melt.StefansExpressionSet',
	 signature= ('StefansExpressionSet' ),
	 function ( dat, groupcol='GroupName', colCol=NULL,  rowCol=NULL, probeNames=NULL, na.rm = FALSE, value.name = "value" ) {
	if ( is.null(probeNames)){
		probeNames <- dat@rownamescol
	}
	ma  <- dat@data[,order(dat@samples[,groupcol] )]
	rownames(ma) <- forceAbsoluteUniqueSample(as.vector(dat@annotation[, probeNames]) )
	melted <- reshape2::melt( cbind(rownames(ma),ma) )
	dat@samples <- dat@samples[order(dat@samples[,groupcol]),]
	if ( length( which ( melted[,2] == '') ) > 0 ){
		melted <- melted[ - which ( melted[,2] == ''),]
	}
	melted[,3] <- as.numeric(as.character(melted[,3]))
	grps <- NULL
	for ( i in as.vector(dat@samples[,groupcol]) ){
		grps <- c( grps, rep( i, nrow(dat@data)))
	}
	cgrps <- NULL
	for ( i in as.vector(dat@samples[,colCol]) ){
		cgrps <- c( cgrps, rep( i, nrow(dat@data)))
	}
	colnames(melted) <- c('ProbeName', 'SampleName', 'Expression')
	melted$Group <- grps
	melted <- addSampleColGroup( dat, melted, colName= colCol)
	melted <- addGeneColGroup( dat, melted, colName= rowCol)
	melted
})

#' @name addSampleColGroup
#' @aliases addSampleColGroup,StefansExpressionSet-method
#' @rdname addSampleColGroup-methods
#' @docType methods
#' @description an internal method adding new color columns to the melted dataset
#' @param x  the StefansExpressionSet object
#' @param melted the already melted matrix
#' @param colName the colnames to add to the matrix (defult=NULL)
#' @title description of function addSampleColGroup
setGeneric('addSampleColGroup', ## Name
		function ( x, melted, colName=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('addSampleColGroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('addSampleColGroup', signature = c ('StefansExpressionSet'),
		definition = function ( x, melted, colName=NULL ) {
			if ( ! is.null(colName) ) {
				datarows = nrow(x@data)
				grps <- NULL
				for ( GNid in 1:length(colName)){
					le <- datarows + GNid -1
					melted_new <- NULL
					for (sid in 1:datarows) {
						melted_new <- rbind( melted_new, as.matrix(melted[(1+le*(sid-1)):(le*sid),]) )
						line <- as.vector( t(melted[le*sid,]))
						melted_new <- rbind(melted_new,  matrix(c(colName[GNid] , as.character(line[2]),
												as.character(x@samples[sid,colName[GNid]]), as.character(line[4]) ), nrow=1) )
					}
					melted <- data.frame(melted_new,row.names= 1:nrow(melted_new))
				}
			}
			melted
		} )
#' @name addGeneColGroup
#' @aliases addGeneColGroup,StefansExpressionSet-method
#' @rdname addGeneColGroup-methods
#' @docType methods
#' @description an internal function to add row level clusters to a merged data table
#' @param x the StefansExpressionSet
#' @param melted the melted table
#' @param colName the annotation column names to use for the coloring default=NULL
#' @title description of function addGeneColGroup
setGeneric('addGeneColGroup', ## Name
		function ( x, melted, colName=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('addGeneColGroup') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('addGeneColGroup', signature = c ('StefansExpressionSet'),
		definition = function ( x, melted, colName=NULL ) {
			if ( ! is.null(colName)) {
				grps <- NULL
				datarows = nrow(x@data)
				ret <- NULL
				for ( GNid in 1:length(n)){
					add  <- as.matrix(melted[1:datarows,])
					add[,3] <- as.vector(x@annotation[,colName[GNid]])
					add[,4] <-'GeneGroup'
					ret <- rbind( ret, add )
				}
			}else {
				ret <- melted
			}
			ret
		} )




#' @name plot.stats.as.heatmap
#' @aliases plot.stats.as.heatmap,StefansExpressionSet-method
#' @rdname plot-methods
#' @docType methods
#' @description  create heatmap for each statistics table and each pvalue cutoff This function uses
#' @description  internally the plot.heatmaps() function for each selected probeset list
#' @param x the StefansExpressionSet
#' @param pvalue a vector of pvalues to set as cut off
#' @param Subset no idear at the moment....
#' @param comp a list of comparisons to restric the plotting to (NULL = all)
#' @param Subset.name no idear what that should do here....
#' @param gene_centered collapse all genes with the same gene symbol into one value
#' @param collaps how to collapse the data if gene_centered
#' @param geneNameCol the name of the gene.symbol column in the annotation
#' @title description of function plot
#' @export 
setGeneric('plot.stats.as.heatmap', ## Name
		function ( x, pvalue=c( 0.1,1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), Subset=NULL , Subset.name= NULL, comp=NULL, gene_centered=F, collaps=NULL,geneNameCol= "mgi_symbol") { ## Argumente der generischen Funktion
			standardGeneric('plot.stats.as.heatmap') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('plot.stats.as.heatmap', signature = c ('StefansExpressionSet'),
		definition = function ( x, pvalue=c( 0.1,1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), Subset=NULL , Subset.name= NULL, comp=NULL, gene_centered=F, collaps=NULL,geneNameCol= "mgi_symbol") {
			if ( !is.null(comp) ){
				print ("At the moment it is not possible to reduce the plotting to one comparison only" )
				return (x)
			}
			add = ''
			orig.name = x@name
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
				x@name = paste( add,x@name, sep='_')
			}
			print ( x@name )
			for (i in match( names(x@sig_genes), pvalue ) ){
				try( plot.heatmaps( 
								x, 
								x@sig_genes[[i]],
								names(x@sig_genes)[i],
								analysis_name = paste (x@name, version), 
								gene_centered = gene_centered,
								Subset = Subset,
								collaps = collaps,
								geneNameCol= geneNameCol
						), silent=F)
			}
			x@name = orig.name
		})

#' @name plot.probeset
#' @aliases plot.probeset,StefansExpressionSet-method
#' @rdname plot.probeset-methods
#' @docType methods
#' @description This function plots the expression data grouped by the GroupName using ggplot2.
#' @description Works only for one probeset.
#' @param x the StefansExpressionSet object
#' @param groupcol the grouping column in the samples data
#' @param colCol the coloring column in the sample data
#' @param probeNames the column in the annotation datacontaining the gene symbol
#' @param x the StefansExpressionSet object
#' @param probeset the probeset name (rownames(x@data))
#' @param boxplot (F or T) create a a dots- or box-plot
#' @param pdf save the file as pdf (default svg)
#' @param geneNameCol the column name for the gene symbol to use in the plots
#' @param sampleGroup the sample table column containing the grouping information
#' @title description of function plot.probeset
#' @examples 
#' p <- plot.probeset( red, rownames(red@data)[20], geneNameCol= 'GeneID', boxplot=T)
#' #produces the file '../../tmp/reducedSet_boxplot_Mybl1_expression.svg'
#' p <- plot.probeset( red, rownames(red@data)[20], geneNameCol= 'GeneID', boxplot=F)
#' #produces the file '../../tmp/reducedSet_points_Mybl1_expression.svg'
#' @export
setGeneric('plot.probeset', ## Name
	function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol", sampleGroup='GroupName' ) { ## Argumente der generischen Funktion
		standardGeneric('plot.probeset') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

####last

setMethod('plot.probeset', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, probeset, boxplot=F, pdf=F, geneNameCol= "mgi_symbol", sampleGroup='GroupName' ) {
	if ( sum(is.na(match(probeset, rownames(x@data)))==F) == 0 ){
		probeset <- rownames(x@data)[match( probeset, x@annotation[,geneNameCol] ) ]
	}
	if ( length(probeset) == 0 ) {
		stop( "I could not find the probeset in the data object!")
	}
	print ( probeset )
	x <- reduce.Obj(x, probeset )
	if ( nrow(x@data) == 1 ) {
		melted <- melt( x, probeNames=geneNameCol, groupcol = sampleGroup , colCol= sampleGroup )
		colnames(melted) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
		collaps <- '_points_'
		if (boxplot){
			collaps <-'_boxplot_'
		}
		if ( pdf ) {
			pdf( file=paste(x@outpath,x@name,collaps,x@annotation[,geneNameCol],"_expression.pdf",sep='') ,width=5,height=4)
		}else {
			devSVG( file=paste(x@outpath,x@name,collaps,x@annotation[,geneNameCol],"_expression.svg",sep='') ,width=5,height=4)
		}
		if ( boxplot ){
			p <- ggplot(melted ,aes(x=Group, y=Expression,color=Group)) + geom_boxplot(shape=19)+theme_bw()
		}else{
			#TODO: adjust jitter width based on sample group numbers
			p <-ggplot(melted ,aes(x=Group, y=Expression,color=Group)) +geom_jitter(shape=19, width=0.2, height=0)+theme_bw()
		}
		print ( p )
		dev.off()
		invisible(p)
	}
	else {
		stop("You should rather use a heatmap to display more than one probeset")
	}
})



#' @name drop.samples
#' @aliases drop.samples,StefansExpressionSet-method
#' @rdname drop.samples-methods
#' @docType methods
#' @description  drops samples from the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @param samplenames which samples to drop (samples like colnames(x@data))
#' @param name the name of the new StefansExpressionSet object
#' @title description of function drop.samples
#' @export 
setGeneric('drop.samples', ## Name
	function ( x, samplenames=NULL, name='dropped_samples' ) { ## Argumente der generischen Funktion
		standardGeneric('drop.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('drop.samples', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, samplenames=NULL, name='dropped_samples' ) {
	if ( ! is.null(samplenames)){
		red  <- new(class(x)[1], name=name )
		for (n in c( slot(x,'simple'), 'annotation') ){
			slot(red,n) <- slot(x,n)
		}
		red@samples <- x@samples[ is.na(match(x@samples[,x@sampleNamesCol], samplenames  ) ) == T ,]
		#print ( paste( "Dropping", length(samplenames), "samples (", paste( samplenames, collapse=", "),")") )
		
		red@data <- x@data[, make.names(as.vector(red@samples[,red@sampleNamesCol]))]
	#	colnames(red@data) <- forceAbsoluteUniqueSample ( as.vector(red@samples[, red@sampleNamesCol ]) )
	#	red@samples[,red@sampleNamesCol] <- colnames(red@data)
	}
	red
})

#' @name restrictSamples
#' @aliases restrictSamples,StefansExpressionSet-method
#' @rdname restrictSamples-methods
#' @docType methods
#' @description Drop the samples, that have been selected!
#' @param x the StefansExpressionSet object
#' @param name the name of the new StefansExpressionSet
#' @param column which column to analyze in the samples table
#' @param value which value to take as 'cut off'
#' @param mode one of 'less', 'more', 'onlyless', 'equals'
#' @title description of function restrictSamples
#' @export 
setGeneric('restrictSamples', ## Name
	function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ) { ## Argumente der generischen Funktion
		standardGeneric('restrictSamples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('restrictSamples', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ) {
		
	S <- NULL
	if ( is.null(value)) {
		stop( "the value must not be NULL!")
	}
	
	switch( mode,
			'less' = S <- x@samples[which ( x@samples[,column] <=  value), ], 
			'more' = S <- x@samples[which ( x@samples[,column] > value ), ], 
			'onlyless' = S <- x@samples[which ( x@samples[,column]  < value ), ],
			'equals' = S <- x@samples[which ( x@samples[,column] ==  value), ]
	)
	browser()
	if ( nrow(S) < nrow(x@samples)){
		x <- drop.samples( x, S[,x@sampleNamesCol], name=name)
	}
	write.table (S, file=paste(name,'_Sample_Description', ".xls", sep=''), sep='\t',  row.names=F,quote=F )
	x
})

####last

#' @name pwd
#' @aliases pwd
#' @rdname pwd-methods
#' @docType methods
#' @description  uses the linux pwd command to determin the working directory 
#' @return A string containing the working directory 
#' @title description of function pwd
#' @export 
setGeneric('pwd', ## Name
	function ( a ) { ## Argumente der generischen Funktion
		standardGeneric('pwd') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('pwd', signature = c () ,
	definition = function ( a ) {
		rm(a)
	system( 'pwd > __pwd' )
	t <- read.delim( file = '__pwd', header=F)
	t <- as.vector(t[1,1])
	t <- paste(t,"/",sep='')
	unlink( '__pwd')
	t
})

#' @name forceAbsoluteUniqueSample
#' @aliases forceAbsoluteUniqueSample,character-method
#' @rdname forceAbsoluteUniqueSample-methods
#' @docType methods
#' @description  force absolute unique names in a vector by adding _<amount of repeats> to each value
#' @description  if there is more than one repeat opf the value @exportMethod
#' @param x the StefansExpressionSet object
#' @param separator '_' or anything you want to set the separator to
#' @title description of function forceAbsoluteUniqueSample
setGeneric('forceAbsoluteUniqueSample', ## Name
	function ( x ,separator='_') { ## Argumente der generischen Funktion
		standardGeneric('forceAbsoluteUniqueSample') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('forceAbsoluteUniqueSample',
		, signature = c ( 'character') ,
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



setMethod('forceAbsoluteUniqueSample',
		, signature = c ( 'factor') ,
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

#' @name addAnnotation
#' @aliases addAnnotation,StefansExpressionSet-method
#' @rdname addAnnotation-methods
#' @docType methods
#' @description  add aditional annotation to the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @param mart the annotation table (data.frame or mart object)
#' @param mart.col which column corresponds to the rownames(x@data)
#' @title description of function addAnnotation
#' @export 
setGeneric('addAnnotation', ## Name
	function (x ,mart, mart.col='refseq_mrna') { ## Argumente der generischen Funktion
		standardGeneric('addAnnotation') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('addAnnotation', signature = c ( 'StefansExpressionSet') ,
	definition = function (x ,mart, mart.col='refseq_mrna') {
	if ( ! class(mart) == 'data.frame' ){
		x@annotation <- cbind(x@annotation, mart[match(rownames(x@data),mart[,mart.col] ), ] )
	}
	else {
		x@annotation <- mart[is.na(match(mart[,mart.col], rownames(x@data)))==F,]
		rownames(x@annotation) <- rownames(x@data)
	}
	x
})

#' @name getAnnotation4probesets
#' @aliases getAnnotation4probesets,StefansExpressionSet-method
#' @rdname getAnnotation4probesets-methods
#' @docType methods
#' @description  subsets the annotation table for a set of probesets
#' @param x the StefansExpressionSet object
#' @param probesets the vector of probesets to annotate (rownames(x@data))
#' @param colname a vector of colnames to annotate
#' @return a vector of annotation values
#' @title description of function getAnnotation4probesets
setGeneric('getAnnotation4probesets', ## Name
	function (x, probesets=c(), colname='Gene.Symbol' ) { ## Argumente der generischen Funktion
		standardGeneric('getAnnotation4probesets') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('getAnnotation4probesets', signature = c ( 'StefansExpressionSet') ,
	definition = function (x, probesets=c(), colname='Gene.Symbol' ) {
	as.vector(x@annotation[match( probesets, rownames(x@data) ), colname ])
})

#' @name ranks
#' @aliases ranks,StefansExpressionSet-method
#' @rdname ranks-methods
#' @docType methods
#' @description  creates a new ranks vector in the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @title description of function ranks
#' @export 
setGeneric('ranks', ## Name
	function (x ) { ## Argumente der generischen Funktion
		standardGeneric('ranks') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('ranks', signature = c ('StefansExpressionSet') ,
	definition = function (x ) {
	if ( ! exists ( 'ranks', where =x ) ){
		x@ranks <- apply( x@data,2,order)
		colnames( x@ranks ) <- colnames(x@data) 
		rownames( x@ranks ) <- rownames(x@data) 
	}
	x
})
#' @name reduce.Obj
#' @aliases reduce.Obj,StefansExpressionSet-method
#' @rdname reduce.Obj-methods
#' @docType methods
#' @description  reduces the dataset based on genes e.g. dropps genes from the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @param probeSets a list of probesets to reduce the data to
#' @param name the new StefansExpressionSet name
#' @title description of function reduce.Obj
#' @export 
setGeneric('reduce.Obj', ## Name
	function ( x, probeSets=c(), name="reducedSet" ) { ## Argumente der generischen Funktion
		standardGeneric('reduce.Obj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reduce.Obj', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, probeSets=c(), name="reducedSet" ) {
	retObj <- new(class(x)[1], name = name)
	useOnly <- match(probeSets, rownames(x@data))
	not.matched <- probeSets[is.na(useOnly)]
	if ( length(not.matched) > 0 ){
		print (paste('Problematic genes:', paste(not.matched,sep=', ')))
		probeSets <- probeSets[ ! is.na(useOnly)]
		useOnly <- useOnly[ ! is.na(useOnly) ]
	}
	for (n in slot(x,'simple')){
		slot(retObj,n) <- slot(x,n)
	}
	retObj@samples <- x@samples
	retObj@data <- data.frame( x@data[ useOnly ,] )
	rownames(retObj@data) <- probeSets
	colnames(retObj@data) <- colnames(x@data)
	retObj@annotation <- data.frame(x@annotation[useOnly,]) ## if I only have one column here
	colnames(retObj@annotation) <- colnames(x@annotation)
	if ( length( names(x@stats)) > 0){
		for ( i in 1:length(names(x@stats))){
			retObj@stats[[i]]= x@stats[[i]][ match(probeSets ,x@stats[[i]][,1] ),]
		}
		names(retObj@stats) <- names(x@stats)
	}
	if ( length(x@ranks) > 0 ){
		retObj@ranks = x@ranks[useOnly,]
	}
	if ( ncol( x@raw ) > 0 ) {
		retObj@raw = x@raw[useOnly,]
	}
	retObj
})



#' @name ggplot.gene
#' @aliases ggplot.gene,StefansExpressionSet-method
#' @rdname ggplot.gene-methods
#' @docType methods
#' @description  Plot one gene in the StefansExpressionSet as boxplot or points plot (using ggplot2)
#' @param dat the StefansExpressionSet object
#' @param gene the gene of interest
#' @param colrs the grouping colors for the x axis (samples)
#' @param groupCol the samples clumn that contains the grouping information
#' @param colCol the sample column that contains the color information
#' @title description of function ggplot.gene
#' @export 
setGeneric('ggplot.gene', ## Name
	function (dat,gene, colrs=NULL, groupCol='GroupID', colCol='GroupID', boxplot=F) { ## Argumente der generischen Funktion
		standardGeneric('ggplot.gene') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('ggplot.gene', signature = c ( 'StefansExpressionSet') ,
	definition = function (dat,gene, colrs=NULL, groupCol='GroupID', colCol='GroupID', boxplot=F) {
	not.in = 'NUKL'
	if ( is.null(colrs)){
		colrs = rainbow( length(unique(dat@samples[,colCol])))
	}
	g1 <- melt.StefansExpressionSet(reduce.Obj ( dat, gene, name=gene ), probeNames=dat@rownamescol, groupcol=groupCol,colCol=colCol)
	#g1 <- subset(dat, Gene.Symbol == gene)
	colnames(g1) <- c( 'Gene.Symbol', 'Sample', 'Expression', 'Group', 'ColorGroup' )
	g1$Group <- factor( g1$Group, levels=unique( as.character(g1$Group)))
	if ( length(g1) == 0){
		not.in = c( gene )
	}
	if ( boxplot ){
		list ( plot=ggplot2::ggplot(g1 ,ggplot2::aes(x=Group, y=Expression,color=Group)) 
						+ ggplot2::geom_boxplot(shape=19)+ggplot2::theme_bw() 
						+ ggplot2::scale_colour_manual( values= colrs, guide=FALSE )
						+ ggplot2::theme( axis.text.x= ggplot2::element_text(angle=90) )
						+ ggplot2::labs(title=gene, x=''),
				not.in = not.in
		)
		
	}else{
		list ( plot=ggplot2::ggplot(g1,ggplot2::aes(x=Group, y=Expression,color=Group)) 
						+ ggplot2::geom_jitter(shape=19)+ggplot2::theme_bw()
						+ ggplot2::scale_colour_manual( values= colrs, guide=FALSE )
						+ ggplot2::theme( axis.text.x= ggplot2::element_text(angle=90) )
						+ ggplot2::labs(title=gene, x=''),
				not.in = not.in
		)
	}
})

#' @name gg.heatmap.list
#' @aliases gg.heatmap.list,StefansExpressionSet-method
#' @rdname gg.heatmap.list-methods
#' @docType methods
#' @description  uses ggplot2 to plot heatmaps
#' @param dat the StefansExpressionSet object
#' @param glist a list of probesets to plot (or all)
#' @param colrs a list of colors for the sample level boxes (or rainbow colors)
#' @param groupCol the column group in the samples table that contains the grouping strings
#' @param colCol the column group in the samples table that contains the color groups
#' @title description of function gg.heatmap.list
#' @export 
setGeneric('gg.heatmap.list', ## Name
	function (dat,glist=NULL, colrs=NULL, groupCol='GroupID', colCol=NULL) { ## Argumente der generischen Funktion
		standardGeneric('gg.heatmap.list') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('gg.heatmap.list', signature = c ( 'StefansExpressionSet') ,
	definition = function (dat,glist=NULL, colrs=NULL, groupCol='GroupID', colCol=NULL) {
	
	if ( ! is.null(glist) ) {
		isect <- reduce.Obj ( dat, glist)
	}else {
		isect <- dat
	}
	if ( is.null(colCol)){
		colCol <- groupCol
	}
	if ( is.null(colrs) ){
		colrs = rainbow( length(unique(isect@samples[,colCol])))
	}
	if ( ! isect@zscored ) {isect <- z.score(isect)}
	dat.ss <- melt.StefansExpressionSet ( isect, probeNames=isect@rownamescol, groupcol=groupCol,colCol=colCol)
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
	samp.cast <- reshape2::dcast(dat.ss,Gene.Symbol~Sample,mean,value.var="z")
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
	brks= c( -20.1, as.vector(quantile(dat.ss$z[which(dat.ss$z != -20)],seq(0,1,by=0.1)) ))
	brks = unique(brks)
	print ( brks )
	brks[length(brks)] = brks[length(brks)] + 0.1
	dat.ss$z <- cut( dat.ss$z, breaks= brks)
	
	p = ( ggplot2::ggplot(dat.ss, ggplot2::aes(x=SampleName,y=ProbeName))
				+ ggplot2::geom_tile(ggplot2::aes(fill=z)) 
				+ ggplot2::scale_fill_manual( values =  c( 'gray', gplots::bluered(length(brks) -2  )) ) 
				+ ggplot2::theme(
						legend.position= 'bottom',
						axis.text.x=ggplot2::element_blank(),
#axis.ticks.x=element_line(color=ss$colrss),
						axis.ticks.length=ggplot2::unit(0.00,"cm")
				)+ labs( y='') )
	if ( ncol(dat.ss) == 5 ){
		p <- p + ggplot2::facet_grid( ColorGroup ~ Group,scales="free", space='free')
	}else if ( ncol(dat.ss) == 4 ) {
		p <- p + ggplot2::facet_grid( . ~ Group,scales="free", space='free')
	}
	list ( plot = p, not.in = setdiff( glist, rownames(isect@data)) )
	
#	list ( plot = ggplot2::ggplot(dat.ss, ggplot2::aes(x=Sample,y=Gene.Symbol))
#					+ ggplot2::geom_tile(ggplot2::aes(fill=z))
#					+ ggplot2::scale_fill_manual( values = c( 'gray', gplots::bluered(length(brks) -2  )) ) 
#					+ ggplot2::theme(
#							legend.position= 'bottom',
#							axis.text.x=ggplot2::element_blank(),
#							axis.ticks.x=ggplot2::element_line(color=ss$colrss),
#							axis.ticks.length=ggplot2::unit(0.6,"cm")
#					)
#					+ ggplot2::labs( y=''),
#			not.in = setdiff( glist, rownames(isect@data)) )
})


#' @name z.score
#' @aliases z.score.matrix,StefansExpressionSet-method
#' @rdname z.score.matrix-methods
#' @docType methods
#' @description  z score the matrix
#' @param m the matrix of column = samples and rows = genes or an StefansExpressionSet
#' @return the z scored matrix
#' @title description of function z.score
#' @export 
setGeneric('z.score', ## Name
		function (m) { ## Argumente der generischen Funktion
			standardGeneric('z.score') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('z.score', signature = c ('matrix'),
		definition = function (m ) {
			rn <- rownames( m )
			me <- apply( m, 1, mean )
			sd <- apply( m, 1, sd )
			sd[which(sd==0)] <- 1e-8
			m <- (m - me) /sd
			rownames(m) <- rn
			m
		})

setMethod('z.score',signature = c ('StefansExpressionSet'),
		definition = function (m) {
			if (! m@zscored ){
			#m$data <- z.score( as.matrix( m$data ))
			rn <- rownames( m@data )
			me <- apply( m@data, 1, mean )
			sd <- apply( m@data, 1, sd )
			sd[which(sd==0)] <- 1e-8
			m@data <- (m@data - me) /sd
			rownames(m@data) <- rn
			m@zscored = TRUE
			}
			m
		})

#' @name plot.heatmaps
#' @aliases plot.heatmaps,StefansExpressionSet-method
#' @rdname plot.heatmaps-methods
#' @docType methods
#' @description  A flexible heatmapping function that depends on heatmap.2 and allows for many data
#' @description  selection/conversion options
#' @param dataOBJ the StefansExpressionSet object
#' @param gene.names the rownames(dataObj@data) level gene names
#' @param pvalue an optional cut off value to select genes from the statistical results tables
#' @param analysis_name the name for the outfiles
#' @param gene_centered in case there are multiple probesets for each gene - sum the data or display each probset
#' @param Subset an optional list of strings matching to the gene symbols used to select genes of interest
#' @param collaps collaps the sample groups into a single column of the heatmap using one of collaps=c('median', 'mean')
#' @title description of function plot.heatmaps
setGeneric('plot.heatmaps', ## Name
	function ( dataOBJ, gene.names=NULL , pvalue=1, analysis_name =NULL, gene_centered = F, Subset=NULL, collaps=NULL,geneNameCol= "mgi_symbol", pdf=F,... ) { ## Argumente der generischen Funktion
		standardGeneric('plot.heatmaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('plot.heatmaps', signature = c ( 'StefansExpressionSet') ,
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
	if ( ! is.null(Subset) ) { ## 
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


#' @name groups.boxplot
#' @aliases groups.boxplot,StefansExpressionSet-method
#' @rdname groups.boxplot-methods
#' @docType methods
#' @description  This function can be used to get an overview of the different gene level groups in
#' @description  an StefansExpressionSet It uses the linux montage command to merge all different boxplots
#' @description  into one figure
#' @param x the ExpresionSet object
#' @param SampleCol the column in the samples table that contains the sample grouping information
#' @param clusters the gene level clusters. The function plots one boxplot for each cluster.
#' @param svg create the single boxplots as svg or (default png) files
#' @param fname the filename extension files are created by paste(x@outpath,fname,<GroupID>,"_boxplot")
#' @param Collapse if set will lead to a collapse of the sample groups into one value per gene. Supports all \code{\link{collaps}} by options.
#' @title description of function groups.boxplot
setGeneric('groups.boxplot', ## Name
	function ( x, SampleCol='GroupName', clusters, svg=F, fname='group_', width=800, height=800,mar=NULL, Collapse=NULL, ...) { ## Argumente der generischen Funktion
		standardGeneric('groups.boxplot') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('groups.boxplot', signature = c ( 'StefansExpressionSet') ,
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

#' @name collaps
#' @aliases collaps,StefansExpressionSet-method
#' @rdname collaps-methods
#' @docType methods
#' @description  This function will collpase the data in the StefansExpressionSet to only contain one value
#' @description  per sample group.
#' @param dataObj the StefansExpressionSet
#' @param by collapsing method c('median','mean','sd','sum', or own function )
#' @param groupCol the sample names you want to group on
#' @title description of function collaps
#' @export 
setGeneric('collaps', ## Name
		function (dataObj, by=c('median','mean','sd','sum' ), groupCol='GroupID' ) { ## Argumente der generischen Funktion
			standardGeneric('collaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('collaps', signature = c ('StefansExpressionSet'),
		definition = function (dataObj, by=c('median','mean','sd','sum', 'var' ), groupCol='GroupID' ) {
			u <- unique(as.vector(dataObj@samples[,groupCol]))
			m <- length(u)
			mm <-  matrix ( rep(0,m * nrow(dataObj@data)), ncol=m)
			colnames(mm) <- u
			rownames(mm) <- rownames(dataObj@data)
			f <- NULL
			if ( is.function(by)){
				f <- by
			}else {
			switch( by,
					median = f<- function (x ) { median(x) },
					mean = f <- function(x) { mean(x) },
					sd = f <- function(x) { sd(x) },
					sum = f <-function(x) { sum(x)},
					var = f <- function(x) { var(x) }
			);
			}
			if ( is.null(f) ) {
				stop("Please set what to one of 'median','mean','sd','sum'" )
			}
			new_samples <- NULL
			for ( i in u ){
				all <- which(as.vector(dataObj@samples[, groupCol]) == i )
				new_samples <- rbind ( new_samples, dataObj@samples[all[1],] )
				mm[,i] <- apply( dataObj@data[ , all],1,f)
			}
			name = paste( unlist(strsplit( paste( dataObj@name, groupCol, by, sep='_') , '\\s')) , collapse='_')
			try ( { ret <- drop.samples ( dataObj,
					samplenames= setdiff( 
							dataObj@samples[,dataObj@sampleNamesCol] , 
							new_samples[,dataObj@sampleNamesCol] 
					), 
					name= name
			) }, silent=TRUE )
			colnames(mm) <- as.vector(new_samples[, dataObj@sampleNamesCol])
			ret@data <- data.frame(mm)
			ret
})


#' @name simpleAnova
#' @aliases simpleAnova,StefansExpressionSet-method
#' @rdname simpleAnova-methods
#' @docType methods
#' @description  This function calculates an annova to identify significant changes in the StefansExpressionSet
#' @description  has a higher sensitivity for multi group analyses to identify group specific changes
#' @description  or general trends in the dataset. This function adds the results into the stats slot
#' @description  of the StefansExpressionSet object.
#' @param x the StefansExpressionSet object
#' @param groupCol the samples table column that contains the grouping information
#' @param padjMethod the p value correction method as described in  \code{\link[stats]{p.adjust}}
#' @title description of function simpleAnova
#' @export 
setGeneric('simpleAnova', ## Name
	function ( x, groupCol='GroupName', padjMethod='BH' ) { ## Argumente der generischen Funktion
		standardGeneric('simpleAnova') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('simpleAnova', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, groupCol='GroupName', padjMethod='BH' ) {
	x <- normalize(x)
	significants <- apply ( x@data ,1, function(a) { anova( lm (a ~ a@Samples[, samples.col]))$"Pr(>F)"[1] } )
	adj.p <- p.adjust( significants, method = padjMethod)
	res <- cbind(significants,adj.p )
	res <- data.frame(cbind( rownames(res), res ))
	colnames(res) <- c('genes', 'pvalue', paste('padj',padjMethod) )
	if ( length (x@stats) == 0 ){
		x@stats <- list ( 'simpleAnova' = res )
	}
	else {
		x@stats[[length(x@stats)+1]] <- res
		names(x@stats)[length(x@stats)] = 'simpleAnova'
	}
	x
})

#' @name createStats
#' @aliases createStats,StefansExpressionSet-method
#' @rdname createStats-methods
#' @docType methods
#' @description  constructor that has to be implemented in the data specific classes
#' @param x the StefansExpressionSet object
#' @param condition the samples column containing the condition of interest
#' @param files whether or not (default) export the stat files using \code{\link{writeStatTables}}
#' @param A an optional condition A to compare to a condition B
#' @param B the condition B
#' @title description of function createStats
setGeneric('createStats', ## Name
	function ( x, condition, files=F, A=NULL, B=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('createStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('createStats', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, condition, files=F, A=NULL, B=NULL ) {
	stop( 'Not implemented' )
})

#' @name normalize
#' @aliases normalize,StefansExpressionSet-method
#' @rdname normalize-methods
#' @docType methods
#' @description  constructor that has to be implemented in the data specific classes
#' @param x the StefansExpressionSet object
#' @title description of function normalize
setGeneric('normalize', ## Name
	function ( x ) { ## Argumente der generischen Funktion
		standardGeneric('normalize') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('normalize', signature = c ('StefansExpressionSet') ,
	definition = function ( x ) {
	x
})

#' @name force.numeric
#' @aliases force.numeric,StefansExpressionSet-method
#' @rdname force.numeric-methods
#' @docType methods
#' @description  The moethod forces the values in the data matrix to be numbers.
#' @param dataObj the StefansExpressionSet object
#' @title description of function force.numeric
#' @export 
setGeneric('force.numeric', ## Name
	function (dataObj ) { ## Argumente der generischen Funktion
		standardGeneric('force.numeric') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('force.numeric', signature = c ('StefansExpressionSet') ,
	definition = function ( dataObj ) {
	for ( i in 1: ncol(dataObj@data) ) { 
		if ( !  paste(as.vector(dataObj@data[,i]), collapse=" ") == paste(as.vector(as.numeric(dataObj@data[,i])), collapse=" ") ) { 
			dataObj@data[,i] <- as.vector(as.numeric(dataObj@data[,i]))
		}
	}
	dataObj
})

#' @name show
#' @aliases show,StefansExpressionSet-method
#' @rdname show-methods
#' @docType methods
#' @description  print the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @return nothing
#' @title description of function show
#' @export 
setMethod('show', signature(object='StefansExpressionSet') ,
	definition = function (object) {
	cat (paste("An object of class", class(object)),"\n" )
	cat("named ",object@name,"\n")
	cat (paste( 'with',nrow(object@data),'genes and', ncol(object@data),' samples.'),"\n")
	cat (paste("Annotation datasets (",paste(dim(object@annotation),collapse=','),"): '",paste( colnames(object@annotation ), collapse="', '"),"'  ",sep='' ),"\n")
	cat (paste("Sample annotation (",paste(dim(object@samples),collapse=','),"): '",paste( colnames(object@samples ), collapse="', '"),"'  ",sep='' ),"\n")
	if ( length(names(object@stats)) > 0 ){
		cat ( "P values were calculated for ", length(names(object@stats)) -1, " condition(s)\n")
	}
})


#' @name export.data
#' @aliases export.data,StefansExpressionSet-method
#' @rdname export.data-methods
#' @docType methods
#' @description  write the StefansExpressionSet data, annotation and samples tables to disk
#' @param x the StefansExpressionSet object
#' @title description of function write.data
#' @export 
setGeneric('export.data', ## Name
		function ( x ) { ## Argumente der generischen Funktion
			standardGeneric('export.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('export.data', signature = c ( 'StefansExpressionSet') ,
		definition = function ( x ) {
			write.table(cbind ( x@annotation, x@data ), file=paste(x@outpath,x@name,"_expressionValues.xls",sep=''),sep='\t', row.names=F,quote=F )
			write.table( x@samples, file=paste( x@outpath, x@name, "_sampleInformation.xls",sep='' ), sep='\t', row.names=F,quote=F )
		})

#' @name write.data
#' @aliases write.data,StefansExpressionSet-method
#' @rdname write.data-methods
#' @docType methods
#' @description  write the StefansExpressionSet data table to disk
#' @param x the StefansExpressionSet object
#' @param annotation a vector of annotation data column names to include in the written table (default=none)
#' @title description of function write.data
#' @export 
setGeneric('write.data', ## Name
	function ( x, annotation=NULL ) { ## Argumente der generischen Funktion
		standardGeneric('write.data') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('write.data', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, annotation=NULL ) {
	if ( !is.null(annotation) ) {
		write.table( cbind( x@annotation[,annotation], x@data), file= paste( x@outpath,x@name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
	else {
		write.table( cbind( rownames(x@data), x@data), file= paste( x@outpath,x@name,"_expressionValues.xls",sep=''), row.names=F, sep="\t",quote=F )
	}
})

#' @name writeStatTables
#' @aliases writeStatTables,StefansExpressionSet-method
#' @rdname writeStatTables-methods
#' @docType methods
#' @description  export the statistic files
#' @param x the StefansExpressionSet
#' @title description of function writeStatTables
#' @export 
setGeneric('writeStatTables', ## Name
	function ( x, annotation=NULL, wData=F ) { ## Argumente der generischen Funktion
		standardGeneric('writeStatTables') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('writeStatTables', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, annotation=NULL, wData=F ) {
	opath = paste(x@outpath,'stats_Pval/', sep='')
	if ( wData ) {
		opath = paste( x@outpath,'stats_wData/',sep='')
	}
	write.table( x@samples, file= paste(opath,'stats_Sample_Description.xls', sep=''),sep='\t', row.names=F,quote=F  )
	Comparisons <- make.names(names(x@stats))
	system ( paste('mkdir ',opath, sep='') )
	for ( i in 1:length(Comparisons )){
		fname <- paste(opath ,Comparisons[i], '.xls',sep='')
		if ( ! is.null(annotation) ) {
			if ( wData==F ) {
				write.table( cbind( x@annotaion[, annotation],  x@stats[[i]] ) , file= fname, row.names=F, sep="\t",quote=F )
			}else {
				write.table( cbind(x@annotation[,annotation], x@stats[[i]], x@data ), file= fname, row.names=F, sep="\t",quote=F )
			}
		}
		else {
			if ( wData==F ) {
				write.table( x@stats[[i]], file= fname, row.names=F, sep="\t",quote=F )
			}else {
				write.table( cbind(x@stats[[i]], x@data ), file= fname, row.names=F, sep="\t",quote=F )
			}
		}
		print ( paste ( "table ",fname," for cmp ",Comparisons[i],"written" ) )
	}
})

#' @name export4GEDI
#' @aliases export4GEDI,StefansExpressionSet-method
#' @rdname export4GEDI-methods
#' @docType methods
#' @description  Convert the values in the StefansExpressionSet to the format GEDI program can import.
#' @param x the StefansExpressionSet object
#' @param fname the filename to export the data to
#' @param tag.col the sample name column in the samples table
#' @param time.col the time column in the samples table
#' @seealso \url{"http://apps.childrenshospital.org/clinical/research/ingber/GEDI/gedihome.htm"}
#' @title description of function export4GEDI
setGeneric('export4GEDI', ## Name
	function ( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) { ## Argumente der generischen Funktion
		standardGeneric('export4GEDI') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('export4GEDI', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) {
	if ( is.null(time.col)){
		stop ( paste( "choose time.col from:", paste( colnames(x@samples), collapse=", ") ) ) 
	}
	if ( is.null(tag.col)){
		stop ( paste( "choose tag.col from:", paste( colnames(x@samples), collapse=", ") ) ) 
	}
	groupnames <- vector('numeric', nrow(x@samples))
	for (i in 1:nrow(x@samples)) {
		groupnames[i] = paste( x@samples[i, tag.col], x@samples[i, time.col] , sep="_" )
	}
	if ( length(which(table(groupnames) > 1)) > 0 ){
		stop ( "Sorry, please calculate the mean expression for the data first ('collaps()')") 
	}
	
	treatments <- unique( as.vector(x@samples[ , tag.col] ))
	
	## now I need to find all possible days
	possible <- NULL
	for ( t in treatments ){
		possible <- c( possible, x@samples[ grep( t, x@samples[ , tag.col]) , time.col] )
	}
	required = names(which(table(possible) > minSample_PerTime) )
	passed <- NULL
	for ( t in treatments ){
		l <- rep ( -1, nrow(x@samples) )
		p <- grep( t, x@samples[ , tag.col])
		if ( length(p) >= length(required) ){
			l[p] <- x@samples[ p, time.col]
			l[ is.na(match(l, required))== T ] <- -1
			passed <- c(passed, paste( c(paste("}",t,sep=""),  l), collapse="\t") )
		}
	}
	
	fileConn<-file(fname)
	open( fileConn, "w" )
	writeLines( paste( "}Dynamic", length(passed), nrow(x@data), x@name, sep="\t"), con=fileConn)
	writeLines( paste(passed , collapse="\n") , con=fileConn )
	
	for ( i in 1:nrow(x@data) ){
		writeLines( paste( c(rownames(x@data)[i], x@data[i,]), collapse="\t" ),  con=fileConn )
	}
	
	close(fileConn)
	print( paste ("created file", fname))
})
