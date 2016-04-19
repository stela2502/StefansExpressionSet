
#' @name getProbesetsFromValues
#' @aliases getProbesetsFromValues,StefansExpressionSet-method
#' @rdname getProbesetsFromValues-methods
#' @docType methods
#' @description Using this function you can select probes based on there expression value in a set of samples
#' @param x the ExpreesionSet object
#' @param v the cut off- or matching value
#' @param sample sample names as used in data column names
#' @param mode one of less, more, onlyless or equals
#' @title description of function getProbesetsFromValues
#' @return a list of probesets that match the requirements
#' @export 
setGeneric('getProbesetsFromValues', ## Name
		function ( x, v='NULL', sample='NULL', mode='less' ) { ## Argumente der generischen Funktion
			standardGeneric('getProbesetsFromValues') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('getProbesetsFromValues', signature = c ('StefansExpressionSet'),
		definition = function ( x, v='NULL', sample='NULL', mode='less' ) {
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
						'less' = probesets <- c(probesets, as.vector(rownames(x@data)[which(x@data[,s] <= v)] ) ) ,
						'more' = probesets <- c(probesets, as.vector(rownames(x@data)[which(x@data[,s] > v)] ) ), 
						'onlyless' = probesets <- c(probesets,  as.vector(rownames(x@data)[which(x@data[,s] < v)] ) ),
						'equals' = probesets <- c(probesets, as.vector(rownames(x@data)[which(x@data[,s] == v)] ) )
				)
			}
			unique(probesets)
		})


#' @name name_4_IDs
#' @aliases name_4_IDs,StefansExpressionSet-method
#' @rdname name_4_IDs-methods
#' @docType methods
#' @description  Select probesets, that show a certain level in expression for a single sample probes
#' @param v The cutoff value
#' @param sample The sample name
#' @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
#' @return a list of probesets that has a expression less than 10 in sample A
#' @title description of function name_4_IDs
#' @export 
setGeneric('name_4_IDs', ## Name
		function ( x, ids=NULL, geneNameCol='mgi_symbol' ) { ## Argumente der generischen Funktion
			standardGeneric('name_4_IDs') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('name_4_IDs', signature = c ('StefansExpressionSet'),
		definition = function ( x, ids=NULL, geneNameCol='mgi_symbol' ) {
			if ( is.null(ids) ) {
				ids <- as.vector(colnames(x@data) )
			}
			as.vector(x@annotation[match( ids,x@annotation[, x@rownamescol]),geneNameCol])
		})


#' @name get_gene_list
#' @aliases get_gene_list,StefansExpressionSet-method
#' @rdname get_gene_list-methods
#' @docType methods
#' @description Query all stat tables and select genes based on multiple p value cutoff values
#' @param x the ExpressionSet
#' @param p_value a list of p value cut offs (0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 )
#' @param geneNameCol the (optional) annotation column defining the gene names
#' @title description of function get_gene_list
#' @return the ExpressionSet object with values in the sig_genes slot
#' @export 
setGeneric('get_gene_list', ## Name
		function (x, p_value = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), geneNameCol=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('get_gene_list') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('get_gene_list', signature = c ('StefansExpressionSet'),
		definition = function (x, p_value = c ( 0.1,  1e-2 ,1e-3,1e-4,1e-5, 1e-6, 1e-7, 1e-8 ), geneNameCol=NULL ) {
			if ( is.null(geneNameCol)){
				geneNameCol = x@rownamescol
			}
			if ( exists(where=x, 'stats')){
				cmps <- names(x@stats)
				p_value <- as.numeric( p_value )
				x@sig_genes <- vector ('list', length(p_value))
				names( x@sig_genes ) <- as.vector(p_value)
				for ( p in 1:length(p_value)){
					x@sig_genes[[p]] <- vector ('list', length(cmps)-1)
					for ( i in 2:length(cmps) ){
						sig.genes <- x@stats[[i]][which(x@stats[[i]]$padj < p_value[p] ), ] 
						sig.names <- name_4_IDs.StefansExpressionSet( x, sig.genes[,1], geneNameCol)
						sig.tab <- cbind(sig.names,sig.genes ) 
						if ( ncol(sig.tab) > 0 ){
							write.table(sig.tab,paste(x@outpath,x@name,'_',cmps[i],p_value[p] ,".xls",sep=''),col.names=T,row.names=F,sep="\t",quote=F)
						}
						x@sig_genes[[p]][[i-1]] = list (id = sig.genes[,1], names=sig.names )
					}
					x@sig_genes[[p]]$all <- list( id = unique(unlist(lapply (x@sig_genes[[p]] , function(a) { a$id } ))), names= unique(unlist(lapply ( x@sig_genes[[p]], function(a) { a$names } ))) )
				}
			}
			x
		})


#' @name describe.probeset
#' @aliases describe.probeset,StefansExpressionSet-method
#' @rdname describe.probeset-methods
#' @docType methods
#' @description print all annotation values + all values + all stat results for one probeset
#' @param x the ExpressionSet object
#' @param probeset the probeset to describe
#' @title description of function describe.probeset
#' @export 
setGeneric('describe.probeset', ## Name
		function ( x, probeset ) { ## Argumente der generischen Funktion
			standardGeneric('describe.probeset') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('describe.probeset', signature = c ('StefansExpressionSet'),
		definition = function ( x, probeset ) {
			ret <- list()
			print ( paste ("Annoataion for probeset ",probeset ) )
			ret$annotation <- x@annotation[ probeset, ] 
			print ( ret$annotation )
			print ( "Values:" ) 
			ret$data <- x@data[probeset, ]
			print ( ret$data )
			ret$stats <- NULL
			#browser()
			if ( length(x@stats)>0 ) {
				for( i in 1:length(x@stats) ) {
					#		browser()
					ret$stats <- rbind(ret$stats, 
							cbind ( 
									rep(names(x@stats)[i],length(probeset) ), 
									x@stats[[i]][is.na(match( x@stats[[i]][,1], probeset))==F,] 
							) 
					) 
				}
				colnames(ret$stats) <- c('comparison', colnames(x@stats[[1]]) )
				print ( ret$stats )
			}
			invisible(ret)
		})


#' @name getProbesetsFromStats
#' @aliases getProbesetsFromStats,StefansExpressionSet-method
#' @rdname getProbesetsFromStats-methods
#' @docType methods
#' @description  getProbesetsFromStats returns a list of probesets (the rownames from the data matrix)
#' @description  for a restriction of a list of stat comparisons probes <- getProbesetsFromStats (
#' @description  x, v=1e-4, pos="adj.P.Val" )
#' @param v The cutoff value
#' @param pos The column in the stats tables to base the selection on
#' @param Comparisons A list of comparisons to check (all if left out)
#' @param mode one of 'less', 'more', 'onlyless' or 'equals' default 'less' ( <= )
#' @return a list of probesets that shows an adjusted p value below 1e-4
#' @title description of function getProbesetsFromStats
#' @export 
setGeneric('getProbesetsFromStats', ## Name
		function ( x, v=1e-4, pos='padj', mode='less', Comparisons=NULL ) { ## Argumente der generischen Funktion
			standardGeneric('getProbesetsFromStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('getProbesetsFromStats', signature = c ('StefansExpressionSet'),
		definition = function ( x, v=1e-4, pos='padj', mode='less', Comparisons=NULL ) {
			if ( is.null(Comparisons)){	Comparisons <- names(x@stats) }
			probesets <- NULL
			for ( i in match(Comparisons, names(x@stats) ) ) {
				switch( mode,
						'less' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] <= v),1] )),
						'more' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] > v),1] )), 
						'onlyless' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] < v),1] )),
						'equals' = probesets <- c( probesets, as.vector(x@stats[[i]][which(x@stats[[i]][,pos] == v),1] ))
				)
			}
			unique(probesets)
		})

