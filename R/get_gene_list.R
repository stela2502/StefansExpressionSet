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


