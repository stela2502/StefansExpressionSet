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

