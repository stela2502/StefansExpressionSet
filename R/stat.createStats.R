#' @name createStats
#' @aliases createStats,NGSexpressionSet-method
#' @rdname createStats-methods
#' @docType methods
#' @description  calculate staistics on all possible groupings using the DEseq nbinomTest test Both
#' @description  together create the group 'A vs. B'
#' @param x the NGSexpressionSet
#' @param condition the grouping column in the samples data.frame
#' @param files write the statistics tables (FALSE)
#' @param A the first component to analyze (Group A)
#' @param B the second component to analyze (Group B)
#' @return the NGSexpressionSet with a set of ststs tables
#' @title description of function createStats
#' @export
setGeneric('createStats', ## Name
	function (x, condition, ..., files=F, A=NULL, B=NULL) { ## Argumente der generischen Funktion
		standardGeneric('createStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('createStats', signature = c ('NGSexpressionSet'),
	definition = function (x, condition, files=F, A=NULL, B=NULL, lib='DESeq2' ) {
		#stop( "Not implemented / broken!")
		if ( nrow(x@data) < 1e+3 ) {
			stop ( "Please calculate the statistics only for the whole dataset!" )
		}
		if ( is.na( match ( condition, colnames(x@samples))) ) {
			stop ( 'Please select a condition from the sample colnames' )
		}
		if ( lib == 'DESeq2'){
			if (!library("DESeq2", quietly = TRUE,logical.return=TRUE )) {
					stop("package 'DESeq2' needed for this function to work. Please install it.",
							call. = FALSE)
			}
			print(paste("Creating stats for", condition, "using the DESeq2 package"))
			id <- match ( condition, colnames(x@samples))
			if ( is.na(id)) { id <- 1 }
			if ( is.null(x@usedObj$cds[[condition]])){
				x<- preprocess ( x, condition )
				x@usedObj$cds[[condition]] <- DESeq2::DESeq(x@usedObj$cds[[condition]], test='Wald')
			}
			conditions <- as.vector(unique(x@samples[,condition]))
			if ( ! is.null(A) && ! is.null(B)){
				x <- add_to_stat ( x,
					stat = as.data.frame(DESeq2::results(x@usedObj$cds[[condition]] ,c(condition, A, B ) )),
					name = paste( 'DESeq2', paste( A, B ,sep=' vs. '))
				)
			}else {
				for ( i in 1:(length(conditions)-1) ){
					for ( a in (i+1):length(conditions) ){
						x <- add_to_stat ( x,
							stat = as.data.frame(DESeq2::results(x@usedObj$cds[[condition]] , c(condition, conditions[i] , conditions[a]) )),
							name = paste( 'DESeq2', paste( conditions[i], conditions[a],sep=' vs. '))
						)
					}
				}
			}
		}
		else if ( lib == 'ROTS' ) {
			if (!library("ROTS", quietly = TRUE,logical.return=TRUE )) {
				stop("package 'ROTS' needed for this function to work. Please install it.",
						call. = FALSE)
		  }
			if ( ! is.null(A) && ! is.null(B)){
				tmp <- run_ROTS ( x, condition, A, B )
				x <- add_to_stat ( x,
					stat = data.frame(
						'GeneID' = rownames(x@data),
						'pvalue' =tmp$pvalue,
						'FDR' = tmp$FDR,
						'logfc' =tmp$logfc,
						'd' = tmp$d
						),
					name = paste( 'ROTS', paste( A, B ,sep=' vs. '))
				)
			}else {
				conditions <- as.vector(unique(as.character(x@samples[,condition])))
				for ( i in 1:(length(conditions)-1) ){
					for ( a in (i+1):length(conditions) ){
						tmp <- run_ROTS ( x, condition, conditions[i], conditions[a] )
						x <- add_to_stat ( x,
							stat = data.frame(
								'GeneID' = rownames(x@data),
								'pvalue' =tmp$pvalue,
								'FDR' = tmp$FDR,
								'logfc' =tmp$logfc,
								'd' = tmp$d
								),
							name = paste( 'DESeq2', paste( conditions[i], conditions[a],sep=' vs. '))
						)
					}
				}
			}
		}else {
			stop( paste(sep='', "please select a lib from 'DESeq2' or 'ROTS', not '", lib,"'") )
		}

		if ( files ) {
			writeStatTables( x )
		}
		x
})

run_ROTS <- function ( x, condition, A, B ) {
	x <- drop.samples ( x, samplenames= colnames(x@data)[
			intersect(
				grep ( A, x@samples[,condition], invert=T ) ,
				grep ( B, x@samples[,condition], invert=T )
				)
		] )
		ROTS::ROTS(
			apply(data@data,2,log2),
			as.numeric(	factor(
				as.character(x@samples[,condition]), levels=c( A, B )
				)
			) -1,
			B = 1000,
			K=500
		)
}

add_to_stat <- function( x, stat, name ) {
	if ( ! is.na( match( name, names(x@stats)))){
		x@stats[[ match( name, names(x@stats)) ]] <- stat
	}else {
		x@stats[[ length( x@stats ) +1 ]] <- stat
		names(x@stats)[length(x@stats) ] <- name
	}
	x
}


setMethod('createStats', signature = c ( 'StefansExpressionSet') ,
		definition = function ( x, condition, files=F, A=NULL, B=NULL ) {
			stop( 'Not implemented' )
		})


setMethod('createStats', signature = c ( 'SingleCellsNGS') ,
		definition = function ( x, condition, files=F, A=NULL, B=NULL ) {
			if (!library("MAST", quietly = TRUE,  logical.return =T) ) {
				stop("MAST needed for this function to work. Please install it.",
						call. = FALSE)
			}
			toM <- function (x) {
				d <- as.matrix(x@data)
				d[which(d==-20)] <- NA
				d[is.na(d)] <- 0
				d
			}
			if ( is.null(x@samples[,condition]) ) {
				stop( paste("the condition",condition, "is not defined in the samples table!"))
			}
			if ( is.null(A) ) {
				name = paste ("SingleCellAssay",condition)
				a <- x
			}else {
				keep <- which( x@samples[,condition] ==A | x@samples[,condition] == B)
				name = paste ("SingleCellAssay",condition, A, B)
				a <- drop.samples ( x, colnames(x@data)[-keep], name=name)
			}
			d <- toM(a)
			sca <- MAST::FromMatrix(class='SingleCellAssay',
					exprsArray=d,
					data.frame(wellKey=colnames(d), GroupName = a@samples[,condition]),
					data.frame(primerid=rownames(d)))

			groups <- colData(sca)$GroupName <- a@samples[,condition]
			zlm.output <- MAST::zlm.SingleCellAssay(~ GroupName, sca, method='glm', ebayes=T)
			zlm.lr <- MAST::lrTest(zlm.output,'GroupName')
			x <- add_to_stat ( x, zlm.lr[,,'Pr(>Chisq)'], name )
			x

		})
