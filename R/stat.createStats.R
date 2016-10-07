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
	function (x, condition, files=F, A=NULL, B=NULL) { ## Argumente der generischen Funktion
		standardGeneric('createStats') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('createStats', signature = c ('NGSexpressionSet'),
	definition = function (x, condition, files=F, A=NULL, B=NULL) {
		#stop( "Not implemented / broken!")
		
		if ( nrow(x@data) < 1e+3 ) {
			stop ( "Please calculate the statistics only for the whole dataset!" )
		}
		if ( is.na( match ( condition, colnames(x@samples))) ) {
			stop ( 'Please select a condition from the sample colnames' )
		}
		id <- match ( condition, names(x@usedObj[['cds']]))
		if ( is.na(id)) { id <- 1 } 
		x<- preprocess ( x, condition )
		conditions <- as.vector(unique(x@samples[,condition]))
		if ( ! is.null(A) && ! is.null(B)) {
			x <- add_to_stat ( x, 
				stat = DESeq::nbinomTest(x@usedObj[['cds']][[id]], A, B ), 
				name = paste( A, B ,sep=' vs. ')
			)
		}
		else {
			for ( i in 1:(length(conditions)-1) ){
				for ( a in (i+1):length(conditions) ){
					x <- add_to_stat ( x, 
						stat = DESeq::nbinomTest(x@usedObj[['cds']][[id]], conditions[i] , conditions[a] ), 
						name = paste( conditions[i], conditions[a],sep=' vs. ')
					)

				}
			}
		}
		if ( files ) {
			writeStatTables( x )
		}
		x
})

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
			toM <- function (x) {
				d <- as.matrix(x@data)
				d[which(d==-20)] <- NA
				d[is.na(d)] <- 0
			}
			if ( is.null(x@samples[,condition]) ) {
				stop( paste("the condition",condition, "is not defined in the samples table!"))
			}
			if ( is.null(A) ) {
				name = paste ("SingleCellAssay",condition)
			}else {
				keep <- which( x@samples[,condition] ==A | x@samples[,condition] == B)
				name = paste ("SingleCellAssay",condition, A, B)
				x <- drop.samples ( x, colnames(x@data)[-keep], name=name)
			}
			browser()
			d <- toM(x)
			sca <- FromMatrix(class='SingleCellAssay', exprsArray=t(d), data.frame(wellKey=rownames(d)), data.frame(primerid=colnames(d)) )
			
			groups <- colData(sca)$GroupName <- x@samples[,condition]
			zlm.output <- zlm.SingleCellAssay(~ GroupName, sca, method='glm', ebayes=T)
			zlm.lr <- lrTest(zlm.output,'GroupName')
			x <- add_to_stat ( x, zlm.lr, name )
			x
			
		})
