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
		}
)

#' @name simpleAnova
#' @aliases simpleAnova,SingleCellsNGS-method
#' @rdname simpleAnova-methods
#' @docType methods
#' @description  This function calculates an annova to identify significant changes in the StefansExpressionSet
#' @description  has a higher sensitivity for multi group analyses to identify group specific changes
#' @description  or general trends in the dataset. This function adds the results into the stats slot
#' @description  of the SingleCellsNGS object. In contrast to the StefansExpressionSet version of the function,
#' @description the cells showing no expression are axcluded from the stats.
#' @param x the SingleCellsNGS object
#' @param groupCol the samples table column that contains the grouping information
#' @param padjMethod the p value correction method as described in  \code{\link[stats]{p.adjust}}
#' @title description of function simpleAnova
#' @export 
setMethod('simpleAnova', signature = c ( 'SingleCellsNGS') ,
		definition = function ( x, groupCol='GroupName', padjMethod='BH' ) {
			if ( is.null(x@stats[[paste('simpleAnova', groupCol)]] )) {
				x <- normalize(x)
				significants <- apply ( x@data ,1, function(a) {
							ids <- which(a > 0 )
							not <- names(which (table(x@samples[ids,groupCol ]) < 10 ))
							ids <- ids[ is.na(match(x@samples[ids,groupCol], not))==T]
							if ( length(table(x@samples[ids,groupCol ]) ) > 1 ) {
								try(anova( lm (a[ids] ~ x@samples[ids,groupCol ]))$"Pr(>F)"[1]) 
							}
							else {
								1
							}
						} )
				adj.p <- p.adjust( significants, method = padjMethod)
				res <- cbind(significants,adj.p )
				res <- data.frame(cbind( rownames(res), res ))
				colnames(res) <- c('genes', 'pvalue', paste('padj',padjMethod) )
				res[,2] <- as.numeric(as.vector(res[,2]))
				res[,3] <- as.numeric(as.vector(res[,3]))
				if ( length (x@stats) == 0 ){
					x@stats <- list ( 'simpleAnova' = res )
				}
				else {
					x@stats[[length(x@stats)+1]] <- res
					names(x@stats)[length(x@stats)] = paste('simpleAnova', groupCol)
				}
			}
			x
		}
)