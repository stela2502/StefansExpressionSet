#' @name as_NGSexpressionSet
#' @alias.NGSexpressionSetes as.NGSexpressionSet,NGSexpressionSet-method
#' @rdname as_NGSexpressionSet-methods
#' @docType methods
#' @description create a NGSexpressionSet from a counts list object
#' @param dat the counts list you get from Rsubread::featureCounts()
#' @title description of function as_NGSexpressionSet
#' @export 
setGeneric('as_NGSexpressionSet', ## Name
		function ( dat ) { ## Argumente der generischen Funktion
			standardGeneric('as_NGSexpressionSet') ## der Aufruf von standardGeneric sorgt für das.NGSexpressionSet Dispatching
		}
)

setMethod('as_NGSexpressionSet', signature = c ('list'),
		definition = function ( dat ) {
			ret = NULL
			if (all.equal( names ( dat), c("counts" ,"annotation", "targets", "stat")  ) ) {
				samples <- data.frame(t(dat$stat))
				colnames(samples) <- as.NGSexpressionSet.character(t(samples[1,]))
				samples$filename <- rownames(samples)
				rownames(samples) <-1:nrow(samples)
				ret <- NGSexpressionSet->new( 
						dat= cbind(dat$annotation, dat$counts), 
						samples = samples, 
						namecol= 'filename', 
						namerow= 'GeneID',
						outpath= ''
				)
			}
			else {
				print ("The list needs to contain the entries counts ,annotation, targets and stat" )
			}
			ret
		} )
