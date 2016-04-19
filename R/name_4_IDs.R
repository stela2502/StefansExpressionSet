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


