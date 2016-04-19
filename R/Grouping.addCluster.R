#' @name addCluster
#' @aliases addCluster,Clusters-method
#' @rdname addCluster-methods
#' @docType methods
#' @description this function add new column to samples and annotation data NOT checking the column order
#' @param x  the StefansExpressionSet object
#' @param userGroups  the userGroups table 
#' @param col the column in the userGroups table that contains the group ids default='groupID'
#' @param what one of 'cells' or 'genes'. Where to add the group.  default= 'cells'
#' @title description of function addCluster
#' @export 
setGeneric('addCluster', ## Name
	function (x, userGroups, col='groupID', what= 'cells' ) { ## Argumente der generischen Funktion
		standardGeneric('addCluster') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('addCluster', signature = c ('StefansExpressionSet'),
	definition = function (x, userGroups, col='groupID', what= 'cells' ) {
	if ( what == 'cells' ) {
		cbind ( x@samples, userGroups[match( x@samples[,x@sampleNamesCol], userGroups[,1] ),col] )
		colnames(x@samples)[ncol(x@samples)] = col;
	}
	else if ( what == 'genes' ) {
		cbind ( x@annotation, userGroups[match( x@annotation[,x@rownamescol], userGroups[,1] ),col] )
		colnames(x@annotation)[ncol(x@annotation)] = col;
	}
	x
} )
