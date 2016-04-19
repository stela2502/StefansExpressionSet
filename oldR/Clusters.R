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
#' @name difference
#' @aliases difference,Clusters-method
#' @rdname difference-methods
#' @docType methods
#' @description this code has to be revised! NOT WORKING
#' @param x the StefansExpressionSet object
#' @param clusters the grouoing strings
#' @param groups.n the max number of groups
#' @title description of function difference
setGeneric('difference', ## Name
	function ( x, clusters, groups.n ) { ## Argumente der generischen Funktion
		standardGeneric('difference') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('difference', signature = c ('data.frame'),
	definition = function ( x, clusters, groups.n ) {
	ret = 0 
	for ( i in 1:groups.n  ) {
		a <- x[which( clusters == i)]
		a <- a[- (is.na(a))==F]
		if ( length(a) > 1 ) {  ret = ret + sum( (a- mean(a) )^2 ) }
	}
	ret
} )
#' @name quality_of_fit
#' @aliases quality_of_fit,Clusters-method
#' @rdname quality_of_fit-methods
#' @docType methods
#' @description this calculated a quality of fit for the expression data
#' @param x the StefansExpressionSet object
#' @param what where to group the data on (cells or genes; default='cells')
#' @param col which samples/annoation column to use for the grouping (default='groupID')
#' @title description of function quality_of_fit
#' @export 
setGeneric('quality_of_fit', ## Name
	function ( x, what='cells', col='groupID' ) { ## Argumente der generischen Funktion
		standardGeneric('quality_of_fit') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('quality_of_fit', signature = c ('StefansExpressionSet'),
	definition = function ( x, what='cells', col='groupID' ) {
	
	test <- x@data
	test[which(test ==  0 ) ] = NA
	if ( what=='cells') {
		clusters <- x@samples[,col]
		ret <- list ( 'single' = apply(test,2, difference, as.character(clusters), max(clusters) ) )
		ret$sum = round(sum(ret$single))
	}
	else if ( what=='genes') {
		clusters <- x@annotation[,col]
		ret <- list ( 'single' = apply(test,1, difference, as.charcter(clusters), max(clusters) ) )
		ret$sum = round(sum(ret$single))
	}
	else {
		stop(paste( what,': what option is not supported!') )
	}
	ret
} )
