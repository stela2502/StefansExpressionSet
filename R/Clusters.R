#' @name addCluster
#' @aliases addCluster,Clusters-method
#' @rdname addCluster-methods
#' @docType methods
#' @description 
#' @param x  the StefansExpressionSet object
#' @param userGroups  the userGroups table 
#' @param col  the column in the userGroups table that contains the group ids default='groupID'
#' @param what one of 'cells' or 'genes'. Where to add the group.  default= 'cells'
#' @title description of function addCluster
#' @export 
setGeneric('addCluster', ## Name
	function (x, userGroups, col='groupID', what= 'cells' ) { ## Argumente der generischen Funktion
		standardGeneric('addCluster') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('addCluster', signature = c ('Clusters'),
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
#' @description 
#' @param x  TEXT MISSING
#' @param clusters  TEXT MISSING
#' @param groups.n  TEXT MISSING
#' @title description of function difference
#' @export 
setGeneric('difference', ## Name
	function ( x, clusters, groups.n ) { ## Argumente der generischen Funktion
		standardGeneric('difference') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('difference', signature = c ('Clusters'),
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
#' @description 
#' @param x  TEXT MISSING
#' @param what  TEXT MISSING default='cells'
#' @param col  TEXT MISSING default='groupID'
#' @title description of function quality_of_fit
#' @export 
setGeneric('quality_of_fit', ## Name
	function ( x, what='cells', col='groupID' ) { ## Argumente der generischen Funktion
		standardGeneric('quality_of_fit') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('quality_of_fit', signature = c ('Clusters'),
	definition = function ( x, what='cells', col='groupID' ) {
	
	test <- x$data
	test[which(test ==  0 ) ] = NA
	if ( what=='cells') {
		clusters <- x$samples[,col]
		ret <- list ( 'single' = apply(test,2, difference, clusters, max(clusters) ) )
		ret$sum = round(sum(ret$single))
	}
	else if ( what=='genes') {
		clusters <- x$annotation[,col]
		ret <- list ( 'single' = apply(test,1, difference, clusters, max(clusters) ) )
		ret$sum = round(sum(ret$single))
	}
	else {
		stop(paste( what,': what option is not supported!') )
	}
	ret
} )
