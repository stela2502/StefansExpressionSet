#' @name reduce.Obj
#' @aliases reduce.Obj,StefansExpressionSet-method
#' @rdname reduce.Obj-methods
#' @docType methods
#' @description  reduces the dataset based on genes e.g. dropps genes from the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @param probeSets a list of probesets to reduce the data to
#' @param name the new StefansExpressionSet name
#' @title description of function reduce.Obj
#' @export 
setGeneric('reduce.Obj', ## Name
	function ( x, probeSets=c(), name="reducedSet" ) { ## Argumente der generischen Funktion
		standardGeneric('reduce.Obj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('reduce.Obj', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, probeSets=c(), name="reducedSet" ) {
	retObj <- new(class(x)[1], name = name)
	useOnly <- match(tolower(probeSets), tolower(rownames(x@data)))
	not.matched <- probeSets[is.na(useOnly)]
	if ( length(not.matched) > 0 ){
		print (paste('Problematic genes:', paste(not.matched,sep=', ')))
		probeSets <- probeSets[ ! is.na(useOnly)]
		useOnly <- useOnly[ ! is.na(useOnly) ]
	}
	for (n in slot(x,'simple')){
		slot(retObj,n) <- slot(x,n)
	}
	retObj@samples <- x@samples
	retObj@data <- data.frame( x@data[ useOnly ,] )
	rownames(retObj@data) <- probeSets
	colnames(retObj@data) <- colnames(x@data)
	retObj@annotation <- data.frame(x@annotation[useOnly,]) ## if I only have one column here
	colnames(retObj@annotation) <- colnames(x@annotation)
	if ( length( names(x@stats)) > 0){
		for ( i in 1:length(names(x@stats))){
			retObj@stats[[i]]= x@stats[[i]][ match(probeSets ,x@stats[[i]][,1] ),]
		}
		names(retObj@stats) <- names(x@stats)
	}
	if ( length(x@ranks) > 0 ){
		retObj@ranks = x@ranks[useOnly,]
	}
	if ( ncol( x@raw ) > 0 ) {
		retObj@raw = x@raw[useOnly,]
	}
	retObj
})



