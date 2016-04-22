#' @name drop.samples
#' @aliases drop.samples,StefansExpressionSet-method
#' @rdname drop.samples-methods
#' @docType methods
#' @description  drops samples from the StefansExpressionSet
#' @param x the StefansExpressionSet object
#' @param samplenames which samples to drop (samples like colnames(x@data))
#' @param name the name of the new StefansExpressionSet object
#' @title description of function drop.samples
#' @export 
setGeneric('drop.samples', ## Name
	function ( x, samplenames=NULL, name='dropped_samples' ) { ## Argumente der generischen Funktion
		standardGeneric('drop.samples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('drop.samples', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, samplenames=NULL, name='dropped_samples' ) {
	if ( ! is.null(samplenames)){
		red  <- new(class(x)[1], name=name )
		for (n in c( slot(x,'simple'), 'annotation') ){
			slot(red,n) <- slot(x,n)
		}
		red@samples <- x@samples[ is.na(match(x@samples[,x@sampleNamesCol], samplenames  ) ) == T ,]
		
		red@data <- x@data[, make.names(as.vector(red@samples[,red@sampleNamesCol]))]
		
		for ( n in c('colorRange') ) {
			red@usedObj[[n]] <- x@usedObj[[n]]
		}
	}
	red
})

