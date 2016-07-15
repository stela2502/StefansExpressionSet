#' @name restrictSamples
#' @aliases restrictSamples,StefansExpressionSet-method
#' @rdname restrictSamples-methods
#' @docType methods
#' @description Drop the samples, that have been selected!
#' @param x the StefansExpressionSet object
#' @param name the name of the new StefansExpressionSet
#' @param column which column to analyze in the samples table
#' @param value which value to take as 'cut off'
#' @param mode one of 'less', 'more', 'onlyless', 'equals' or 'grep'
#' @title description of function restrictSamples
#' @export 
setGeneric('restrictSamples', ## Name
	function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ) { ## Argumente der generischen Funktion
		standardGeneric('restrictSamples') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('restrictSamples', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, column='Analysis', value=NULL, name='newSet', mode= 'equals' ) {
		
	S <- NULL
	if ( is.null(value)) {
		stop( "the value must not be NULL!")
	}
	
	switch( mode,
			'less' = S <- x@samples[which ( as.numeric(x@samples[,column]) >  value), ], 
			'more' = S <- x@samples[which ( as.numeric(x@samples[,column]) < value ), ], 
			'onlyless' = S <- x@samples[which ( as.numeric(x@samples[,column])  >= value ), ],
			'equals' = S <- x@samples[which ( x@samples[,column] ==  value), ],
			'grep' = S <- x@samples[ grep(value,x@samples[,column]) ,]
	)
	if ( nrow(S) < nrow(x@samples)){
		x <- drop.samples( x, S[,x@sampleNamesCol], name=name)
	}
	write.table (S, file=paste(name,'_Sample_Description', ".xls", sep=''), sep='\t',  row.names=F,quote=F )
	x
})

####last

