#' @name export4GEDI
#' @aliases export4GEDI,StefansExpressionSet-method
#' @rdname export4GEDI-methods
#' @docType methods
#' @description  Convert the values in the StefansExpressionSet to the format GEDI program can import.
#' @param x the StefansExpressionSet object
#' @param fname the filename to export the data to
#' @param tag.col the sample name column in the samples table
#' @param time.col the time column in the samples table
#' @seealso \url{"http://apps.childrenshospital.org/clinical/research/ingber/GEDI/gedihome.htm"}
#' @title description of function export4GEDI
setGeneric('export4GEDI', ## Name
	function ( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) { ## Argumente der generischen Funktion
		standardGeneric('export4GEDI') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('export4GEDI', signature = c ( 'StefansExpressionSet') ,
	definition = function ( x, fname="GEDI_output.txt", tag.col = NULL, time.col=NULL, minSample_PerTime=1 ) {
	if ( is.null(time.col)){
		stop ( paste( "choose time.col from:", paste( colnames(x@samples), collapse=", ") ) ) 
	}
	if ( is.null(tag.col)){
		stop ( paste( "choose tag.col from:", paste( colnames(x@samples), collapse=", ") ) ) 
	}
	groupnames <- vector('numeric', nrow(x@samples))
	for (i in 1:nrow(x@samples)) {
		groupnames[i] = paste( x@samples[i, tag.col], x@samples[i, time.col] , sep="_" )
	}
	if ( length(which(table(groupnames) > 1)) > 0 ){
		stop ( "Sorry, please calculate the mean expression for the data first ('collaps()')") 
	}
	
	treatments <- unique( as.vector(x@samples[ , tag.col] ))
	
	## now I need to find all possible days
	possible <- NULL
	for ( t in treatments ){
		possible <- c( possible, x@samples[ grep( t, x@samples[ , tag.col]) , time.col] )
	}
	required = names(which(table(possible) > minSample_PerTime) )
	passed <- NULL
	for ( t in treatments ){
		l <- rep ( -1, nrow(x@samples) )
		p <- grep( t, x@samples[ , tag.col])
		if ( length(p) >= length(required) ){
			l[p] <- x@samples[ p, time.col]
			l[ is.na(match(l, required))== T ] <- -1
			passed <- c(passed, paste( c(paste("}",t,sep=""),  l), collapse="\t") )
		}
	}
	
	fileConn<-file(fname)
	open( fileConn, "w" )
	writeLines( paste( "}Dynamic", length(passed), nrow(x@data), x@name, sep="\t"), con=fileConn)
	writeLines( paste(passed , collapse="\n") , con=fileConn )
	
	for ( i in 1:nrow(x@data) ){
		writeLines( paste( c(rownames(x@data)[i], x@data[i,]), collapse="\t" ),  con=fileConn )
	}
	
	close(fileConn)
	print( paste ("created file", fname))
})
