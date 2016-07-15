#' @name add2usedObj
#' @aliases add2usedObj,add2usedObj-method
#' @rdname add2usedObj-methods
#' @docType methods
#' @description Add values to the usedObj list including some checks
#' @param dataObj the StefansExpressionSet object
#' @param list.name the toplevel name in the usedObj slot to store the data in
#' @param data.name the second level data name (the first level becomes a list; default = NULL)
#' @param object the object to store in the usedObj
#' @title description of function add2usedObj
#' @export 
setGeneric('add2usedObj', ## Name
	function ( dataObj,list.name, data.name, object ) { ## Argumente der generischen Funktion
		standardGeneric('add2usedObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('add2usedObj', signature = c('StefansExpressionSet') ,
	definition = function ( dataObj,list.name, data.name=NULL, object ) {
	if ( is.null(dataObj@usedObj[[list.name]] ) ) {
		dataObj@usedObj[[list.name]] <- list()
	}
	if ( is.null(data.name)){
		dataObj@usedObj[[list.name]] <- object
	}else {
		
	}
	dataObj
})
		


