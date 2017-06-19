#' @name as_BioData
#' @rdname as_BioData-methods
#' @docType methods
#' @description create a BioData object from a StefansExpressionSet object
#' @param dat the (old) StefansExpressionSet object
#' @title description of function as
#' @export 
setGeneric('as_BioData', ## Name
		function ( dat ) { ## Argumente der generischen Funktion
			standardGeneric('as_BioData') ## der Aufruf von standardGeneric sorgt für das_BioData Dispatching
		}
)

setMethod('as_BioData', signature = c ('StefansExpressionSet'),
		definition = function ( dat ) {
			#cbind(annotation,dat), Samples=samples, name="testObject",namecol='sname', outpath = ""
			ret <- BioData$new( cbind( dat@annotation, dat@data), Samples=dat@samples, name= dat@name, namecol=dat@sampleNamesCol, namerow=dat@rownamescol, outpath=dat@outpath )
			ret$usedObj <- dat@usedObj
			ret
		} )
