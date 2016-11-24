#' @name saveObj
#' @aliases saveObj,StefansExpressionSet-method
#' @rdname saveObj-methods
#' @docType methods
#' @description This function saves the object either as analysis.RData or norm_data.RData if the analysi.RData has not been produced before
<<<<<<< HEAD
#' @param obj the Rscexv object
#' @param file theoutfile default='SCExV_Grps.txt'
#' @title description of function analyse.data
=======
#' @param data the StefansExpressionSet object
#' @param file theoutfile default=paste(str_replace_all(data@name, "\\s+", "_" ), '.RData',sep='')
#' @title description of function saveObj
>>>>>>> branch 'master' of git@github.com:stela2502/StefansExpressionSet.git
#' @export 
setGeneric('saveObj', ## Name
		function ( data, file=NULL ){	## Argumente der generischen Funktion
			standardGeneric('saveObj') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('saveObj', signature = c ('StefansExpressionSet'),
		definition = function ( data, file=NULL ){
			if ( is.null(file)){
				file=paste(str_replace_all(data@name, "\\s+", "_" ), '.RData',sep='')
			}
			print ( paste('data exported to', file.path(data@outpath,file) ) )
			save(data , file=file.path(data@outpath, file) )
			
		}
)


