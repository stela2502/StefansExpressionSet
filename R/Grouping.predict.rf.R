#' @name predict.rf
#' @aliases 'predict.rf,SingleCellsNGS-method
#' @docType methods
#' @description simple prediction of groups using a random forest trained during the bestGrouping process
#' @param x the single cells ngs object
#' @param rf the random forst model to use for the classification
#' @param bestColname the column name to store the results in
#' @title description of function predict.rf
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('predict.rf',
		function ( x, rf,  bestColname='predicted group using random forest'){
			standardGeneric('predict.rf')
		}
)
setMethod('predict.rf', signature = c ('SingleCellsNGS'),
		definition = function (x, rf, bestColname='predicted group using random forest') {
			predicted2 <-predict( rf, t(as.matrix(workingSet@data)) )
			
		}
)


