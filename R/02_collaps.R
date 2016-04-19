#' @name collaps
#' @aliases collaps,StefansExpressionSet-method
#' @rdname collaps-methods
#' @docType methods
#' @description  This function will collpase the data in the StefansExpressionSet to only contain one value
#' @description  per sample group.
#' @param dataObj the StefansExpressionSet
#' @param by collapsing method c('median','mean','sd','sum', or own function )
#' @param groupCol the sample names you want to group on
#' @title description of function collaps
#' @export 
setGeneric('collaps', ## Name
		function (dataObj, by=c('median','mean','sd','sum' ), groupCol='GroupID' ) { ## Argumente der generischen Funktion
			standardGeneric('collaps') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('collaps', signature = c ('StefansExpressionSet'),
		definition = function (dataObj, by=c('median','mean','sd','sum', 'var' ), groupCol='GroupID' ) {
			u <- unique(as.vector(dataObj@samples[,groupCol]))
			m <- length(u)
			mm <-  matrix ( rep(0,m * nrow(dataObj@data)), ncol=m)
			colnames(mm) <- u
			rownames(mm) <- rownames(dataObj@data)
			f <- NULL
			if ( is.function(by)){
				f <- by
			}else {
			switch( by,
					median = f<- function (x ) { median(x) },
					mean = f <- function(x) { mean(x) },
					sd = f <- function(x) { sd(x) },
					sum = f <-function(x) { sum(x)},
					var = f <- function(x) { var(x) }
			);
			}
			if ( is.null(f) ) {
				stop("Please set what to one of 'median','mean','sd','sum'" )
			}
			new_samples <- NULL
			for ( i in u ){
				all <- which(as.vector(dataObj@samples[, groupCol]) == i )
				new_samples <- rbind ( new_samples, dataObj@samples[all[1],] )
				mm[,i] <- apply( dataObj@data[ , all],1,f)
			}
			name = paste( unlist(strsplit( paste( dataObj@name, groupCol, by, sep='_') , '\\s')) , collapse='_')
			try ( { ret <- drop.samples ( dataObj,
					samplenames= setdiff( 
							dataObj@samples[,dataObj@sampleNamesCol] , 
							new_samples[,dataObj@sampleNamesCol] 
					), 
					name= name
			) }, silent=TRUE )
			colnames(mm) <- as.vector(new_samples[, dataObj@sampleNamesCol])
			ret@data <- data.frame(mm)
			ret
})




