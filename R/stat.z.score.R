#' @name z.score
#' @aliases z.score,StefansExpressionSet-method
#' @rdname z.score-methods
#' @docType methods
#' @description  z score the matrix
#' @param m the matrix of column = samples and rows = genes or an StefansExpressionSet
#' @return the z scored matrix
#' @title description of function z.score
#' @export 
setGeneric('z.score', ## Name
		function (m) { ## Argumente der generischen Funktion
			standardGeneric('z.score') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)


setMethod('z.score', signature = c ('SingleCellsNGS'),
		definition = function ( m ) {
			
			if ( ! m@zscored ) {
				m@raw <- m@data
				ma  <- as.matrix(m@data)
				i = 0
				ret <- t(
						apply(ma,1, function (x) {
									i = i+1
									n <- which(x==0)
									if ( length(x) - length(n) > 1 ){
										if (length(n) == 0 ){
											x <-  scale(as.vector(t(x)))
										}
										else {
											x[-n] <- scale(as.vector(t(x[-n])))
											x[n] <- -20
										}
										
									}
									else {
										x[] = -20
									}
									x}
						)
				)
				m@data <- data.frame(ret)
				colnames(m@data)<- colnames(m@raw)
				m@zscored = TRUE
			}
			m
		})

setMethod('z.score', signature = c ('matrix'),
		definition = function (m ) {
			rn <- rownames( m )
			me <- apply( m, 1, mean )
			sd <- apply( m, 1, sd )
			sd[which(sd==0)] <- 1e-8
			m <- (m - me) /sd
			rownames(m) <- rn
			m
		})

setMethod('z.score',signature = c ('StefansExpressionSet'),
		definition = function (m) {
			if (! m@zscored ){
			m@data <- data.frame(z.score( as.matrix( m@data )))
			m@zscored = TRUE
			}
			m
		})

