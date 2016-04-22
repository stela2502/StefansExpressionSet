#' @name StefansExpressionSet
#' @aliases StefansExpressionSet,data.frame-method
#' @rdname StefansExpressionSet-methods
#' @docType methods
#' @description  this file contains all generic fnction for data export and ploting Create an StefansExpressionSet
#' @description  object (S3) This object is mainly used for subsetting of the data and plotting @export
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
#' @param outpath Where to store the output from the analysis
#' @param annotation The annotation table from e.g. affymetrix csv data
#' @param newOrder The samples column name for the new order (default 'Order')
#' @title description of function StefansExpressionSet
#' @exportMethod  StefansExpressionSet
setGeneric("StefansExpressionSet", ## Name
		function( dat, Samples, class='StefansExpressionSet',  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = ''){ ## Argumente der generischen Funktion
			standardGeneric("StefansExpressionSet") ## der Aufruf von standardGeneric sorgt für das Dispatching
		})

setMethod("StefansExpressionSet", signature = c ('data.frame'), 
		definition = function ( dat, Samples, class='StefansExpressionSet',  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = '' ) {
			S <- Samples
			if ( ! is.null(Analysis) ){
				S <- Samples[which ( Samples$Analysis == Analysis ), ]
			}
			if ( ! is.null(usecol) ) {
				S <- Samples[which ( Samples[, usecol] == 1 ),]
			}
			if ( exists('filename',S) ) {
				n <- make.names(as.vector(S$filename))
				mat <- match( as.vector(S$filename), colnames(dat))
				if ( sum(is.na(mat)) > 0 ) {
					stop(paste( 'The files', 
						paste( as.vector(S$filename)[is.na(mat)], collapse=', '),
						'Do not have data column in the "dat" data.frame' ) 
					)
				}
				ret <- dat[, mat ]
				annotation <- dat[, is.na(match( colnames(dat), as.vector(S$filename) ))==T ]
			}else{
				n <- make.names(as.vector(S[,namecol]))
				mat <- match( as.vector(S[,namecol]), colnames(dat))
				if ( sum(is.na(mat)) > 0 ) {
					stop(paste( 'The samples', 
						paste( as.vector(S[,namecol])[is.na(mat)], collapse=', '),
						'Do not have data column in the "dat" data.frame' ) 
					)
				}
				ret <- dat[, mat ]
				annotation <- dat[, is.na(match( colnames(dat), as.vector(S[,namecol]) ))==T ]
			}
			if ( class(annotation) == 'factor'){
				annotation <- data.frame( annotation )
				colnames(annotation) <- namerow
			}
			
			if ( exists( 'Order', S)){
				ret <- ret[, order(S$Order)  ]
				S <- S[order(S$Order), ]
			}
			
			if ( outpath == '' ){
				outpath = pwd()
			}
			if ( ! file.exists(outpath)){
				dir.create( outpath )
			}
			ret <- data.frame(ret)
			if ( is.null(dim(annotation))){
				## this xcould be a problem... hope we at least have a vector
				if ( length(annotation) == nrow(ret)) {
					rownames(ret) <- annotation
				}
				else {
					stop ( "Sorry, please cbind the rownames to the expression values before creating this object!")
				}
			}else{
				rownames(ret) <- annotation[,namerow]
			}
			
			data <- list ( 'data' = ret, samples = S, name= name, annotation = annotation, rownamescol = namerow )
			data$outpath <- outpath
			
			colnames(ret) <- make.names(forceAbsoluteUniqueSample ( as.vector(S[, namecol]) ))
			S$SampleName <- colnames(ret)
			
			write.table (cbind(rownames(ret), ret ), file=paste(name, '_DataValues',".xls", sep=''), sep='\t',  row.names=F,quote=F )
			write.table (S, file=paste(name,'_Sample_Description', ".xls", sep=''), sep='\t',  row.names=F,quote=F )
			data$sampleNamesCol <- namecol
			data$batchRemoved=0
			r <- new ( class, data = data$data, samples = data$samples, name = data$name, annotation = data$annotation, rownamescol= data$rownamescol,sampleNamesCol = data$sampleNamesCol , outpath= data$outpath )
			force.numeric(r)
		}
)



#' @name NGSexpressionSet
#' @aliases NGSexpressionSet,NGSexpressionSet-method
#' @rdname NGSexpressionSet-methods
#' @docType methods
#' @description  create a new NGSexpressionSet object This object is mainly used for plotting the
#' @description  data
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
#' @param outpath Where to store the output from the analysis
#' @param annotation The annotation table from e.g. affymetrix csv data
#' @param newOrder The samples column name for the new order (default 'Order')
#' @return A new NGSexpressionSet object
#' @title description of function NGSexpressionSet
#' @export 
setGeneric('NGSexpressionSet', ## Name
		function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) { ## Argumente der generischen Funktion
			standardGeneric('NGSexpressionSet') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('NGSexpressionSet', signature = c ('data.frame'),
		definition = function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = NULL) {
			x <- StefansExpressionSet( dat, Samples,  Analysis = Analysis, name=name, namecol=namecol, namerow= namerow, usecol= usecol, outpath =  outpath)
			x <- as(x,'NGSexpressionSet')
			x
		} 
)


#' @name SingleCellsNGS
#' @aliases SingleCellsNGS,SingleCellsNGS-method
#' @rdname SingleCellsNGS-methods
#' @docType methods
#' @description  create a new SingleCellsNGS object This object is mainly used for plotting the
#' @description  data
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param Analysis If the samples table contains a Analysis column you can subset the data already here to the Analysis here (String). This table has to contain a column filenames that is expected to connect the sample table to the dat column names.
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID 
#' @param outpath Where to store the output from the analysis
#' @param annotation The annotation table from e.g. affymetrix csv data
#' @param newOrder The samples column name for the new order (default 'Order')
#' @return A new SingleCellsNGS object
#' @title description of function SingleCellsNGS
#' @export 
setGeneric('SingleCellsNGS', ## Name
		function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = '') { ## Argumente der generischen Funktion
			standardGeneric('SingleCellsNGS') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('SingleCellsNGS', signature = c ('data.frame'),
		definition = function ( dat, Samples,  Analysis = NULL, name='WorkingSet', namecol='GroupName', namerow= 'GeneID', usecol='Use' , outpath = '') {
			x <- StefansExpressionSet( dat, Samples,  Analysis = Analysis, name=name, namecol=namecol, namerow= namerow, usecol= usecol, outpath =  outpath)
			as(x,'SingleCellsNGS')
		} )


