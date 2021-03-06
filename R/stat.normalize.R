#' @name normalize
#' @aliases normalize,StefansExpressionSet-method
#' @rdname normalize-methods
#' @docType methods
#' @description  normalize the expression data (sample wise)
#' This generic method is implemented in a object specific way - more options available
#' for each obejct specific function.
#' @param x The NGSexpressionSet
#' @param name the new name of the object (deafule old name + normalized)
#' @return the normalized data set (original data stored in NGS$raw
#' @title description of function normalize
#' @export 
setGeneric('normalize', ## Name
	function ( object , ..., name=NULL) { ## Argumente der generischen Funktion
		standardGeneric('normalize') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

#' @name normalize
#' @aliases normalize,NGSexpressionSet-method
#' @rdname normalize-methods
#' @docType methods
#' @description  normalize the expression data (sample wise)
#' @param x The NGSexpressionSet
#' @param readCounts The number of reads from each bam file or another value you want to normalize the data to
#' @param to_gene_length FALSE whether or not to normalize the data to gene length
#' @param geneLengthCol the column in the annotation data.frame to (in addition) normalize the genes to (e.g. trancript length)
#' @param name the new name of the object (deafule old name + normalized)
#' @return the normalized data set (original data stored in NGS$raw
#' @title description of function normalize
#' @export 
setMethod('normalize', signature = c ('NGSexpressionSet'),
		definition = function (  object, ..., readCounts=NULL, to_gene_length=FALSE, geneLengthCol='transcriptLength', name=NULL ) {
			if ( ! object@snorm ){
				browser()
				if ( is.null( readCounts ) ) {
					readCounts <- as.vector( DESeq2::estimateSizeFactorsForMatrix ( as.matrix(object@data)) )
				}
				object@samples$SizeFactor <- readCounts
				object@raw <- object@data
				object@data =  data.frame(t(apply(object@data,1, function(a) { a / readCounts } ) ))
				colnames(object@data) = colnames(object@raw)
				rownames(object@data) = rownames(object@raw)
				if (to_gene_length){
					for ( i in 1:nrow(object@data)) {
						object@data[i,] <- object@data[i,]/ (object@annotation[i,geneLengthCol ] / 1000)
					}
				}
				
				object@snorm=TRUE
				if(is.null(name)){
					name = paste( object@name ,'normalized' )
				}
				object@name = name
			}
			object
		})


#' @name normalize
#' @aliases normalize,SingleCellsNGS-method
#' @rdname normalize-methods
#' @docType methods
#' @description  
#' normalize the expression data by subsampling as described in PMID 24531970
#' @param x The SingleCellsNGS object
#' @param reads the required read depth
#' @param name the name of the new object
#' @return the normalized data set (original data stored in slot 'raw'
#' @title description of function normalize
#' @export 
setMethod('normalize', signature = c ('SingleCellsNGS'),
		definition = function (  object, ..., reads=600, name='normalized' ) {
			if ( ! object@snorm ) {
			if ( length( object@samples$counts ) == 0 ) {
				object@samples$counts <- apply( object@data, 2, sum)
			}
			object <- drop.samples( object, object@samples[which(object@samples$counts < reads), object@sampleNamesCol ] 
					, name=name )
			
			if ( ! object@snorm ){
				object@raw <- object@data
			}
			## resample the data
			n <- nrow(object@raw)
			object@data[] <- 0
			for ( i in 1:ncol(object@raw) ) {
				d <- sample(rep ( 1:n, object@raw[,i]) , reads)
				t <- table(d)
				object@data[ as.numeric(names(t)),i] <- as.numeric(t)
			}
			}
			object@snorm = TRUE
			object
		}
)




#' @name normalize
#' @aliases normalize,StefansExpressionSet-method
#' @rdname normalize-methods
#' @docType methods
#' @description  constructor that has to be implemented for a generic StefansExpressionSet
#' This generic version was meant for array data and I have not had the need nor time to implement this part.
#' @param x the StefansExpressionSet object
#' @title description of function normalize
setMethod('normalize', signature = c ('StefansExpressionSet') ,
	definition = function ( object , ..., name=NULL) {
	x
})

