#' @name StefansExpressionSet
#' @aliases StefansExpressionSet,data.frame-method
#' @rdname StefansExpressionSet-methods
#' @docType methods
#' @description  this file contains all generic fnction for data export and ploting Create an StefansExpressionSet
#' @description  object (S3) This object is mainly used for subsetting of the data and plotting @export
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID
#' @param outpath Where to store the output from the analysis
#' @title description of function StefansExpressionSet
#' @exportMethod  StefansExpressionSet
setGeneric("StefansExpressionSet", ## Name
		function( dat, Samples, class='StefansExpressionSet', name='WorkingSet', namecol=NULL, namerow= 'GeneID', outpath = ''){ ## Argumente der generischen Funktion
			standardGeneric("StefansExpressionSet") ## der Aufruf von standardGeneric sorgt für das Dispatching
		})

setMethod("StefansExpressionSet", signature = c ('data.frame'),
		definition = function ( dat, Samples, class='StefansExpressionSet', name='filename', namecol=NULL, namerow= 'GeneID', outpath = '' ) {
			S <- Samples
			
			if ( is.null(namecol)){
				r <-  apply ( S,2, function( x, coln ) { length( which( is.na(match(x, coln))==F))  }, colnames( dat) )
				namecol = names( which( r == max(r) ) )
			}
			
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
			
			if ( class(annotation) == 'factor'){
				annotation <- data.frame( annotation )
				colnames(annotation) <- namerow
			}
			if ( class(annotation) == 'character'){
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
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID
#' @param outpath Where to store the output from the analysis
#' @return A new NGSexpressionSet object
#' @title description of function NGSexpressionSet
#' @export
setGeneric('NGSexpressionSet', ## Name
		function ( dat, ... ) { ## Argumente der generischen Funktion
			standardGeneric('NGSexpressionSet') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('NGSexpressionSet', signature = c ('data.frame'),
		definition = function ( dat, Samples,  name='WorkingSet', namecol='GroupName', namerow= 'GeneID', outpath = NULL) {
			x <- StefansExpressionSet( dat, Samples,  name=name, namecol=namecol, namerow= namerow, outpath =  outpath)
			x <- as(x,'NGSexpressionSet')
			x
		}
)

setMethod('NGSexpressionSet', signature = c ('list'),
		definition = function ( dat ) {
			ret = NULL
			if (all.equal( names ( dat), c("counts" ,"annotation", "targets", "stat")  ) ) {
				samples <- data.frame(t(dat$stat))
				colnames(samples) <- as.NGSexpressionSet.character(t(samples[1,]))
				samples$filename <- rownames(samples)
				rownames(samples) <-1:nrow(samples)
				ret <- NGSexpressionSet( 
						dat= cbind(dat$annotation, dat$counts), 
						samples = samples, 
						namecol= 'filename', 
						namerow= 'GeneID',
						outpath= ''
				)
			}
			else {
				print ("The list needs to contain the entries counts ,annotation, targets and stat" )
			}
			ret
		} )

#' @name SingleCellsNGS
#' @aliases SingleCellsNGS,SingleCellsNGS-method
#' @rdname SingleCellsNGS-methods
#' @docType methods
#' @description  create a new SingleCellsNGS object This object is mainly used for plotting the
#' @description  data
#' @param dat data frame or matrix containing all expression data
#' @param Samples A sample description table
#' @param name The name of this object is going to be used in the output of all plots and tables - make it specific
#' @param namecol The samples table column that contains the (to be set) names of the samples in the data matrix
#' @param namerow This is the name of the gene level annotation column in the dat file to use as ID
#' @param outpath Where to store the output from the analysis
#' @return A new SingleCellsNGS object
#' @title description of function SingleCellsNGS
#' @export
setGeneric('SingleCellsNGS', ## Name
		function ( dat, Samples,  name='WorkingSet', namecol='GroupName', namerow= 'GeneID',  outpath = '') { ## Argumente der generischen Funktion
			standardGeneric('SingleCellsNGS') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('SingleCellsNGS', signature = c ('data.frame'),
		definition = function ( dat, Samples,  name='WorkingSet', namecol='GroupName', namerow= 'GeneID', outpath = '') {
			x <- StefansExpressionSet( dat, Samples, name=name, namecol=namecol, namerow= namerow, outpath =  outpath)
			as(x,'SingleCellsNGS')
		} )


#' @name FromCountsObj
#' @aliases FromCountsObj,FromCountsObj-method
#' @rdname FromCountsObj-methods
#' @docType methods
#' @description  create a new StefansExpressionSet object based on a DEseq counts object
#' @param dat the counts object
#' @param type the classnam 'StefsnExpressionSet', 'NGSexpressionSet' (default) or 'SingleCellNGS'
#' @param outpath the outpath for this object
#' @param name the name of the new obejct
#' @title description of function FromCountsObj
#' @export 
setGeneric('FromCountsObj' ,
	function ( dat, type="NGSexpressionSet", outpath="", name="WorkingSet" ) {
		standardGeneric('FromCountsObj')
	} )

setMethod('FromCountsObj', signature = c ('list'),
		definition = function ( dat, type="NGSexpressionSet", outpath=NULL, name='WorkingSet') {
			if ( is.null(outpath) ){
				outpath = paste(pwd(),"../output/", sep='')
			}
			if ( ! all.equal( sort(names(dat)) ,c("annotation","counts","stat","targets") ) ) {
				stop ( "dat is not a counts object containing the right entries" )
			}
			if ( is.na(match(type, c('StefansExpressionSet', 'NGSexpressionSet', 'SingleCellsNGS')) ) ){
				stop( "the type has to be one of c('StefansExpressionSet', 'NGSexpressionSet', 'SingleCellsNGS')" )
			}
						
			samples <- data.frame(t(dat$stat))
			colnames(samples) <- as.character(t(samples[1,]))
			samples <- samples[-1,]
			samples$filename <- rownames(samples)
			
			dat$annotation[,2] <- unlist(lapply(dat$annotation[,2],function(x) { t <- str_split(x,';' ); t[[1]][1] } ))
			dat$annotation[,3] <- unlist(lapply(dat$annotation[,3],function(x) { t <- unlist(str_split(x,';' )); min(t) } ))
			dat$annotation[,4] <- unlist(lapply(dat$annotation[,4],function(x) { t <- unlist(str_split(x,';' )); max(t) } ))
			
			x <- StefansExpressionSet(  cbind(dat$annotation, dat$counts),
					samples ,
					namecol='filename',
					namerow=colnames(dat$annotation)[1],
					outpath=outpath,
					name=name
			)			
			as(x, type )
		} )

