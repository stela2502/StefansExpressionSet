#' @name mds
#' @aliases mds,StefansExpressionSet-method
#' @rdname mds-methods
#' @docType methods
#' @description Calculates the MDS for a given MDS type and stores the 3 dimensions in the object for later use.
#' @param dataObj the StefansExpressionSet object
#' @param mds.type Which MDS function should be called default="PCA"
#' @param onwhat condense which dataset at the moment only Expression is supported default='expression'
#' @param genes do it on genes not on samples (default = F)
#' @title description of function mds.and.clus
#' @export 
setGeneric('mds', ## Name
	function ( dataObj, ..., mds.type="PCA" , onwhat ='Expression', genes=F ) { ## Argumente der generischen Funktion
		standardGeneric('mds') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('mds', signature = c ('StefansExpressionSet'),
	definition = function ( dataObj, ..., mds.type="PCA", onwhat ='Expression', genes=F ) {
	## the code is crap re-code!!
	if(onwhat=="Expression"){
		tab <- as.matrix(t(dataObj@data))
	} 
	else {
		stop( paste("Sorry, the option onwhat",onwhat,"is not supported") )
	}
	if ( genes ) {
		tab <- t(tab)
	}
	this.k <- paste(onwhat,mds.type)
	if ( is.null( dataObj@usedObj$MDS) ) {
		dataObj@usedObj$MDS = list()
	}
	if ( (is.null(dataObj@usedObj$MDS[[this.k]])) ||  all.equal( rownames(dataObj@usedObj$MDS[[this.k]]), colnames(dataObj@data) )==F ) {
		mds.proj <- NULL
		pr <- NULL
		#system ( 'rm loadings.png' )
		if(mds.type == "PCA"){
			pr <- prcomp(tab)
			mds.proj <- pr$x[,1:3]
			png ( file=file.path( dataObj@outpath,'loadings.png'), width=1000, height=1000 )
			plot (  pr$rotation[,1:2] , col='white' );
			text( pr$rotation[,1:2], labels= rownames(pr$rotation), cex=1.5 )
			dev.off()
			write.table( cbind( Genes = rownames(pr$rotation), pr$rotation[,1:2] ), 
					file=file.path( dataObj@outpath,'gene_loadings.xls') , row.names=F, sep='\t',quote=F )
			#	mds.trans <- prcomp(t(tab))$x[,1:3]
		} else if ( mds.type=='DM') {
			if (!library("destiny", quietly = TRUE,logical.return=TRUE )) {
				stop("package 'destiny' needed for this function to work. Please install it.",
						call. = FALSE)
			}
			if ( ! exists('sigma', mode='numeric') ){
				if ( exists('find_sigmas', where='package:destiny', mode='function') ){
					sigmas <- destiny::find_sigmas(tab, verbose=F)
					sigma <- destiny::optimal_sigma(sigmas)
				}else {
					sigmas <- destiny::find.sigmas(tab, verbose=F)
					sigma <- destiny::optimal.sigma(sigmas)
				}
				
			}
			if ( !exists('distance', mode='character')){
				distance = 'cosine'
			}
			dm <- destiny::DiffusionMap(tab, distance = distance, sigma = sigma)
			mds.proj <- destiny::as.data.frame(dm)[,1:3]
		}
		else if ( mds.type == "LLE"){
			mds.proj <- LLE( tab, dim = 3, k = as.numeric(LLEK) )
			#	mds.trans <- LLE( t(tab), dim = 3, k = as.numeric(LLEK) )
		}else if ( mds.type == "ISOMAP"){
			mds.proj <- Isomap( tab, dim = 3, k = as.numeric(LLEK) )$dim3
			#	mds.trans <- Isomap( t(tab), dim = 3, k = as.numeric(LLEK) )$dim3
		}else if ( mds.type == "ZIFA" ) {
			stop( "Sorry ZIFA has to be double checked - are you working on normalized data - than ZIFA can not be applied!")
			print ( "Running external python script to apply ZIFA dimensional reduction (PCR data only)" )
			if ( genes ) {
				ZIFA <- 40.000001 - dataObj@raw
			}
			else {
				ZIFA <- 40.000001 - t(dataObj@raw)
			}
			
			write.table( ZIFA, file="ZIFA_input.dat", sep=" ", col.names=F, row.names=F , quote=F)
			write( c("from ZIFA import ZIFA","from ZIFA import block_ZIFA", "import numpy as np",
							"Y = np.loadtxt('ZIFA_input.dat')", "Z, model_params = ZIFA.fitModel( Y, 3 )", 
							"np.savetxt('TheMDS_ZIFA.xls', Z )" ), 
					file= 'ZIFA_calc.py' )
			system( "python ZIFA_calc.py" )
			Sys.sleep(5)
			mds.proj <- read.delim( "TheMDS_ZIFA.xls", sep=' ', header=F)
			if ( genes ) {
				rownames(mds.proj) <- colnames(ZIFA)
			}else {
				rownames(mds.proj) <- rownames(ZIFA)
			}
			colnames(mds.proj) <- c( 'x','y','z')
			
		} else if ( mds.type == "DDRTree" ) {
			tab <- t(tab) ## mds works the other way round!
			DDRTree_res <- DDRTree( tab, dimensions=3)
			mds.proj <- t(DDRTree_res$Z)
			rownames(mds.proj) <- colnames(tab)
			
			if ( genes ) {	
				dataObj@usedObj$DRRTree.genes <- DDRTree_res
			}else {
				dataObj@usedObj$DRRTree <- DDRTree_res
			}
		}
		else {
			print( paste("Sory I can not work on the mds.type option",mds.type) )
		}
		if ( genes ) {
			if ( is.null(dataObj@usedObj$MDSgenes)){
				dataObj@usedObj$MDSgenes <- list()
			}
			dataObj@usedObj$MDSgenes[[mds.type]] <- mds.proj
		}else{
			if ( is.null(dataObj@usedObj$MDS)){
				dataObj@usedObj$MDS <- list()
			}
			dataObj@usedObj$MDS[[mds.type]] <- mds.proj
		}
		

	}
#	dataObj <- clusters ( dataObj, onwhat=onwhat, clusterby=clusterby, groups.n = groups.n,
#			ctype = ctype, cmethod=cmethod, useGrouping=useGrouping )
	
	dataObj
} )
