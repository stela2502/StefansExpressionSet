#' @name rfCluster_col
#' @aliases 'rfCluster_col,StefansExpressionSet-method
#' @title rfCluster_col
#' @name rfCluster_col-methods
#' @docType methods
#' @description This fucntion uses the RFclust.SGE to create fandomForest based unsupervised clusters on a subset of the data.
#' @description Default is on 200 cells using all (provided) genes with 500 forests and 500 trees per forest for 5 repetitions.
#' @description You are asked to give a k numer of expected clusters (better too many than too little), classifies the total 
#' @description data using the 5 different unsupervised runs and all cluster ids from these runs are merged into the final cluster id.
#' @description This <summaryCol> will be part of the return objects samples table, together with a <usefulCol> where
#' @description all clusters with less than 10 cells have been merged into the 'gr. 0'.
#' @description The final results will be reported as new columns in the samples table containing the 'name'
#' @param x the single cells ngs object
#' @param email your email to use together with the SGE option
#' @param SGE whether to use the sun grid engine to calculate the rf grouping
#' @param rep how many repetitions for the random forest grouping should be run (default = 5)
#' @param slice how many processes should be started for each random forest clustering (default = 30)
#' @param bestColname the column name to store the results in
#' @param k the numer of expected clusters (metter more than to view)
#' @param subset how many cells should be randomly selected for the unsupervised clustering (default = 200)
#' @param name if you want to run multiple RFclusterings on e.g. using different input genes you need to specify a name (default ='RFclust')
#' @param nforest the numer of forests to grow for each rep (defualt = 500)
#' @param ntree the numer of trees per forest (default = 500)
#' @param settings a list of slurm settings (A t and p are used) that will make the RFclust.SGE module use a slurm backend for calculations.
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('rfCluster_col',
		function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200,nforest=500, ntree=500, name='RFclust', settings=list()){
			standardGeneric('rfCluster_col')
		}
)
setMethod('rfCluster_col', signature = c ('StefansExpressionSet'),
		definition = function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200 ,nforest=500, ntree=1000, name='RFclust', settings=list()) {
			summaryCol=paste( 'All_groups', name,sep='_')
			usefulCol=paste ('Usefull_groups',name, sep='_')
			n= paste(x@name, name,sep='_')
			m <- max(k)
			OPATH <- paste( x@outpath,"/",str_replace( x@name, '\\s', '_'), sep='')
			opath = paste( OPATH,"/RFclust.mp",sep='' )
			
			if ( ! dir.exists(OPATH)){
				dir.create( OPATH )
			}
			processed = FALSE
			single_res_col <- paste('RFgrouping',name)
			for ( i in 1:rep) {
				tname = paste(n,i,sep='_')
				
				if ( is.null(x@usedObj[['rfExpressionSets']][[tname]]) ){
					## start the calculations!
					if ( dir.exists(opath)){
						if ( opath == '' ) {
							stop( "Are you mad? Not giving me an tmp path to delete?")
						}
						system( paste('rm -f ',opath,"/*",tname,'*', sep='') )
					}else {
						dir.create( opath )
					}
					total <- ncol(x@data)
					if ( total-subset <= 20  && rep > 1) {
						stop( paste( 'You have only', total, 'samples in this dataset and request to draw random',subset, "samples, which leaves less than 20 cells to draw on random!") )
					}
					else if ( total < subset ){
						stop ( paste("You can not ask for more than the max of",total, "samples in the test dataset!") )			
					}
					if ( is.null(x@usedObj[['rfExpressionSets']])){
						x@usedObj[['rfExpressionSets']] <- list()
						x@usedObj[['rfObj']][[ i ]] <- list()
					}
					
					if ( length( x@usedObj[['rfExpressionSets']] ) < i  ) {
						x@usedObj[['rfExpressionSets']][[ i ]] <- drop.samples( x, colnames(x@data)[sample(c(1:total),total-subset)], tname )
						if ( length(settings) > 0 ){
							 x@usedObj[['rfObj']][[ i ]] <- 
								 RFclust.SGE ( 
									dat=x@usedObj[['rfExpressionSets']][[ i ]]@data, 
									SGE=F, 
									slices=slice, 
									email=email, 
									tmp.path=opath, 
									name= tname,
								        slurm=T,
									settings=settings
									)
						}else {
							x@usedObj[['rfObj']][[ i ]] <- 
								RFclust.SGE ( 
									dat=x@usedObj[['rfExpressionSets']][[ i ]]@data, 
									SGE=SGE, 
									slices=slice, 
									email=email, 
									tmp.path=opath, 
									name= tname 
									)
						}

					}
					names(x@usedObj[['rfExpressionSets']]) [i] <- tname
					names(x@usedObj[['rfObj']]) [i] <- tname
					x@usedObj[['rfObj']][[ i ]] <- runRFclust ( x@usedObj[['rfObj']][[ i ]] , nforest=nforest, ntree=ntree, name=tname )
					if ( SGE){
						print ( "You should wait some time now to let the calculation finish! check: system('qstat -f') -> re-run the function")
					}
					else {
						print ( "You should wait some time now to let the calculation finish! -> re-run the function")
						print ( "check: system( 'ps -Af | grep Rcmd | grep -v grep')")
					}
				}
				else {
					
					## read in the results
					try ( x@usedObj[['rfObj']][[ i ]] <- runRFclust ( x@usedObj[['rfObj']][[ i]] , nforest=nforest, ntree=ntree, name=tname ) )
					if ( ! is.null(x@usedObj[['rfObj']][[ i ]]@RFfiles[[tname]]) ){
						stop( "please re-run this function later - the clustring process has not finished!")
					}
					for ( a in k ){
						x@usedObj[["rfExpressionSets"]][[i]]@samples <- 
								x@usedObj[["rfExpressionSets"]][[i]]@samples[ ,
										is.na(match ( colnames(x@usedObj[["rfExpressionSets"]][[i]]@samples), paste('group n=',a) ))==T 
								]
					}
					x <- createRFgrouping_col( x, RFname=tname,  k=10, single_res_col = paste( single_res_col, i) )
					
					print ( paste("Done with cluster",i))
					processed = TRUE
				}
			}
			x		
		}
)




#' @name createRFgrouping_col
#' @aliases createRFgrouping_col,StefansExpressionSet-method
#' @rdname createRFgrouping_col-methods
#' @docType methods
#' @description Create a sample grouping data from one RFclust.SGE object
#' @param x the StefansExpressionSet object
#' @param RFname the name of the RFclust.SGE object in the StefansExpressionSet object. This object has to be populized with data!
#' @param k the number of wanted groups ( default = 10)
#' @param single_res_col the new column in the samples table default= paste('RFgrouping', RFname)
#' @param colFunc a function giving the colours back for the grouping (gets the amount of groups) default = function(x){rainbow(x)}
#' @title description of function createRFgrouping_col
#' @export 
setGeneric('createRFgrouping_col', ## Name
		function ( x, RFname, k=10, single_res_col = paste('StefansExpressionSet',RFname), colFunc=NULL) { ## Argumente der generischen Funktion
			standardGeneric('createRFgrouping_col') ## der Aufruf von standardGeneric sorgt für das Dispatching
		}
)

setMethod('createRFgrouping_col', signature = c ('StefansExpressionSet'),
		definition = function ( x, RFname, k=10, single_res_col = paste('RFgrouping',RFname), colFunc=NULL) {
			if ( is.na( match( RFname, names(x@usedObj[['rfObj']])))){
				stop( paste("the RFname",RFname,"is not defined in this object; defined grouings are:",paste(names(x@usedObj[['rfObj']]), collapse=" ",sep=', ') ) )
			}
			browser()
			groups <- createGroups( x@usedObj[['rfObj']][[RFname]], k=k, name=RFname )
			x@usedObj[['rfExpressionSets']][[RFname]]@samples <- 
					cbind ( x@usedObj[['rfExpressionSets']][[RFname]]@samples, groups[,3:(2+length(k))] )
			le <- ncol(x@usedObj[['rfExpressionSets']][[RFname]]@samples)
			colnames(x@usedObj[['rfExpressionSets']][[RFname]]@samples)[(le-length(k)+1):le] <- 
					paste('group n=',k)
			m <- max(k)
			## create the predictive random forest object
			## this is not necessary if the whole object has been analyzed!
			if ( all.equal( colnames(x@usedObj[['rfObj']][[RFname]]@dat), colnames(x@data) ) ) {
				## use the column in grouping
				for ( id in 1:length(k) ){
					x@samples[, paste( single_res_col, ' n=', k[id], sep="") ] = factor(groups[,2+id], levels=c(1:k[id]))
					x <- colors_4( x, paste( single_res_col, ' n=', k[id], sep="")  )
				}
			}else {
				#predict based on the RFdata
				x@usedObj[['rfExpressionSets']][[RFname]] <- 
					bestGrouping( x@usedObj[['rfExpressionSets']][[RFname]], group=paste('group n=', m), bestColname = paste('OptimalGrouping',m ,RFname))
			
				x@samples[, paste( single_res_col) ] <-
					predict( x@usedObj[['rfExpressionSets']][[RFname]]@usedObj[[paste( 'predictive RFobj group n=',m) ]], t(as.matrix(x@data)) )
				x@samples[, paste( single_res_col) ] <- factor( x@samples[, paste( single_res_col) ], levels= 1:m )
				x <- colors_4( x, single_res_col )
			}
			x
		} 
)


