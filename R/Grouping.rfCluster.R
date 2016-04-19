#' @name rfCluster
#' @aliases 'rfCluster,SingleCellsNGS-method
#' @title rfCluster
#' @name rfCluster-methods
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
#' @param pics create a heatmap for each grouping that has been accessed (in the outpath folder; default = FALSE)
#' @param nforest the numer of forests to grow for each rep (defualt = 500)
#' @param ntree the numer of trees per forest (default = 500)
#' @return a SingleCellsNGS object including the results and storing the RF object in the usedObj list (bestColname)
#' @export 
setGeneric('rfCluster',
		function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, pics=F ,nforest=500, ntree=500, name='RFclust'){
			standardGeneric('rfCluster')
		}
)
setMethod('rfCluster', signature = c ('SingleCellsNGS'),
		definition = function ( x, rep=5, SGE=F, email, k=16, slice=30, subset=200, pics=F ,nforest=500, ntree=1000, name='RFclust') {
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
					if ( total-subset <= 20 ) {
						stop( paste( 'You have only', total, 'samples in this dataset and request to draw random',subset, "samples, which leaves less than 20 cells to draw on random!") )
					}
					if ( is.null(x@usedObj[['rfExpressionSets']])){
						x@usedObj[['rfExpressionSets']] <- list()
						x@usedObj[['rfObj']][[ i ]] <- list()
					}
					
					if ( length( x@usedObj[['rfExpressionSets']] ) < i  ) {
						x@usedObj[['rfExpressionSets']][[ i ]] <- drop.samples( x, colnames(x@data)[sample(c(1:total),total-subset)], tname )
						x@usedObj[['rfObj']][[ i ]] <- RFclust.SGE ( dat=x@usedObj[['rfExpressionSets']][[ i ]]@data, SGE=SGE, slice=slice, email=email, tmp.path=opath, name= tname )
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
					groups <- createGroups( x@usedObj[['rfObj']][[i]], k=k, name=tname )
					x@usedObj[['rfExpressionSets']][[i]]@samples <- cbind ( x@usedObj[['rfExpressionSets']][[i]]@samples, groups[,3:(2+length(k))] )
					
					le <- ncol(x@usedObj[['rfExpressionSets']][[i]]@samples)
					colnames(x@usedObj[['rfExpressionSets']][[i]]@samples)[(le-length(k)+1):le] <- paste('group n=',k)
					
					## create the required RF object
					m <- max(k)
					x@usedObj[['rfExpressionSets']][[i]] <- bestGrouping( x@usedObj[['rfExpressionSets']][[i]], group=paste('group n=', m), bestColname = paste('OptimalGrouping',m ,name) )
					## the 'predictive RFobj group n=' object is created by the bestGrouping call
					x@samples[, paste( single_res_col, i) ] <-
							predict( x@usedObj[['rfExpressionSets']][[i]]@usedObj[[paste( 'predictive RFobj group n=',m) ]], t(as.matrix(x@data)) )
					if ( pics ){
						fn <- paste(OPATH,'/heatmap_rfExpressionSets_',i,'.png', sep='')
						png ( file=fn, width=800, height=1600 )
						gg.heatmap.list( x, groupCol=paste( single_res_col , i) )
						dev.off()
						print ( paste('heatmap stored in', fn) )
					}
					print ( paste("Done with cluster",i))
					processed = TRUE
				}
			}
			if ( processed ) {
				try ( {combine <- identifyBestGrouping( x, c( paste(single_res_col, 1:rep)) )} , silent =T)
				if ( all.equal(as.vector(combine$res), rep('', rep)) ) {
					print( 'No really usful grouping of the data obtained - I recommend re-run with more trees/forests and a new name')
					x <- combine$x
				}
				else {
					x <- combine$x
					colnames(x@samples)[which( colnames(x@samples) == names(combine$res)[1] )] <- usefulCol
					if ( pics ){
						fn <- paste(OPATH,'/heatmap_',str_replace( usefulCol, '\\s', '_'),'.png', sep='')
						png ( file=fn, width=800, height=1600 )
						gg.heatmap.list( x, groupCol= usefulCol )
						dev.off()
						print ( paste('heatmap stored in', fn ))
					}
				}
				x@usedObj$combinationAnalysis <- list ( 'initial_significants' = combine$names, 'merged_significants' = combine$res )				
			}
			x		
		}
)

