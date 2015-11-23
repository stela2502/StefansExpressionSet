# TODO: Add comment
# 
# Author: stefan
###############################################################################



readGMT <- function (inputFile) 
{
	f <- readLines(inputFile)
	lst <- sapply(f, function(x) unlist(strsplit(x, "\t", fixed = TRUE)))
	names(lst) <- sapply(lst, function(x) x[1])
	link <- unlist(t((lapply( lst, function(x) x[2] ))))
	lst <- lapply(lst, function(x) x[-(1:2)])
	attributes(lst)$link <- link
	return(lst)
}

## read all gmt files in a path into the object and return the pathwaySearch object

start <- function ( path=pwd(), background=c() ) {
	ret <- list()
	if ( length(background) > 0 ){
		ret$background <- background
		ret$bgl <- length(ret$background)
	}else {
		stop( "I need a list of of all analyzed genes in the dataset" )
	}
	ret$fn <- list.files( path=path, pattern= ".gmt$",  all.files = T, full.names = T)
	ret$pw <- list()
	i =1
	for (f in  ret$fn ){
		ret$pw[[i]] <- readGMT ( f )
		i = i+1
	}
	names( ret$pw ) <- ret$fn
	
	class(ret) <- 'pathwaySearch'
	ret
}

testGeneList <- function (x, lst=c() ){
	UseMethod('testGeneList', x)
}

testGeneList.pathwaySearch <- function (x, lst=c() ){
	ret <- matrix(nrow=0, ncol=8 )
	colnames(ret) <- c('file', 'pathway', 'descr', 'pVal','adj.P BH', 'match', 'total','pw length' )
	for ( f in 1:length(names(x$pw)) ) {
		for (p in 1:length( names(x$pw[[f]]) ) ) {
			pos <- x$background[is.na(match( x$background, x$pw[[f]][[p]]))==F]
			if ( length(pos) > 5 ) {
				OK <- lst[is.na(match( lst, x$pw[[f]][[p]]))==F]
				if ( length(OK) > 2 ) {
					ret <- rbind(ret, c(
									x$fn[[f]], 
									names(x$pw[[f]])[p], 
									attributes(x$pw[[f]])$link[p],  
									phyper( length(OK), length(pos), x$bgl - length(pos), length(lst),lower.tail=FALSE  ), 
									0,
									length(OK), 
									length(lst), 
									length(pos), 
									paste( OK, collapse=";") )
					) 		
				}
			}
		}
	}
	ret[,5] <- p.adjust( ret[,4], method ='BH' )
	ret <- ret[order(ret[,5]),]
	ret
}	
t <- testGeneList( PW, toupper(geneNames_4_Group ( HSC_single_cells, userGroups,  gId=10, nameCol='Gene.Symbol' )))




