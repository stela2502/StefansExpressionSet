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

## PW = PathwaysGMT obj
## l  = amount of genes from the original list that match to the background

bootstrap_phyper <- function(PW, l, n=1000, v, pw_genes ) {
	## I need to get l random genes from the background (no lay back)
	f <- function ( bg, l, pw_genes ){
		lst <- sample( bg, l )
		pos <- bg[is.na(match( bg, pw_genes))==F]
		OK <- lst[is.na(match( lst, pw_genes ))==F]
		phyper( length(OK), length(pos), x$bgl - length(pos), l,lower.tail=FALSE  )
	}
	
	tmp <- boot(PW$background,f, R=500, l=l,  pw_genes=pw_genes )
	tmp
}


testGeneList.pathwaySearch <- function (x, lst=c() ){
	ret <- matrix(nrow=0, ncol=9 )
	colnames(ret) <- c('file', 'pathway', 'descr', 'pVal','adj.P BH', 'bootstrap', 'match', 'total','pw length' )
	for ( f in 1:length(names(x$pw)) ) {
		for (p in 1:length( names(x$pw[[f]]) ) ) {
			pos <- x$background[is.na(match( x$background, x$pw[[f]][[p]]))==F]
			if ( length(pos) > 5 ) {
				OK <- lst[is.na(match( lst, x$pw[[f]][[p]]))==F]
				if ( length(OK) > 2 ) {
					t<- phyper( length(OK), length(pos), x$bgl - length(pos), length(lst),lower.tail=FALSE  )
					bootst <- NA
					if ( t < 0.05 ) {
						
					}
					ret <- rbind(ret, c(
									x$fn[[f]], 
									names(x$pw[[f]])[p], 
									attributes(x$pw[[f]])$link[p],  
									t, 
									0,
									bootst,
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




