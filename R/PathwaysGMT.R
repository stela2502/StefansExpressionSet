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

bootstrap_phyper <- function(PW, l, n=3000, v, f, p  ) {
	## I need to get l random genes from the background (no lay back)
	f <- function ( bg, original, l, pw_genes ){
		#print ( paste( bg, l, pw_genes ) )
		lst <- sample( bg, l )
		pos <- bg[is.na(match( bg, pw_genes))==F]
		OK <- lst[is.na(match( lst, pw_genes ))==F]
		
		if ( length(OK) < 2 ) {
			c(1)
		}
		else {
			c(phyper( length(OK), length(pos), length(bg) - length(pos), l,lower.tail=FALSE  ))
		}
	}
	if ( ! exists( 'boot', where=PW) ) {
		PW$boot <- list()
		for ( i in 1:length(names(PW$pw))){
			PW$boot[[i]] <- list()
			PW$boot[[i]][[length(names(PW$pw[[i]]))+1]] = -1
		}
	}
	if ( is.null(PW$boot[[f]][[p]]) ) {
		PW$boot[[f]][[p]] <-  boot(data=PW$background,statistic=f, R=n, l=l, pw_genes=PW$pw[[f]][[p]] )
	}
	list( p= length(which(PW$boot[[f]][[p]]$t < v)) )
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
						bootst <- bootstrap_phyper ( x, length(lst), v=t, f=f, p=p )
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




