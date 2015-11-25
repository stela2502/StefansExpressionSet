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
## v how many were mapped against this pathway?

bootstrap_phyper <- function(PW, l, n=2000, v, file.id, pathway.id  ) {
	## I need to get l random genes from the background (no lay back)
	
	f <- function ( bg, original, l, pw_genes ){
		#print ( paste( bg, l, pw_genes ) )
		lst <- sample( bg, l )
		pos <- sum(is.na( pw_genes )==F)
		OK <- sum( is.na( match(pw_genes, lst ) )==F)
		OK
	}
	if ( ! exists( 'boot', where=PW) ) {
		PW$boot <- list()
		for ( i in 1:length(names(PW$pw))){
			PW$boot[[i]] <- list()
			PW$boot[[i]][[length(names(PW$pw[[i]]))+1]] = -1
		}
	}
	if ( is.null(PW$boot[[file.id]][[pathway.id]]) ) {
		r <- list('R' = n, t <- vector('numeric', n) )
		range <- 1:PW$bgl
		pw_genes=match(PW$pw[[file.id]][[pathway.id]], PW$background)
		for ( i in 1:n ) { r$t[i] <- f( bg=range, l=l, pw_genes =  pw_genes ) }
		PW$boot[[file.id]][[pathway.id]] <- r
		#PW$boot[[file.id]][[pathway.id]] <-  boot(data=1:PW$bgl ,statistic=f, R=n, l=l, pw_genes=match(PW$pw[[file.id]][[pathway.id]], PW$background) )
		#PW$boot[[file.id]][[pathway.id]] <-  boot(data=PW$background,statistic=f, R=n, l=l, pw_genes=PW$pw[[file.id]][[pathway.id]] )
	}
	r <- length(which(PW$boot[[file.id]][[pathway.id]]$t >= v))
	if ( r == 0 ) {
		r =1
		#list( p= paste(">",( 1 /  PW$boot[[file.id]][[pathway.id]]$R ), sep='' ) ) 
	}
	#else {
		list( p= r /  PW$boot[[file.id]][[pathway.id]]$R )
	#}
}


testGeneList.pathwaySearch <- function (x, lst=c(),  bootstrap =F ){
	ret <- matrix(nrow=0, ncol=9 )
	colnames(ret) <- c('file', 'pathway', 'descr', 'pVal','adj.P BH', 'bootstrap', 'match', 'total','pw length' )
	for ( f in 1:length(names(x$pw)) ) {
		for (p in 1:length( names(x$pw[[f]]) ) ) {
			pos <- x$background[is.na(match( x$background, x$pw[[f]][[p]]))==F]
			if ( length(pos) > 5 ) {
				OK <- lst[is.na(fmatch( lst, x$pw[[f]][[p]]))==F]
				if ( length(OK) > 2 ) {
					#t <- bootstrap_phyper ( x, length(lst), v=length(OK), file.id=f, pathway.id=p )$p
					t<- (phyper( length(OK), length(pos), x$bgl - length(pos), length(lst),lower.tail=FALSE  ) * 20) ## emperical check using the bootsrap phyper is too positive
					bootst <- NA
					if ( t < 0.01 && bootstrap ) {
						bootst <- bootstrap_phyper ( x, length(lst), v=length(OK), file.id=f, pathway.id=p )$p
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




