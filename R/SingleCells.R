


z.score.SingleCells <- function ( x ) {
	ma  <- as.matrix(x$data)
	i = 0
	ret <- t(
			apply(ma,1, function (x) {
						i = i+1
						n <- which(x==0)
						if ( length(x) - length(n) > 1 ){
							x[-n] <- scale(as.vector(t(x[-n])))
						}
						else {
							x[] = -20
						}
						x[n] <- -20
						x}
			)
	)
	x$data <- ret
	x
}

gg.heatmap <- function (dat, colrs, groupCol='GroupName', colCol='GroupName' ) {
	UseMethod('gg.heatmap',dat )
}


gg.heatmap.SingleCells <- function(dat, colrs, groupCol='GroupName', colCol='GroupName' ){
	
	dat <- z.score ( dat )
	melted <- melt(dat)
browser()
	samp.cast <- dcast(melted,ProbeName~SampleName,mean,value.var="Expression")
	samp.mat <- as.matrix(samp.cast[,2:ncol(samp.cast)])
	ord.genes <-
			as.vector(samp.cast[hclust(dist(samp.mat),method="ward.D")$order,1])
	
	for ( n in c('ProbeName', 'SampleName', 'Group', 'ColorGroup' )){
		melted[,n] <- factor(n, levels=unique(melted[,n]))
	}
	
	melted$colrss <- colrs[melted$ColorGroups]
	
	ss <-dat.ss[which(dat.ss$ProbeName==dat.ss$ProbeName[1]),]
	
	brks= c( -20.1, quantile(melted$Expression[which(melted$Expression != -20)],seq(0,1,by=0.1)) )
	brks[length(brks)] = brks[length(brks)] + 0.1
	melted$Expression <- cut( melted$Expression, breaks= brks)
	plot = ggplot(dat.ss, aes(x=SampleName,y=ProbeName))
	+ geom_tile(aes(fill=Expression))
	+ scale_fill_manual( values = c( 'gray', bluered(10)) ) 
	+ theme(
			legend.position= 'bottom',
			axis.text.x=element_blank(),
			axis.ticks.x=element_line(color=ss$colrss),
			axis.ticks.length=unit(0.6,"cm")
	)
	+ labs( y='')
	print ( plot )
	invisible (plot)
}