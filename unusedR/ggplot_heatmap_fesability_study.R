library('ggplot2')
library('reshape2')
library(grid)

colorpanel <- function(n, low, mid, high) 
{
    if (missing(mid) || missing(high)) {
            low <- col2rgb(low)
        if (missing(high)) 
            high <- col2rgb(mid)
        else high <- col2rgb(high)
        red <- seq(low[1, 1], high[1, 1], length = n)/255
        green <- seq(low[3, 1], high[3, 1], length = n)/255
        blue <- seq(low[2, 1], high[2, 1], length = n)/255
    }
    else {
            isodd <- n == round(n/2) *2
            if (isodd) {
                n <- n + 1
            }
            low <- col2rgb(low)
            mid <- col2rgb(mid)
            high <- col2rgb(high)
            lower <- floor(n/2)
            upper <- n - lower
            red <- c(seq(low[1, 1], mid[1, 1], length = lower), seq(mid[1, 
                1], high[1, 1], length = upper))/255
        green <- c(seq(low[3, 1], mid[3, 1], length = lower), 
               seq(mid[3, 1], high[3, 1], length = upper))/255
        blue <- c(seq(low[2, 1], mid[2, 1], length = lower), 
              seq(mid[2, 1], high[2, 1], length = upper))/255
        if (isodd) {
            red <- red[-(lower + 1)]
            green <- green[-(lower + 1)]
                blue <- blue[-(lower + 1)]
            }
    }
        rgb(red, blue, green)
}

bluered <- function(n){
colorpanel(n, "blue", "white", "red")
}
obj <- list (data = matrix ( sample( 0:12, 10000, replace=T), nrow=100 ) )
print ( dim(obj$data))
obj$samples <- data.frame( SampleName = paste( rep('Sample', 100), 1:100, sep='_') )
obj$samples$GroupName<- paste( 'Group',kmeans( t(obj$data), centers=10)$cluster)
obj$samples$GroupName<- factor( obj$samples$GroupName, levels= paste( 'Group',1:10))
colnames(obj$data) <- obj$samples$SampleName
rownames(obj$data) <-  paste('Gene', 1:100,sep='')
obj$annotation <- data.frame( ProbeName = rownames(obj$data), clusters = kmeans( obj$data, centers=10)$cluster )

ma  <- obj$data[,order(obj$samples[,'GroupName'] )]
melted <- melt( cbind(rownames(ma),ma) )
if ( length( which ( melted[,2] == '') ) > 0 ){
melted <- melted[ - which ( melted[,2] == ''),]
}
melted[,3] <- as.numeric(as.character(melted[,3]))
print(dim(melted))
grps<-NULL
for ( i in as.vector(obj$samples[,'GroupName']) ){
grps <- c( grps, rep( i, nrow(obj$data)))
}
cgrps <- NULL
for ( i in as.vector(obj$samples[,'GroupName']) ){
cgrps <- c( cgrps, rep( i, nrow(obj$data)))
}
colnames(melted) <- c('ProbeName', 'SampleName', 'Expression')

melted$Group <- factor(grps)
melted$ColorGroup <- factor(cgrps)

for ( n in c('ProbeName', 'SampleName' )){
melted[,n] <- factor(melted[,n], levels=unique(melted[,n]))
}


colrs <- rainbow( length(unique(melted$ColorGroup)) )
melted$colrss <- colrs[as.numeric(melted$ColorGroup)]
ss <-melted[which(melted$ProbeName==melted$ProbeName[1]),]
brks= c( -20.1, quantile(melted$Expression[which(melted$Expression != -20)],seq(0,1,by=0.1)) )
brks[length(brks)] = brks[length(brks)] + 0.1
melted$Expression <- cut( melted$Expression, breaks= brks)

## now I need to add the sample colors
melted_new <- NULL
melted_new <- as.matrix(melted[1:nrow(obj$data),])
x <- as.vector( t(melted[nrow(obj$data),]))
melted_new <- rbind( melted_new, matrix(c('SampleGroup', as.character(x[2]), as.character(x[5]), as.character(x[4]), as.character(x[5]), as.character(x[6]) ), nrow=1) )

for ( i in 2:ncol(obj$data)) {
melted_new <- rbind( melted_new, as.matrix(melted[(1+nrow(obj$data)*(i-1)):(nrow(obj$data)*i),]) )
x <- as.vector( t(melted[nrow(obj$data)*i,]))
melted_new <- rbind( melted_new, matrix(c('SampleGroup', as.character(x[2]), as.character(x[5]), as.character(x[4]), as.character(x[5]), as.character(x[6]) ), nrow=1) )
}

## and now add a gene level color code
add  <- as.matrix(melted[1:nrow(obj$data),])
add[,3] <- paste('Group', 1:10) [obj$annotation$clusters]
add[,4] <- add[,5] <-'GeneGroup'


rownames(add) <- (nrow(melted_new)+1):(nrow(melted_new)+nrow(add))
melted_new <- rbind( melted_new, add)

melted_df <- data.frame( melted_new)
melted_df$Expression <- factor( melted_df$Expression, levels = c( levels(melted$Expression), levels(obj$samples$GroupName)))
melted_df$clusters <- rep( c(obj$annotation$clusters, 'SG'), nrow(obj$data)+1 ) [1:nrow(melted_df)]

#melted_df$ProbeName <- factor( melted_df$ProbeName, levels= c( unique(melted[,n]), 'SampleGroup'))
melted_df$ProbeName <- factor( melted_df$ProbeName, levels= c('SampleGroup',  as.vector(obj$annotation$ProbeName)[order(obj$annotation$clusters)]))
melted_df$Group <- factor( melted_df$Group, levels= c(levels(obj$samples$GroupName), 'GeneGroup'))

for ( n in c('SampleName', 'ColorGroup' )){
	melted_df[,n] <- factor(melted_df[,n], levels=unique(melted[,n]))
}



p = ( ggplot(melted_df, aes(x=SampleName,y=ProbeName))
     	+ geom_tile(aes(fill=Expression)) 
+ scale_fill_manual( values = c( 'gray', bluered(10), rainbow(10) )) 
+ theme(
axis.text.x=element_blank(),
#axis.ticks.x=element_line(color=ss$colrss),
axis.ticks.length=unit(0.00,"cm")
)+ labs( y='') )
p <- p + facet_grid( clusters ~ Group,scales="free", space='free') 
print(p)

gg2tree <- function(gg) {
  la <- lapply(unlist(renquote(gg)), eval)
  # capture the output of the list values as one long character string
  vals <- lapply(la, function(x) paste(utils::capture.output(x), collapse = "<br>"))
    #names(vals) <- gsub("\\.", "_", names(vals))
    # preallocate matrix
    lvls <- strsplit(names(vals), "\\.")
    d <- max(sapply(lvls, length)) #maximum depth of the list
      m <- matrix(NA, nrow = length(lvls), ncol = d)
      for (i in seq_len(d)) m[,i]  <- sapply(lvls, function(x) x[i])
        m <- data.frame(m, value = as.character(vals))
        list(name = "ggplot", children = makeList(m))
}
#t <- gg2tree(p)


