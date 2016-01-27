#library('ExpressionSet')

dat = data.frame( 'a' =c( 0:10), 'b' = c(10:20) ,ProbeSetID = c( paste ( 'probe_',c(1:11),sep='')), Crap=c(10:20) )
rownames(dat) <- dat$ProbeSetID
t <- ExpressionSet (
		dat = dat , 
		Samples= data.frame( 'GroupName' = paste('Group',c('a', 'b')), 'filename' = c('a','b')),
		namerow = 'ProbeSetID',
		outpath='./nowhere/',
		usecol = NULL
)

t
expect_equal(class(t)[1], 'ExpressionSet')

r <- melt(t, groupcol='GroupName', colCol='GroupName', probeNames="ProbeSetID")

