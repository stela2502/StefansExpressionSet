
PMID25158935 <- ExpressionSet( PMID25158935exp, PMID25158935samples,  Analysis = NULL, name='PMID25158935', namecol='Sample', namerow= 'GeneID', usecol=NULL , outpath = '../../tmp/')

expect_equal(class(PMID25158935)[1], 'ExpressionSet' )
expect_equal( dim(PMID25158935@data), c(24062,15) )

expect_equal( dim(PMID25158935@samples), c( 15,21) )
expect_equal( PMID25158935@rownamescol, "GeneID" )
expect_equal( PMID25158935@sampleNamesCol , "Sample" )
expect_equal( PMID25158935@name, 'PMID25158935' )

red <- reduce.Obj ( PMID25158935, rownames(PMID25158935@data)[1:100], name='minimal' )

expect_equal( dim(red@data), c(100,15) )
expect_equal( red@name, 'minimal' )

expect_equal( as.matrix(melt(red, probeNames='GeneID')[1,]) , 
		as.matrix(data.frame( ProbeName = 23343, SampleName= 'ERR420375', Expression = 0, Group = 'HSC', ColorGroup='HSC')[1,]) )

if ( file.exists( red@outpath, 'reducedSet_points_Mybl1_expression.svg', sep='/' ) ){
 system ( paste ('rm ', red@outpath, 'reducedSet_points_Mybl1_expression.svg',sep='' ) )
}


## test the first plot.probeset function

p <- plot.probeset( red, rownames(red@data)[20], geneNameCol= 'GeneID')

expect_true( file.exists( paste(red@outpath, 'reducedSet_points_Mybl1_expression.svg', sep='/' ) ))

if ( file.exists( red@outpath, 'reducedSet_points_Mybl1_expression.svg', sep='/' ) ){
	system ( paste ('rm ', red@outpath, 'reducedSet_points_Mybl1_expression.svg',sep='' ) )
}

if ( file.exists( red@outpath, 'reducedSet_boxplot_Mybl1_expression.svg', sep='/' ) ){
	system ( paste ('rm ', red@outpath, 'reducedSet_boxplot_Mybl1_expression.svg',sep='' ) )
}

p <- plot.probeset( red, rownames(red@data)[20], geneNameCol= 'GeneID', boxplot=T)

expect_true( file.exists( paste(red@outpath, 'reducedSet_boxplot_Mybl1_expression.svg', sep='/' ) ))

if ( file.exists( red@outpath, 'reducedSet_boxplot_Mybl1_expression.svg', sep='/' ) ){
	system ( paste ('rm ', red@outpath, 'reducedSet_boxplot_Mybl1_expression.svg',sep='' ) )
}

dropS <- drop.samples( red, samplenames=red@samples[1:5, red@sampleNamesCol])
expect_equal( dim(dropS@data), c(100,10) )
expect_equal( dim(dropS@samples), c( 10,21) )
expect_equal( dim(dropS@annotation), c( 100,2) )

p <- gg.heatmap.list(dropS, groupCol='GroupName')





