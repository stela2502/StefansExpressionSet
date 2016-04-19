#' @name StefansExpressionSet
#' @title StefansExpressionSet
#' @description  An S4 class to visualize Expression data.
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot zscored genes are normalized?
#' @slot snorm samples normalized?
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot stats the stats list for all stats created in the object
#' @slot usedObj here a set of used and probably lateron important objects can be saved. Be very carful using any of them!
#' @exportClass StefansExpressionSet
setClass(
		Class='StefansExpressionSet', 
		representation = representation ( 
				data='data.frame',
				samples='data.frame',
				ranks='numeric',
				raw="data.frame",
				annotation='data.frame',
				snorm='logical',
				outpath='character',
				name='character',
				rownamescol='character',
				sampleNamesCol='character',
				stats = 'list',
				sig_genes = 'list',
				zscored = 'logical',
				simple = 'character',
				usedObj = 'list'
		),
		prototype(outpath ='', name = 'StefansExpressionSet',
				sampleNamesCol=NA_character_, 
				stats=list(),
				snorm=F,
				zscored=F,
				sig_genes=list(),
				usedObj= list(),
				simple= c( 'outpath', 'rownamescol', 'sampleNamesCol', 'simple', 'snorm', 'zscored') )
)


#' @name PMID25158935exp
#' @title Read counts for the expression data described in PMID25158935
#' @description The data was re-mapped against mouse mm10 using HISAT
#' @description and quantified using the R subreads package.
#' @docType data
#' @usage PMID25158935exp
#' @format data.frame
'PMID25158935exp'

#' @name PMID25158935samples
#' @title Read counts for the sample data for the expression information described in PMID25158935
#' @description The data was collected from the NCBI SRA archive
#' @docType data
#' @usage PMID25158935samples
#' @format data.frame
'PMID25158935samples'

#' @name red
#' @title reduced PMID25158935 dataset to a 100x15 StefansExpressionSet
#' @description Reduced StefansExpressionSet from the PMID25158935exp + PMID25158935samples dataset
#' @docType data
#' @usage red
#' @format StefansExpressionSet
'red'

#' @name NGSexpressionSet
#' @title NGSexpressionSet
#' @docType package
#' @description  An S4 class to visualize Expression data from a NGS experiment.
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot stats the stats list for all stats created in the object
#' @importClassesFrom StefansExpressionSet StefansExpressionSet
#' @exportClass NGSexpressionSet
setClass( 
		Class='NGSexpressionSet', 
		representation = representation ( 
		##	NGS = 'binary'
		),
		contains='StefansExpressionSet',
		prototype(outpath ='', name = 'NGSexpressionSet', 
				rownamescol=NA_character_, 
				sampleNamesCol=NA_character_, 
				stats=list() )
)

#' @name SingleCellsNGS
#' @title SingleCellsNGS
#' @docType package
#' @description  An S4 class to visualize Expression data from a single cell NGS experiment.
#' @slot data a data.frame containing the expression values for each gene x sample (gene = row)
#' @slot samples a data.frame describing the columnanmes in the data column
#' @slot annotation a data.frame describing the rownames in the data rows
#' @slot outpath the default outpath for the plots and tables from this package
#' @slot name the name for this package (all filesnames contain that)
#' @slot rownamescol the column name in the annotation table that represents the rownames of the data table
#' @slot sampleNamesCol the column name in the samples table that represents the colnames of the data table
#' @slot stats the stats list for all stats created in the object
setClass( 
		Class='SingleCellsNGS', 
		representation = representation ( 
		##	NGS = 'binary'
		),
		contains='NGSexpressionSet',
		prototype(outpath ='', name = 'SingleCellsNGS', 
				rownamescol=NA_character_, 
				sampleNamesCol=NA_character_, 
				stats=list() )
)


