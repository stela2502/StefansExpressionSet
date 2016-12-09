#' @name read.bams
#' @aliases read.bams,NGSexpressionSet-method
#' @rdname read.bams-methods
#' @docType methods
#' @description  Use most of'StefansExpressionSet'functionallity with minor changes to NGS data (like normalizing)
#' @description  This package is mainly meant to plot the data - all important quantification or gene
#' @description  annotation is performed using DEseq functionallity. Read a set of bam files and perform
#' @description  the quantification (better do that without using this function!)
#' @param bamFiles a file containing one file name per line
#' @param annotation the gene level annotation to use
#' @param as.obj if true return the couts object instead of a matrix (default =F)
#' @title description of function read.bams
#' @export 
setGeneric('read.bams', ## Name
	function ( bamFiles, annotation, GTF.featureType='exon', GTF.attrType = "gene_id", isPairedEnd = FALSE, nthreads = 2, as.obj=F) { ## Argumente der generischen Funktion
		standardGeneric('read.bams') ## der Aufruf von standardGeneric sorgt für das Dispatching
	}
)

setMethod('read.bams', signature = c ('character') ,
	definition = function ( bamFiles, annotation, GTF.featureType='exon', GTF.attrType = "gene_id", isPairedEnd = FALSE, nthreads = 2, as.obj=F) {
	if (file.exists(bamFiles)){
		bamFiles <- readLines(bamFiles)
	}
	# break that into 500 samples chunks

	counts <- featureCounts(files =bamFiles,annot.ext = annotation ,isGTFAnnotationFile = TRUE,GTF.featureType = GTF.featureType,
		GTF.attrType = GTF.attrType,allowMultiOverlap=T, isPairedEnd =isPairedEnd , nthreads = nthreads)
	save(counts, file="count_object.RData" )
	if ( ! as.obj ) {
		counts
	}else {
		counts.tab <- cbind(counts$annotation,counts$counts)  # combine the annotation and counts
		counts.tab
	}
})


