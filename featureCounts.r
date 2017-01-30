"""Christa Caggiano 
Jan 2017 
Feature Counts command parameters for read counts on whole brain 
atac data
"""


library("Rsubread") 

read_files <- c("FE165_nodup.sort.bam", "FE167_nodup.sort.bam", "FE177_nodup.sort.bam", "FE178_nodup.sort.bam")

fc = featureCounts(read_files,

	# annotation
	annot.ext="intact_peaks.saf",
	isGTFAnnotationFile=FALSE,
	GTF.featureType="exon",
	GTF.attrType="gene_id",
	chrAliases=NULL,
	
	# level of summarization
	useMetaFeatures=TRUE,
	
	# overlap between reads and features
	allowMultiOverlap=FALSE,

	# multi-mapping reads
	countMultiMappingReads=FALSE,

	# read filtering
	minMQS=0,
	
	# strandness
	strandSpecific=0,
	
	# parameters specific to paired end reads
	isPairedEnd=FALSE,
	requireBothEndsMapped=FALSE,
	checkFragLength=FALSE,
	minFragLength=50,
	maxFragLength=600,
	countChimericFragments=TRUE,	
	
	# miscellaneous
	nthreads=1,
	# maxMOp=10,
	reportReads=FALSE)

write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],fc$counts,stringsAsFactors=FALSE),file="counts.txt",quote=FALSE,sep="\t",row.names=FALSE)

