# Read in and calculate Cnvrd2 on bam files



require(getopt)
require(CNVrd2)

bamHeader=read.table("bam_header_file.txt",header=T)
results=list()

read_matrices = list()
segment_scores = list()
polymorphic_regions = list()

bamHeader$Chromosome = as.character(bamHeader$Chromosome)

for(i in 1:length(bamHeader[,1])){
results[[i]] = new("CNVrd2", windows = 1000, st=0,en=as.numeric(bamHeader$Length[i]),chr = bamHeader$Chromosome[i], dirBamFile= paste('../by_chromosome/',bamHeader$Chromosome[i],sep=''),genes=c(0,1))
read_matrices[[i]] = countReadInWindow(Object = results[[i]])
segment_scores[[i]] = segmentSamples(Object=results[[i]],stdCntMatrix = read_matrices[[i]])
polymorphic_regions[[i]] = identifyPolymorphicRegion(Object = results[[i]],segmentObject = segment_scores[[i]])
}
for (i in 1:length(read_matrices)){
	png(paste("Lg",i+1,".polymorphic.png",sep=""))
	plotPolymorphicRegion(Object = results[[i]], polymorphicRegionObject = polymorphic_regions[[i]],drawThresholds=TRUE)
	dev.off()
}

for (i in 1:length(read_matrices)){
	polymorphic_regions[[i]] = identifyPolymorphicRegion(Object = results[[i]],segmentObject = segment_scores[[i]])
}
lg2=polymorphic_regions[[1]]

segment_scores[[1]] = segmentSamples(Object=results[[1]],stdCntMatrix = read_matrices[[1]])
polymorphic_regions[[1]] = identifyPolymorphicRegion(Object = results[[1]],segmentObject = segment_scores[[1]])


# filter the CNVrd2 results

pdf("Lg1polymorphicHist.pdf")
filtered_5000length = polymorphic_regions[[1]]$SSofPolymorphicRegions[polymorphic_regions[[1]]$SSofPolymorphicRegions[,2] - polymorphic_regions[[1]]$SSofPolymorphicRegions[,1] > 5000,]
apply(filtered_5000length,1,function(x) hist(x[-c(1,2) ],breaks=100))
dev.off()

require(rtracklayer)
malus_gff = import.gff("Malus_x_domestica.v1.0.consensus.gff")
# Look for genes that fall in the regions
total_CDS_count=0
total_MRNA_count=0

