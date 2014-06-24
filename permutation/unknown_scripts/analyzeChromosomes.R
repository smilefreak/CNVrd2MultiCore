# Read in and calculate Cnvrd2 on bam files


require(getopt)
require(CNVrd2)

bamHeader=read.table("bam_header_file.txt",header=T)
results=list()

read_matrices = list()
segment_scores = list()
polymorphic_regions = list()

bamHeader$Chromosome = as.character(bamHeader$Chromosome)

segmentSamplesWrapper = function(x,y){
	return(segmentSamples(Object=x,stdCntMatrix=y))
}
identifyPolymorphicRegionWrapper = function (x,y){
	return(identifyPolymorphicRegion(Object = x,segmentObject = y))
}


for(i in 1:length(bamHeader[,1])){
results[[i]] = new("CNVrd2", windows = 1000, st=0,en=as.numeric(bamHeader$Length[i]),chr = bamHeader$Chromosome[i], dirBamFile= paste('../by_chromosome/',bamHeader$Chromosome[i],sep=''),genes=c(0,1))
}
read_matrices = lapply(results,FUN=function(x){return(countReadInWindow(Object = x))}) 
segment_scores= mapply(segmentSamplesWrapper,results,read_matrices)
polymorphic_regions = mapply(identifyPolymorphicRegionWrapper,results,segment_scores)
for (i in 1:length(read_matrices)){
	png(paste("Lg",i+1,".polymorphic.png",sep=""))
	plotPolymorphicRegion(Object = results[[i]], polymorphicRegionObject = polymorphic_regions[[i]],drawThresholds=TRUE)
	dev.off()
}

for (i in 1:length(read_matrices)){
	segment_scores[[i]] = segmentSamples(Object=results[[i]],stdCntMatrix = read_matrices[[i]])
	polymorphic_regions[[i]] = identifyPolymorphicRegion(Object = results[[i]],segmentObject = segment_scores[[i]],plotLegend =F)
}
for (i in 1:length(read_matrices)){
	polymorphic_regions[[i]] = identifyPolymorphicRegion(Object = results[[i]],segmentObject = segment_scores[[i]],plotLegend =F)
}
lg2=polymorphic_regions[[1]]

segment_scores[[1]] = segmentSamples(Object=results[[1]],stdCntMatrix = read_matrices[[1]])
polymorphic_regions[[1]] = identifyPolymorphicRegion(Object = results[[1]],segmentObject = segment_scores[[1]])


# filter the CNVrd2 results

pdf("Lg1polymorphicHist.pdf")
filtered_5000length = polymorphic_regions[[1]]$SSofPolymorphicRegions[polymorphic_regions[[1]]$SSofPolymorphicRegions[,2] - polymorphic_regions[[1]]$SSofPolymorphicRegions[,1] > 5000,]
apply(filtered_5000length,1,function(x) hist(x[-c(1,2) ],breaks=100))
dev.off()


# gets all the genes that toally fit in the predicted CNV region
get_genes_in_regions = function(x,gff){
	gff_ranges =as.data.frame(ranges(gff))
	CDS_in_region=list()
	no_items=1
	for ( i in 1:nrow(gff_ranges)){
		if(x[1] <= gff_ranges$start[i] && gff_ranges$end[i] <= x[2]){
						CDS_in_region[[no_items]]=gff[i]
						no_items = no_items + 1
		}
	}
	return(CDS_in_region)
}


require(rtracklayer)
malus_gff =import.gff3("Malus_x_domestica.v1.0.consensus.fixed.gff",asRangedData=F)
# Look for genes that fall in the regions
total_CDS_count=0
#Start with CDS.
malus_gff_cds=malus_gff[which(malus_gff$type=="CDS"),]
one_region=(polymorphic_regions[[1]]$putativeBoundary)
ranged_regions=(GRanges("chr2",one_region),strand="+")
intersect(ranged_regions,chr2_gff_cds)
regions_data_frame=as.data.frame(polymorphic_regions[[1]]$putativeBoundary)
matrix_putative_boundaries=apply(as.data.frame(polymorphic_regions[[1]]$putativeBounday),1,get_genes_in_regions,gff=chr2_gff_cds) 
list_putative_boundaries=list()

#create chromosome by chromosome #

list_CDS=list()
for (i in 2:17){
	list_CDS[[i-1]]=malus_gff_cds[which(as.vector(seqnames(malus_gff_cds)) == paste('chr',i,sep='')),]
}
list_putative_boundaries_CDS=list()
for ( i in 1:length(polymorphic_regions)){
	chr_data_frame=as.data.frame(polymorphic_regions[[i]]$putativeBoundary)
	list_putative_boundaries_CDS[[i]]=apply(chr_data_frame,1,get_genes_in_regions,gff=list_CDS[[i]]) 
}
total_MRNA_count=0
cds = get_genes_in_regions(one_region,chr2_gff_cds)

#Do MRNA
malus_gff_mrna=malus_gff[which(malus_gff$type=="mRNA"),]
list_mRNA=list()
for (i in 2:17){
	list_mRNA[[i-1]]=malus_gff_mrna[which(as.vector(seqnames(malus_gff_mrna)) == paste('chr',i,sep='')),]
}

list_putative_boundaries_mrna=list()
for ( i in 1:length(polymorphic_regions)){
	chr_data_frame=as.data.frame(polymorphic_regions[[i]]$putativeBoundary)
	list_putative_boundaries_mrna[[i]]=apply(chr_data_frame,1,get_genes_in_regions,gff=list_mRNA[[i]]) 
}

#remove nothing matches
#off set by one chr2 is position 1 in all the lists
library(GenomicRanges)
mrna_hits=list()
do_plots  <- function(list_putative_boundaries,polymorphic_regions,postfix){
for(i in 1){
pdf(paste("Lg",i+1,postfix,'.pdf',sep=''))
for ( j in 1:length(list_putative_boundaries[[i]])){
			if(length(list_putative_boundaries[[i]][[j]]) !=0){
				limits=c(start(polymorphic_regions[[i]]$putativeBoundary[j])-20000,end(polymorphic_regions[[i]]$putativeBoundary[j])+20000)
				hi = ggplot(data.frame(x=as.numeric(polymorphic_regions[[i]]$SSofPolymorphicRegions[j,-c(1,2) ])),aes(x=x)) +geom_histogram(bin=0.05)
				plot(hi)
				print(j)
				temp_list = clone(list_putative_boundaries_mrna[[i]][[j]][[1]])
				if(length(list_putative_boundaries[[i]][[j]]) > 1){
			
				for( z in 2:length(list_putative_boundaries[[i]][[j]])){
						temp_list = append(temp_list,list_putative_boundaries[[i]][[j]][[z]])
				}	
				}
				granges_putative=GRanges(seqnames=(paste("chr",i+1,sep='')),strand=c("*"),polymorphic_regions[[i]]$putativeBoundary[j])
				rec = ggplot(list_putative_boundaries[[i]][[j]][[1]],aes(color=factor(list_putative_boundaries[[i]][[j]][[1]]$ID)))
				rec = rec + geom_rect(data=polymorphic_regions[[i]]$putativeBoundary[j])
				rec = rec + geom_rect(data=list_putative_boundaries[[i]][[j]][[1]],aes(color=factor(list_putative_boundaries[[i]][[j]][[1]]$ID),alpha=0.5,fill=factor(list_putative_boundaries[[i]][[j]][[1]]$ID)))
				rec = rec + xlim(c(start(polymorphic_regions[[i]]$putativeBoundary[j])-20000,end(polymorphic_regions[[i]]$putativeBoundary[j])+20000))
				print(i)
				plot(rec)
				p = plotPolymorphicRegion(Object = results[[i]], polymorphicRegionObject = polymorphic_regions[[i]],drawThresholds=TRUE,xlim=limits)
		}
	}
dev.off()
}
}
library(ggplot)


print(i)
plot_all_chromosomes = function(polymorphic_regions,results){
pdf("all.chromosomes.polymorphicplot.pdf")
for(i in 1:length(polymorphic_regions)){
plotPolymorphicRegion(Object = results[[i]], polymorphicRegionObject = polymorphic_regions[[i]],drawThresholds=TRUE,thresholdForPolymorphicRegions=c(0.75,0.25))
}
dev.off()
}
do_plots(list_putative_boundaries_mrna,polymorphic_regions,"mRNAPlots")
do_plots(list_putative_boundaries_CDS,polymorphic_regions,"CDSPlots")
plot_all_chromosomes(polymorphic_regions,results)
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#get mrNA plots
for(i in 1:length(polymorphic_regions)){
pdf(paste("Lg",i+1,"mrna",".pdf",sep=""))
for ( j in 1:length(list_putative_boundaries_mrna[[i]])){
			if(length(list_putative_boundaries_mrna[[i]][[j]]) !=0){
				limits=c(start(polymorphic_regions[[i]]$putativeBoundary[j])-20000,end(polymorphic_regions[[i]]$putativeBoundary[j])+20000)
				hi = ggplot(data.frame(x=as.numeric(polymorphic_regions[[i]]$SSofPolymorphicRegions[j,-c(1,2) ])),aes(x=x)) +geom_histogram(bin=0.05)
				temp_list = list_putative_boundaries_mrna[[i]][[j]][[1]]
				if(length(list_putative_boundaries_mrna[[i]][[j]]) > 1){
				
				for( z in 2:length(list_putative_boundaries_mrna[[i]][[j]])){
						temp_list = append(temp_list,list_putative_boundaries_mrna[[i]][[j]][[z]])
				}	
				}
				granges_putative=GRanges(seqnames=(paste("chr",i+1,sep='')),strand=c("*"),polymorphic_regions[[i]]$putativeBoundary[j])
				rec = ggplot(temp_list,aes(color=factor(temp_list$ID)))
				rec = rec + geom_rect(data=polymorphic_regions[[i]]$putativeBoundary[j])
				rec = rec + geom_rect(data=temp_list,aes(color=factor(temp_list$ID),alpha=0.5,fill=factor(temp_list$ID)))
				rec = rec + xlim(c(start(polymorphic_regions[[i]]$putativeBoundary[j])-20000,end(polymorphic_regions[[i]]$putativeBoundary[j])+20000))
				multiplot(rec,hi)
				p = plotPolymorphicRegion(Object = results[[i]], polymorphicRegionObject = polymorphic_regions[[i]],drawThresholds=TRUE,xlim=limits)
		}
	}
dev.off()
}

# get CDS plots
for(i in 1:length(polymorphic_regions)){
pdf(paste("Lg",i+1,"cds",".pdf",sep=""))
for ( j in 1:length(list_putative_boundaries_CDS[[i]])){
			if(length(list_putative_boundaries_CDS[[i]][[j]]) !=0){
				limits=c(start(polymorphic_regions[[i]]$putativeBoundary[j])-20000,end(polymorphic_regions[[i]]$putativeBoundary[j])+20000)
				hi = ggplot(data.frame(x=as.numeric(polymorphic_regions[[i]]$SSofPolymorphicRegions[j,-c(1,2) ])),aes(x=x)) +geom_histogram(bin=0.05)
				temp_list = list_putative_boundaries_CDS[[i]][[j]][[1]]
				if(length(list_putative_boundaries_CDS[[i]][[j]]) > 1){
				
				for( z in 2:length(list_putative_boundaries_CDS[[i]][[j]])){
						temp_list = append(temp_list,list_putative_boundaries_CDS[[i]][[j]][[z]])
				}	
				}
				granges_putative=GRanges(seqnames=(paste("chr",i+1,sep='')),strand=c("*"),polymorphic_regions[[i]]$putativeBoundary[j])
				rec = ggplot()
				rec = rec + geom_rect(data=polymorphic_regions[[i]]$putativeBoundary[j])
				rec = rec + geom_rect(data=temp_list,aes(color=factor(temp_list$ID),alpha=0.5,fill=factor(temp_list$ID)))
				rec = rec + xlim(c(start(polymorphic_regions[[i]]$putativeBoundary[j])-20000,end(polymorphic_regions[[i]]$putativeBoundary[j])+20000))
				multiplot(rec,hi)
				p = plotPolymorphicRegion(Object = results[[i]], polymorphicRegionObject = polymorphic_regions[[i]],drawThresholds=TRUE,xlim=limits)
		}
	}
dev.off()
}


# Script to cluster the results from the CNVrd2 manual

objectCluster <- new("clusteringCNVs", x = segment_scores[[1]]$segmentationScores[,1], k = 3, EV =TRUE)
copyNumberGroups  <- groupCNVs(Object = objectCluster, rightLimit = 1.5)
allGroups  <-  copyNumberGroups$allGroups
duplicatedSamples  <- rownames(allGroups[allGroups[,2] > 2,])

dir.create("readCounts")

for (i in 1:length(results)){
for (ii in rownames(allGroups)){
	dir.create(paste("readCounts/Lg",i+1,sep=''))
	png(paste("readCounts","/","Lg",i+1,'/',ii,".png",sep=""))
	plotCNVrd2(Object = results[[1]],segmentObject = segment_scores[[1]], stdCntMatrix = read_matrices[[1]], sampleName = ii,data1000Genomes=F)
	dev.off()
}
command_image_magick=paste("composite readCounts/Lg",i + 1,'/',"*.png",sep="") 
system(command_image_magick)
}


#!?!?!?!?? PERMUTE THE SAMPLES
