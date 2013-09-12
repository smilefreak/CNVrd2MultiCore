require(CNVrd2)
require(parallel)
cores=1

permuteSample <- function(segmentationResultsForSample,st,windows){
		ranks=seq(1,nrow(segmentationResultsForSample))	
		reorder_ranks=sample(ranks,replace=F)
		segmentationResultsForSample = segmentationResultsForSample[reorder_ranks,]	
		# Improve speed of this horrible for loop lol
		for( i in 1:nrow(segmentationResultsForSample)){
			dif = segmentationResultsForSample[i,4] - segmentationResultsForSample[i,3] + windows
			segmentationResultsForSample[i,3] = st
			segmentationResultsForSample[i,4] = st  + dif - windows
			st = st + dif
		}	
		return(segmentationResultsForSample)
}

kstestwrapper = function(x,y){
		polyMorphicResampling[[2]]$subRegionMatrix[1,]
		mapply(ks.test,split(x$subRegionMatrix,col(x$subRegionMatrix)),split(y$subRegionMatrix,col(y$subRegionMatrix)))
}

# Permutes the regions that have been found by Binary segmentation.
permuteSegmentationScores  <-  function(segment_scores,cnvrd2Object){
	ranks = seq(1,length(segment_scores[[1]]))
	samples = rownames(segment_scores$segmentationScores)
	#pull out samples
	sample_list=lapply(samples, function(x) segment_scores[[1]]$segmentResults[[which(rownames(segment_scores[[1]]$segmentationScores)==as.character(x))]],cnvrd2Object=results[[1]])
  polyMorphicResampling = list()
	new_samples = list()
	# do permutation 1000 times
	pdf('permutation_testing.pdf')
	real_polyMorphic = identifyPolymorphicRegion(Object = results[[1]],segmentObject=segment_scores[[1]])
	for ( i in 1:100){
		permuted_samples=lapply(sample_list,permuteSample,st = results[[1]]@st,results[[1]]@windows)
		new_samples[[i]]=list(segmentResults=(permuted_samples))
		polyMorphicResampling[[i]]=identifyPolymorphicRegion(Object=results[[1]],segmentObject=new_samples[[i]])
	}
	merged_regions = real_polyMorphic$subRegion
	# Get new subregions likely to just be the window size
	for ( i in 1:length(polyMorphicResampling)){	
	merged_regions=rbind(merged_regions,polyMorphicResampling[[i]]$subRegion)	
}
	sub_regions=sort(unique(c(merged_regions[,1],merged_regions[,2])))

# Hoang code stolen from his project
fix_subregions  <-  function (segment_data,polymorphic_region,sub_regions)	{
	print("Started Fixing Subregions")
	segmentResults <- segment_data$segmentResults
  sampleid <- sapply(segmentResults, function(x) x[1, 1])
  cna.out <- do.call(rbind, segmentResults)

  cna.out[, 4] <- cna.out[, 4] + results[[1]]@windows

  ##Obtain sub-regions#########
	subRegionData <- data.frame(sub_regions[-length(sub_regions)], sub_regions[-1])
	subRegionMatrix <- matrix(0, nrow = length(sampleid),
                          ncol = dim(subRegionData[1]))
		
	  for (ii in 1:length(sampleid)){
      subCNA <- cna.out[grep(sampleid[ii], cna.out[, 1]), c(3, 4, 6)]
			for (jj in 1:nrow(subCNA)){
				first_index = 1
				start_position = subCNA[jj,1]
				end_position = subCNA[jj,2]
				startSeq = findInterval(x=start_position,vec=subRegionData[,1])
				endSeq = findInterval(x=end_position,vec=subRegionData[,2])
				while(subRegionData[endSeq,2] <= end_position && endSeq <= length(subRegionData[endSeq,2])){
					endSeq = endSeq+1
				}
				indexes = seq(from=startSeq,to=endSeq)
        subRegionMatrix[ii, indexes] <-subCNA[jj,3] 
        }}
	polymorphic_region$subRegion=subRegionData
	polymorphic_region$subRegionMatrix=subRegionMatrix
	print(class(polymorphic_region))
	return(polymorphic_region)
}
polyMorphicResamplingNew=mapply(fix_subregions,new_samples,polyMorphicResampling,MoreArgs=list(sub_regions=sub_regions),SIMPLIFY=FALSE)

real_polyMorphic = fix_subregions(segment_scores[[1]],real_polyMorphic,sub_regions)


#calculate KS test for each sample against the null distribution if 95% of the time we get a 
# significant result we can be confident that we have a polymorphic copy number region.
#

# Ignore the start of the chromosome
 get_min_and_max = function(permutations){
		subRegions = permutations$subRegionMatrix
		min_and_max = apply(subRegions,2,function(x) return(c(min(x),max(x))))
		return(min_and_max)
}
 get_quantiles = function(permutations,probs){
		subRegions = permutations$subRegionMatrix
		return( apply(subRegions,2,function(x) return(quantile(x,probs=probs)[2] - quantile(x,probs=probs)[1])))
	}
	
	full_window_ks_test=lapply(polyMorphicResamplingNew,kstestwrapper,y=real_polyMorphic)
	ks.results = matrix(nrow=ncol(full_window_ks_test[[1]]),ncol=length(full_window_ks_test))
	anovaResluts = matrix(nrow=ncol(full_window_ks_test[[1]]),ncol=length(full_window_ks_test))
	for(j in 1:ncol(full_window_ks_test[[1]])){
	for(i in 1:length(full_window_ks_test)){
				ks.results[j,i] = as.numeric(full_window_ks_test[[i]][2,j])
				anovaResluts[j,i] = as.numeric(full_window_ks_test[[i]][2,j])
		}
	}
	ks.summary = apply(ks.results,1,function(x) sum(x<.05)/length(x))
	max_and_min = lapply(polyMorphicResamplingNew,get_min_and_max)
	max_and_min_mat=as.matrix(max_and_min)
	max_results_per_region=max_and_min[[1]][,2]
	min_results_per_region=max_and_min[[1]][,1]
	for(i in 2:length(max_and_min)){
			max_results_per_region=rbind(max_results_per_region,max_and_min[[i]][,2])
			min_results_per_region=rbind(min_results_per_region,max_and_min[[i]][,1])
	}
	max_results_per_region=max_and_min[[1]][,2]
	min_results_per_region=max_and_min[[1]][,1]
	

	permutedCI=function(permuted_row,alpha){
			permuted_row= sort(permuted_row)
			perc_ci = c()
			perc_ci[1] = permuted_row[length(permuted_row)*alpha/2]
			perc_ci[2] = permuted_row[length(permuted_row)*(1-alpha/2)]
			return(perc_ci)
	}
	perc_ci_max=apply(max_results_per_region,2,permutedCI,alpha=0.05)
	perc_ci_min=apply(min_results_per_region,2,permutedCI,alpha=0.05)
	min_breaks=apply(min_results_per_region,2,function(x) mean(x))
	max_breaks=apply(max_results_per_region,2,function(x) mean(x))	
	ci=c()
	ci[1] = min_breaks[floor(length(min_breaks))*(0.05/2)]
	ci[2] = max_breaks[floor(length(max_breaks))*(1-0.05/2)]
	abline(h=ci)	
	#construct ci.
		
	dev.off()
	# Stitch data back together
	# parse back together into something that looks like an actual format we need for identify polymorphic Regions
}

rangeDifference=function(x,y){
	#print('x')
	#print(length(x))
	#print('y')
	#print(length(y))
	return ((max(y) - min(y)) - (max(x)-min(x)))
}
rangeDifferenceQuantile=function(x,y){
	#print('x')
	#print(length(x))
	#print('y')
	#print(length(y))
	x_quants = quantile(x,probs=c(.95,.05))
	y_quants = quantile(y,probs=c(0.95,.05)) 
	return ((y_quants[1] - y_quants[2]) -(x_quants[1] - x_quants[2]) )
}
rangeDifferenceIQR=function(x,y){
	#print('x')
	#print(length(x))
	#print('y')
	#print(length(y))
	x_quants = quantile(x,probs=c(.75,.25))
	y_quants = quantile(y,probs=c(0.75,.25)) 
	return ((y_quants[1] - y_quants[2]) -(x_quants[1] - x_quants[2]) )
}
min_and_max_real_data = get_min_and_max(real_polyMorphic)
ranges=lapply(max_and_min,function(x) return(mapply(rangeDifference,split(x,row(x)),split(min_and_max_real_data,col(min_and_max_real_data)))))
samples = rownames(segment_scores[[1]]$segmentationScores)
permuteSegmentationScores(segment_scores[[1]],results[[1]])
#bartlett test
ranges_per_region = as.matrix(ranges)
ranges_per_region = apply(q_ranges_per_region,1,function(x) {return(x[[1]])})
#quantile tests

q_permutationdata = lapply(polyMorphicResamplingNew,get_quantiles,probs=c(.05,.95))
q_realdata = get_quantiles(real_polyMorphic,probs=c(.05,.95)) 

p_values  = c()
for ( i in 1:nrow(q_ranges_per_region)){
		pvalue[i] = sum(q_ranges_per_region[i,] < q_realdata[i])/length(q_realdata[i])
		if(pvalue[i] > 95){
			pvalue[i] = 1
		}
		else{
			pvalue[i] = 0
		}
}
wilcox.test
#range_quantile=lapply(q_permutationdata,function(x) return(mapply(rangeDifference,split(x,row(x)),split(q_realdata,col(q_realdata)))))
#range_quantile=lapply(q_permutationdata,function(x) return(mapply(rangeDifferenceQuantile,split(x,row(x)),split(min_and_max_real_data,col(min_and_max_real_data)))))
q_ranges_per_region = as.matrix(q_permutationdata)
q_ranges_per_region = apply(q_ranges_per_region,1,function(x) {return(x[[1]])})
quantile_per_region = apply(q_ranges_per_region,1,function(x) return(quantile(x,probs=c(.05,.95))))
dif_in_quantile = apply(quantile_per_region,2,function(x) return(x[2] - x[1]))
thousand_row_difference=c()
counts = c()
poly_morphic_regions = function(ranges_per_region,distance_in_std_deviation)
for ( i in 1:nrow(ranges_per_region)){
	print(i)
	counts[i] = 0
	for( j in 1:ncol(ranges_per_region)){
#				print(ranges_per_region[1000,j])
				if((mean(ranges_per_region) + distance_in_std_deviation* sd(ranges_per_region[,j])) <= ranges_per_region[i,j]){
					counts[i] = counts[i] + 1		
				}
		}
}
count2 = c()
for ( i in 1:nrow(ranges_per_region)){
	print(i)
	count2[i] = 0
	for( j in 1:ncol(ranges_per_region)){
#				print(ranges_per_region[1000,j])
				if((mean(ranges_per_region) + 1.5* sd(ranges_per_region[,j])) <= ranges_per_region[i,j]){
		count2[i] = counts[i] + 1		
		}
	}
}
count2 = c()

new_counts = c()
count_indices = c()
index = 1
for( i in 1:length(counts)){
	if( counts[i] > 95){
		count_indices[index] = i
		index = index + 1
		print(i)
		new_counts[i] = 1
	}else{
		new_counts[i] = 0
	}
}

#range difference

#cut out windows with less than 2 in length


truncate_polymorphic_regions = function(new_counts){
	start_index = c()
	end_index = c()
	index = 1
	set_start = T
	window_size=0
	for ( i in 1:length(new_counts)){
		if(new_counts[i] == 1){
				if(set_start == T){
					start_index[index] = (i * 1000) - 1000	
					set_start = F
					window_size=0
				}
			window_size = window_size +  1
		}else{
				# going to be one longer than it needs to be if 
				# what we are testing really doesn't match up
				if(set_start == F && window_size >= 5){
					end_index[index] = (i * 1000) 
					set_start =T 
					index= index + 1
				}else{
					set_start = T
				}
		}
	}
	if(length(start_index) == (length(end_index) + 1)){
		start_index=start_index[1:length(start_index)- 1]
	}
	print(length(start_index))
	print(length(end_index))
	return(cbind(start_index,end_index))
}
ninetyseven=quantile(apply(q_ranges_per_region,1,quantile,probs=.975),probs=.975)
putative_regions = truncate_polymorphic_regions(sapply(q_realdata,function(x) ifelse(x >= ninetyseven,1,0)))
putative_regions = truncate_polymorphic_regions(new_counts)
putative_regions = truncate_polymorphic_regions(pvalue)
temp_x = identifyPolymorphicRegion(Object = results[[1]],segmentObject=segment_scores[[1]],xlim=limits)
temp_y = identifyPolymorphicRegion(Object = results[[1]],segmentObject=new_samples[[2]],xlim=limits)
pdf('putative_regions.pdf')
for ( i in 1:10){
limits=c(putative_regions[i,1] ,putative_regions[i,2] )
plotPolymorphicRegion(Object = results[[1]],polymorphicRegionObject=x,xlim=limits)
plotPolymorphicRegion(Object = results[[1]],polymorphicRegionObject=y,xlim=limits)
}
dev.off()
