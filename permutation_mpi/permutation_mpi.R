library(Rmpi)
#library(CNVrd2)
library(Rcpp)
args  <- commandArgs(T)
load('/home/james.boocock/disease_genes_project/apple_segment.RData')
load('/home/james.boocock/disease_genes_project/cnvrd2_objects_nesi.Rdata')
write.table(date(),file="start.txt")
mpi.spawn.Rslaves(nslaves=200)
   # In case R exits unexpectedly, have it automatically clean up
    # resources taken up by Rmpi (slaves, memory, etc...)
permuteChromosome2 <- function(task){
    chr  <- task[2]
    i <- task[1]
    st  <- task[3]
    en = task[4]
    windows  <- task[5]
    sample_list  <- sample_list_per_chromosome[[chr]]
    result <- lapply(sample_list,permuteSample,st=st,windows=windows)
    result <- list(segmentResults=result)
    polyMorphicResampling = getSubRegionMatrixFromSegmentScores(segmentResults=result$segmentResults,window_size = windows,st,en)
    return(list(i=i,chr=chr,result=polyMorphicResampling))
}
permuteChromosome <- function(task){
    chr  <- task[2]
    i <- task[1]
    st  <- task[3]
    en  <- task[4]
    windows  <- task[5]
    sample_list  <- sample_list_per_chromosome[[chr]]
    result <- lapply(sample_list,permuteSample,st=st,windows=windows)
    result <- list(segmentResults=result)
    polyMorphicResampling=identifyPolymorphicRegion(Object=cnvrd2_objects[[chr]],segmentObject=result,plotPolymorphicRegion=F)
    return(list(i=i,chr=chr,result=polyMorphicResampling))
}

names(cnvrd2_objects)  <- lapply(cnvrd2_objects,function(x){ x@chr})
require(generateSubRegions)

run_slave_permutation <- function(){
    require(CNVrd2)
    require(generateSubRegions)
#    sourceCpp("/home/james.boocock/CNVrd2MultiCore/Cppextensions/polyMorphicStretch.cpp")
    task  <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    task_info  <- mpi.get.sourcetag()
    tag <- task_info[2]
    while(tag != 2){
       result <- permuteChromosome2(task)
       mpi.send.Robj(result,0,1)
       task  <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
       task_info  <- mpi.get.sourcetag()
       tag <- task_info[2]
    }
    junk <- 0
    print("Finished processing DATA")
    mpi.send.Robj(junk,0,2)
}

#collect results

#`sourceCpp("generate_permutations.cpp")

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
samples=rownames(segment_results[[1]]$segmentationScores)
sample_list_per_chromosome  <- list()
chromosomes = c(10:17,1:9)
for( i in 1:length(segment_results)){
    sample_list_per_chromosome[[i]]=lapply(samples, function(x) segment_results[[i]]$segmentResults[[which(rownames(segment_results[[i]]$segmentationScores)==as.character(x))]])
}

mpi.bcast.Robj2slave(permuteSample)
#mpi.bcast.Robj2slave(segment_results)
#mpi.bcast.Robj2slave(permuteChromosome)
mpi.bcast.Robj2slave(permuteChromosome2)
mpi.bcast.Robj2slave(run_slave_permutation)
mpi.bcast.Robj2slave(sample_list_per_chromosome)
mpi.bcast.Robj2slave(cnvrd2_objects)
permutations <- 10000
tasks <- list()
for(j in 1:length(sample_list_per_chromosome)){
    start = ((j*permutations) %/% permutations -1) * permutations
    for(i in 1:permutations){
        tasks[[start + i]] <- c(start + i,which(cnvrd2_objects[[j]]@chr==names(cnvrd2_objects)),cnvrd2_objects[[j]]@st,cnvrd2_objects[[j]]@en,cnvrd2_objects[[j]]@windows)
    }
}

# Temp shorten task list
#tasks  <- tasks[1:5]
mpi.bcast.cmd(run_slave_permutation())
n_slaves  <- mpi.comm.size()-1
task_assignees <- rep(1:n_slaves, length=length(tasks))
for (i in 1:length(tasks)){
    slave_id <- task_assignees[i]
    mpi.send.Robj(tasks[[i]],slave_id,1)    
}

#mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
result_list <- list()
j= 1
file_count = 1
for(i in 1:length(tasks)){
    print(i)
    messages  <-  mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    result_list[[j]]  <- messages
    if ( j > 100){
        save(result_list,file=paste0('results',file_count,'.Rdata'))
        j = 1
        file_count = file_count + 1
    }
    result_list[[j]]  <- messages
    j = j + 1
}
save(result_list,file=paste0('results',file_count,'.Rdata'))
for(i in 1:n_slaves){
    junk <- 0
    mpi.send.Robj(junk,i,2)
}
for(i in 1:n_slaves){
    mpi.recv.Robj(mpi.any.source(),2)
}
mpi.close.Rslaves()
write.table(date(),file='finish_time.txt')
mpi.quit(save='no')

