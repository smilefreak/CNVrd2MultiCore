library(Rmpi)
library(CNVrd2)
args  <- commandArgs(T)
load('/home/james.boocock/disease_genes_project/.RData')
mpi.spawn.Rslaves(nslaves=args[1])

   # In case R exits unexpectedly, have it automatically clean up
    # resources taken up by Rmpi (slaves, memory, etc...)

permuteChromosome <- function(task){
    chr  <- task[2]
    i <- task[1]
    st  <- task[3]
    windows  <- task[4]
    sample_list  <- sample_list_per_chromosome[[chr]]
    result <- lapply(sample_list,permuteSample,st=st,windows=windows)
    result <- list(segmentResults=result)
    polyMorphicResampling=identifyPolymorphicRegion(Object=cnvrd2_objects[[chr]],segmentObject=result,plotPolymorphicRegion=F)
    return(list(i=i,chr=chr,result=polyMorphicResampling))
}

names(cnvrd2_objects)  <- lapply(cnvrd2_objects,function(x){ x@chr})

run_slave_permutation <- function(){
    require(CNVrd2)
    task  <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    task_info  <- mpi.get.sourcetag()
    tag <- task_info[2]
    while(tag != 2){
       result <- permuteChromosome(task) 
       mpi.send.Robj(result,0,1)
       task  <- mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
       task_info  <- mpi.get.sourcetag()
       tag <- task_info[2]
    }
    junk <- 0
    mpi.send.Robj(junk,0,2)
}

#collect results

#`sourceCpp("generate_permutations.cpp")

permuteSample <- function(segmentationResultsForSample,st,windows){
        ranks=seq(1,nrow(segmentationResultsForSample))
        reorder_ranks=sample(ranks,replace=F)
        
        #segmentationResultsForSample = segmentationResultsForSample[reorder_ranks,]
        # Improve speed of this horrible for loop lol
        for( i in 1:nrow(segmentationResultsForSample)){
            dif = segmentationResultsForSample[i,4] - segmentationResultsForSample[i,3] + windows
            segmentationResultsForSample[i,3] = st
            segmentationResultsForSample[i,4] = st  + dif - windows
            st = st + dif
        }
        return(segmentationResultsForSample)
}
samples=rownames(segment_scores[[1]]$segmentationScores)
sample_list_per_chromosome  <- list()
for( i in 1:length(segment_scores)){
    sample_list_per_chromosome[[i]]=lapply(samples, function(x) segment_scores[[i]]$segmentResults[[which(rownames(segment_scores[[i]]$segmentationScores)==as.character(x))]])
}

mpi.bcast.Robj2slave(permuteSample)
#mpi.bcast.Robj2slave(segment_scores)
mpi.bcast.Robj2slave(permuteChromosome)
mpi.bcast.Robj2slave(run_slave_permutation)
mpi.bcast.Robj2slave(sample_list_per_chromosome)
mpi.bcast.Robj2slave(cnvrd2_objects)
permutations <- 10000
tasks <- list()
for(j in 1:length(sample_list_per_chromosome)){
    start = ((j*permutations) %/% permutations -1) * permutations
    for(i in 1:permutations){
        tasks[[start + i]] <- c(start + i,which(cnvrd2_objects[[j]]@chr==names(cnvrd2_objects)),cnvrd2_objects[[j]]@st,cnvrd2_objects[[j]]@windows)
    }
}

# Temp shorten task list
#tasks  <- tasks[1:5]
print(length(tasks))
mpi.bcast.cmd(run_slave_permutation())
n_slaves  <- mpi.comm.size()-1
task_assignees <- rep(1:n_slaves, length=length(tasks))
for (i in 1:length(tasks)){
    print(length(tasks))
    slave_id <- task_assignees[i]
    mpi.send.Robj(tasks[[i]],slave_id,1)    
}
result_list <- list()
j= 1
file_count = 1
for(i in 1:length(tasks)){
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
for(i in 1:n_slaves){
    junk <- 0
    mpi.send.Robj(junk,i,2)
}
for(i in 1:n_slaves){
    mpi.recv.Robj(mpi.any.source(),i,2)
}
mpi.close.Rslaves()
mpi.quit(save='no')

