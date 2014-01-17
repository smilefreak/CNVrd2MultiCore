library(Rmpi)
library(CNVrd2)
args  <- commandArgs(T)
load('/home/james.boocock/disease_genes_project/.RData')
mpi.spawn.Rslaves(nslaves=as.numeric(args[1]))


permuteChromosome <- function(task){
    chr  <- task[1]
    i <- task[2]
    result <- permuteSample(segment_scores[[chr]]$segmentationScores)
    return(list(i=i,chr=chr,result=result))
}


run_slave_permutation <- function(){
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
samples=rownames(segment_scores[[1]]$segmentationScores)
sample_list_per_chromosome  <- list()
for( i in 1:length(segment_scores)){
    sample_list_per_chromosome[[i]]=lapply(samples, function(x) segment_scores[[i]]$segmentResults[[which(rownames(segment_scores[[i]]$segmentationScores)==as.character(x))]])
}

mpi.bcast.Robj2slave(permuteSample)
mpi.bcast.Robj2slave(segment_scores)
mpi.bcast.Robj2slave(permuteChromosome)

tasks <- list()
for(i in 1:5){
    for(j in 1:length(sample_list_per_chromosome)){
        tasks[[i*j]] <- c(i,cnvrd2_objects[[j]]@chr)
    }
}

n_slaves  <- mpi.comm.size()-1
task_assignees <- rep(1:n_slaves, length=length(tasks))
for (i in 1:length(tasks)){
    slave_id <- task_assignees[i]
    mpi.send.Robj(tasks[[i]],slave_id,1)    
}

result_list <- list()
for( i in 1:length(tasks)){
    messages  <-  mpi.recv.Robj(mpi.any.source(),mpi.any.tag())
    result_list[[i]]  <- messages
    save(result_list[[i]],paste0('results',i,'.Rdata'))
}
save(result_list,'result_list.Rdata')
for(i in 1:n_slaves){
    junk <- 0
    mpi.send.Robj(junk,i,2)
}
for(i in 1:n_slaves){
    mpi.recv.Robj(mpi.any.source(),i,2)
}
mpi.close.Rslaves()
mpi.quit(save='no')

