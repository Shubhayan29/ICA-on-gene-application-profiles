library(fastICA)
library(foreach)
library(doParallel)

#THE BELOW FUNCTION WILL GIVE IC's FOR ANY DATASET.
ICs=function(X,N_1,N_2,r){
  if (N_1 + N_2 < ncol(X) |  N_1 + N_2 == ncol(X)){
    if (N_1 + N_2 < ncol(X)){
      
      #PCA
      PCs = prcomp(t(X), scale = TRUE)
      
      Combined_indep_comps = matrix(, nrow = nrow(X), ncol = 0)
      
      cores <- detectCores()
      cl <- makeCluster(cores[1]-1) #not to overload your computer
      registerDoParallel(cl)
      
      Combined_indep_comps = foreach(x = 1:r, .combine=cbind,.packages = c("fastICA")) %dopar%{
        
        #sampling
        T=sample(c((N_1 + 1):ncol(X)),N_2,replace=FALSE)
        
        #Applying ICA on the combined PCs
        Combined_PCs = cbind(PCs$rotation[, 1:N_1], PCs$rotation[, T])
        
        Indep_components = fastICA(
          Combined_PCs,
          ncol(Combined_PCs),
          alg.typ = "parallel",
          fun = "logcosh",
          alpha = 1,
          method = "C",
          row.norm = FALSE,
          maxit = 200,
          tol = 0.0001,
          verbose = FALSE
        )
        
        Indep_components$S
      }
      
      #stop cluster
      stopCluster(cl)
    }
    else if (N_1 + N_2 == ncol(X)){
      #PCA
      PCs = prcomp(t(X), scale = TRUE)
      
      cores <- detectCores()
      cl <- makeCluster(cores[1]-1) #not to overload your computer
      registerDoParallel(cl)
      
      Combined_indep_comps = foreach(x = 1:r, .combine=cbind,.packages = c("fastICA")) %dopar%  {
        Combined_PCs = PCs$rotation
        Indep_components = fastICA(
          Combined_PCs,
          ncol(Combined_PCs),
          alg.typ = "parallel",
          fun = "logcosh",
          alpha = 1,
          method = "C",
          row.norm = FALSE,
          maxit = 200,
          tol = 0.0001,
          verbose = FALSE
        )
        
        Indep_components$S
      }
      
      #stop cluster
      stopCluster(cl)
    }
    
    return(Combined_indep_comps)
  }
}


#NOW WE WILL GENERATE THE FIRST DATASET.

#install packages required (once per computer)
install.packages("graphsim")

#load required packages (once per R instance)
library("graphsim")

#load packages for examples
library("igraph"); library("gplots"); library("scales")

#generate graph structure and inhibitions
edges=c()
state=c()

for (i in 0:4500){
  edges=rbind(edges,c(5*i+1,5*i+2),c(5*i+2,5*i+3),c(5*i+3,5*i+4),c(5*i+4,5*i+5))
  edges=rbind(edges,c(5*i+5,5*i+1))
  edges=rbind(edges,c(5*i+1,5*i+6),c(5*i+2,5*i+7),c(5*i+3,5*i+8),c(5*i+4,5*i+9))
  edges=rbind(edges,c(5*i+5,5*i+10))
  state=c(state,c(1,1,-1,1,-1,1,1,-1,1,-1))
}

edges=matrix(as.character(edges),nrow(edges),ncol(edges))
graph=graph.edgelist(edges, directed = TRUE)

#plot graph structure with inhibitions
plot_directed(graph, state=state, layout = layout.kamada.kawai,
              cex.node = 2, cex.arrow = 4, arrow_clip = 0.2)

#simulate data
expr_1=generate_expression(10000, graph, cor = 0.8, mean = 0,
                           comm = FALSE, dist =TRUE, absolute = FALSE, state = state)

p_1=ncol(expr_1)

N1_1=c(5000,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500,7750,8000)
N2_1=c(250,500,1000,1250,1500,1750,2000,2250,2500)
RUN_1=c(100,250,500,750,900,1000,2000,2500,3000,4500,5000)

#List that will contain all the independent components of the first dataset.
IC_lst1=c()

for (i in 21:45){
  for (N_1 in N1_1){
    for (run in RUN_1){
      for (N_2 in N2_1){
        if (N_2>=floor((p_1-N_1)/run)){
          IC_lst1=c(IC_lst1,list(ndatapoints=(500*i),runs=run,N1=N_1,N2=N_2,ICs(expr_1[1:(500*i),],N_1,N_2,run)))
        }
      }
    }
  }
}



#NOW WE WILL GENERATE THE SECOND DATASET.

install.packages("SeqNet")
library("SeqNet")

#number of genes
n_genes=20000

#number of persons
n_persons=8000

#generate the underlying network structure over which the RNA sequence dataset will be generated.
network=random_network(n_genes)

#generate data
generated_data=gen_rnaseq(n_persons, network)

#Simulated data
expr_2=t(generated_data$x)

p_2=ncol(expr_2)

N1_2=c(4000,4250,4500,4750,5000,5250,5500,5750,6000,6250,6500)
N2_2=c(250,500,1000,1250,1500,2000,2500)
RUN_2=c(100,250,500,750,900,1000,2000,2500,3000,4500,5000)

#List that will contain all the independent components of the second dataset.
IC_lst2=c()

for (i in 10:20){
  for (N_1 in N1_2){
    for (run in RUN_2){
      for (N_2 in N2_2){
        if (N_2>=floor((p_2-N_1)/run)){
          IC_lst2=c(IC_lst2,list(ndatapoints=(1000*i),runs=run,N1=N_1,N2=N_2,ICs(expr_2[1:(1000*i),],N_1,N_2,run)))
        }
      }
    }
  }
}



