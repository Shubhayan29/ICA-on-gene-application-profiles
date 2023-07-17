library(fastICA)
library(foreach)
library(doParallel)

#THE BELOW FUNCYTION WILL GIVE IC's FOR ANY DATASET.
ICs=function(X,Pi1,Pi2,r){
  if (Pi1 + Pi2 < 1 |  Pi1 + Pi2 == 1){
    if (Pi1 + Pi2 < 1) {
      #PCA
      PCs = prcomp(t(X), scale = TRUE)
      #print(PCs)
      #Cumulative explained variance vector
      Total_explained = summary(PCs)$importance[3, ]
      
      #Proportion of Variance Vector
      Proportion_of_Variance = summary(PCs)$importance[2, ]
      
      #Number of PCs requuired for explaining Pi1 proportion
      N_1 = min(which(Total_explained > Pi1 |
                        Total_explained == Pi1))
      #PCs$rotation[,1:N_1]
      
      if ((1 - Total_explained[N_1]) < Pi2) {
        Combined_PCs = PCs$rotation[, 1:N_1]
        Combined_indep_comps = matrix(, nrow = nrow(X), ncol = 0)
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
        
        Combined_indep_comps = cbind(Combined_indep_comps, Indep_components$S)
      }else{
        #Number of PCs required for explaining Pi2 proportion
        Combined_indep_comps = matrix(, nrow = nrow(X), ncol = 0)
        
        cores <- detectCores()
        cl <- makeCluster(cores[1]-1) #not to overload your computer
        registerDoParallel(cl)
        
        N_tot=min(which(Total_explained > Pi1+Pi2 |
                          Total_explained == Pi1+Pi2))
        
        N_2=N_tot-N_1
        
        Combined_indep_comps = foreach(x = 1:r, .combine=cbind,.packages = c("fastICA")) %dopar%  {
          
          T=sample(c((N_1 +1):ncol(PCs$rotation)),N_2,replace=FALSE)
          
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
    }else if (Pi1 + Pi2 == 1) {
      #PCA
      PCs = prcomp(t(X), scale = TRUE)
      
      #Cumulative explained variance vector
      Total_explained = summary(PCs)$importance[3, ]
      
      #Proportion of Variance Vector
      Proportion_of_Variance = summary(PCs)$importance[2, ]
      
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
      stopCluster(cl)}
    
  }
  return(Combined_indep_comps)
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

for (i in 0:3000){
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



#NOW WE WILL GET THE IC's ONE BY ONE FROM EACH DATASET.

RUN=c(100,250,500,750,900,1000,2000,2500,3000,4500,5000,6000,7000,7500,9000,9500,10000)
P1=c(0.75,0.8,0.85,0.9)
P2=c(0.05,0.08,0.1,0.12,0.15)

#List that will contain all the independent components of the first dataset.
#IC_lst1=c()


for (Pi1 in P1){
  for (Pi2 in P2){
    for (i in 1:30){
      for (r in RUN){
        if (Pi1+Pi2<1 || Pi1+Pi2==1){
          
          #IC_lst1=c(IC_lst1,list(ndatapoints=(500*i),runs=r,PI1=Pi1,PI2=Pi2,ICs(expr_1[1:(500*i),],Pi1,Pi2,r)))
          save(list(ndatapoints=(500*i),runs=r,PI1=Pi1,PI2=Pi2,ICs(expr_1[1:(500*i),],Pi1,Pi2,r)),file=file_name1,ascii = FALSE,oldstyle=FALSE)
        }
      }
    }
  }
}

#List that will contain all the independent components of the second dataset.
#IC_lst2=c()


for (Pi1 in P1){
  for (Pi2 in P2){
    for (i in 1:20){
      for (r in RUN){
        if (Pi1+Pi2<1 || Pi1+Pi2==1){
          
          #IC_lst2=c(IC_lst2,list(ndatapoints=(1000*i),runs=r,PI1=Pi1,PI2=Pi2,ICs(expr_2[1:(1000*i),],Pi1,Pi2,r)))
          save(list(ndatapoints=(1000*i),runs=r,PI1=Pi1,PI2=Pi2,ICs(expr_2[1:(1000*i),],Pi1,Pi2,r)),file=file_name2,ascii=FALSE,oldstyle=FALSE)
        }
      }
    }
  }
}

