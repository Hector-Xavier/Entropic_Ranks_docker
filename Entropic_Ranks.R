arguments <- commandArgs(TRUE)
data_under_analysis <- read.table(arguments[1],sep="\t",dec=".",header=TRUE)
population_vector <- read.table(arguments[2],sep="\t",dec=".",header=TRUE)[,1]
if (arguments[3]!="null")
{data_origin <- read.table(arguments[3],sep="\t",dec=".",header=TRUE)[,1]} else {data_origin <- NULL}
granularity <- as.integer(arguments[4])
supervised <- as.logical(arguments[5])
process_log <- as.logical(arguments[6])
export_plots <- as.logical(arguments[7])
create_output_files <- as.logical(arguments[8])
is_logged <- as.logical(arguments[9])
logbase <- as.integer(arguments[10])
huge_feature_list <- as.logical(arguments[11])

library("RankProd")
library("entropy")
library("factoextra")
library("plotly")

setwd("/data")

entropic_ranks <- function(data_under_analysis,population_vector,data_origin=NULL,granularity=1,supervised,process_log=FALSE,export_plots=FALSE,create_output_files=FALSE,is_logged,logbase,huge_feature_list=FALSE)
{
  if (is.null(data_origin))
    data_origin <- rep(1,length(population_vector))
  
  message("Calculating Rank Products. May take a long time, depending on data set size.")
  comparison <- RPadvance(data_under_analysis,cl=population_vector,origin=data_origin,logged=is_logged,na.rm=FALSE,gene.names=rownames(data_under_analysis),plot=process_log,huge=TRUE)
  if (huge_feature_list){
    message("Investigating only the first 20000 features.")
    rank_product_lists <- topGene(comparison,num.gene=20000,logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
  } else {
    rank_product_lists <- topGene(comparison,cutoff=0.99,method="pfp",logged=is_logged,logbase=logbase,gene.names=rownames(data_under_analysis))
  }
  if (is.null(rownames(rank_product_lists$Table1)) || is.null(rownames(rank_product_lists$Table2)))
    stop("Rank Products is unable to adequately identify differentially expressed features.")
  
  path_down <- NULL
  path_up <- NULL
  if (export_plots)
  {
    if(!file.exists(paste(getwd(),"Entropic Ranks plots",sep="/")))
      dir.create(paste(getwd(),"Entropic Ranks plots",sep="/"))
    path_down <- paste(getwd(),"Entropic Ranks plots","Downregulated",sep="/")
    path_up <- paste(getwd(),"Entropic Ranks plots","Upregulated",sep="/")
  }
  
  if (supervised)
  {
    message("Calculating suggested cutoff points for the downregulated feature list...")
  } else {
    message("Trimming down the list of downregulated features in unsupervised mode...")
  }
  downregulated_suggestions <- isolate_significant_elements(rank_product_lists$Table1[,2],granularity,supervised,process_log,export_plots,path=path_down)
  if (supervised)
  {
    message("Calculating suggested cutoff points for the upregulated feature list...")
  } else {
    message("Trimming down the list of upregulated features in unsupervised mode...")
  }
  upregulated_suggestions <- isolate_significant_elements(rank_product_lists$Table2[,2],granularity,supervised,process_log,export_plots,path=path_up)
  
  if (supervised || export_plots || process_log)
  {
    downregulated_suggestions$entropy_classification_topography <- downregulated_suggestions$entropy_classification_topography %>% layout(title="Topography of entropy clustering in downregulation")
    upregulated_suggestions$entropy_classification_topography <- upregulated_suggestions$entropy_classification_topography %>% layout(title="Topography of entropy clustering in upregulation")
  }
  
  if (!supervised && create_output_files)
  {
    if (!is.null(rownames(rank_product_lists$Table1)))
      write.table(file="Downregulated list [original].txt",rank_product_lists$Table1[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
    if (!is.null(rownames(rank_product_lists$Table2)))
      write.table(file="Upregulated list [original].txt",rank_product_lists$Table2[,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
  }
  if (!supervised)
  {
    if (downregulated_suggestions$suggested_cutoff_surface==1)
    {
      single_feature <- as.table(t(rank_product_lists$Table1[1:1,]))
      rownames(single_feature) <- rownames(rank_product_lists$Table1)[1]
      rank_product_lists$Table1 <- single_feature
      rm(single_feature)
    } else {
      rank_product_lists$Table1 <- rank_product_lists$Table1[1:downregulated_suggestions$suggested_cutoff_surface,]
    }
    if (upregulated_suggestions$suggested_cutoff_surface==1)
    {
      single_feature <- as.table(t(rank_product_lists$Table2[1:1,]))
      rownames(single_feature) <- rownames(rank_product_lists$Table2)[1]
      rank_product_lists$Table2 <- single_feature
      rm(single_feature)
    } else {
      rank_product_lists$Table2 <- rank_product_lists$Table2[1:upregulated_suggestions$suggested_cutoff_surface,]
    }
  }
  if (!supervised && create_output_files)
  {
    if (downregulated_suggestions$suggested_cutoff_surface > 0)
    {
      if (downregulated_suggestions$suggested_cutoff_surface > 1)
      {
        write.table(file="Downregulated list [information-dense].txt",rank_product_lists$Table1[1:downregulated_suggestions$suggested_cutoff_surface,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
      } else {
        single_feature <- as.table(t(rank_product_lists$Table1[,3:5]))
        rownames(single_feature) <- rownames(rank_product_lists$Table1)
        write.table(file="Downregulated list [information-dense].txt",single_feature,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
        rm(single_feature)
      }
    }
      
    if (upregulated_suggestions$suggested_cutoff_surface > 0)
    {
      if (upregulated_suggestions$suggested_cutoff_surface > 1)
      {
        write.table(file="Upregulated list [information-dense].txt",rank_product_lists$Table2[1:upregulated_suggestions$suggested_cutoff_surface,3:5],sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
      } else {
        single_features <- as.table(t(rank_product_lists$Table2[,3:5]))
        rownames(single_feature) <- rownames(rank_product_lists$Table2)
        write.table(file="Upregulated list [information-dense].txt",single_feature,sep="\t",quote=FALSE,row.names=TRUE,col.names=TRUE)
        rm(single_feature)
      }
    }
    message("Output files created successfully.")
  }
  
  return(list(ranked_lists=rank_product_lists,downregulated_surface_suggestions=downregulated_suggestions$suggested_cutoff_surface,downregulated_classification_topography=downregulated_suggestions$entropy_classification_topography,upregulated_surface_suggestions=upregulated_suggestions$suggested_cutoff_surface,upregulated_classification_topography=upregulated_suggestions$entropy_classification_topography))
}

isolate_significant_elements <- function(ordered_vector,granularity=1,supervised=FALSE,process_log=FALSE,export_plots=FALSE,path=NULL)
{
  bin_min <- 15
  bin_max <- 35
  bin_increment <- 5
  window_min <- 60
  window_max <- 250
  window_increment <- 10
   
  if (process_log)
  {
    par(mfcol=c(4,2))
  } else {
    progress <- txtProgressBar(max=((bin_max-bin_min)/bin_increment +1) * ((window_max-window_min)/window_increment +1),char="=",style=3)
  }
  suggested_cutoffs <- c()
  entropy_chartography <- matrix(nrow=(bin_max-bin_min+bin_increment)/bin_increment*(window_max-window_min+window_increment)/window_increment,ncol=3+length(ordered_vector))
  counter <- 0
  for (i in seq(from=bin_min, to=bin_max, by=bin_increment))
    for (j in seq(from=window_min, to=window_max, by=window_increment))
    {
      if (process_log)
      {
        cat("Calculating using ",i," bins and a sliding window of ",j," features","\n",sep="")
      } else {
        setTxtProgressBar(progress,(i-bin_min)/bin_increment*((window_max-window_min)/window_increment + 1) + (j-window_min)/window_increment + 1)
      }
      temp_entropy <- entropic_analysis(ordered_vector,step_up=granularity,bins=i,window_size=j,verbose=process_log,export_plots,path)
      suggested_cutoffs <- c(suggested_cutoffs,temp_entropy$cutoff_surface)
      counter <- counter + 1
      entropy_chartography[counter,1] <- i
      entropy_chartography[counter,2] <- j
      entropy_chartography[counter,3:(length(temp_entropy$entropy_classification)+2)] <- temp_entropy$entropy_classification
    }
  if (process_log)
  {
    par(mfrow=c(1,1))
  } else {
    close(progress)
  }
  
  if (supervised || export_plots || process_log)
  {
    if (supervised || process_log)
    {
      print(table(suggested_cutoffs))
      message("Most consistent cutoff point: feature no ",as.integer(rownames(table(suggested_cutoffs))[table(suggested_cutoffs) == max(table(suggested_cutoffs))])[1],".")
    }
    entropy_cluster <- entropy_chartography[sort(entropy_chartography[,2],index.return=TRUE)$ix,c(3,3:1003)]
    entropy_topography <- plot_ly(z = ~(entropy_cluster))
    entropy_topography <- entropy_topography %>% add_surface()
    entropy_topography <- entropy_topography %>% layout(title="Topography of entropy clustering",scene=list(xaxis=list(title="Feature index"),yaxis=list(title="Window size (% of maximum)"),zaxis=list(title="Entropy cluster (high vs low)"),camera=list(eye=list(x=-1.30,y=-1.4,z=1)),aspectratio=list(x=1, y=1, z=0.7)))
  } else {
  entropy_topography <- NULL
  }
  
  if (!supervised)
  {
    return(list(suggested_cutoff_surface=as.integer(rownames(table(suggested_cutoffs)))[table(suggested_cutoffs) == max(table(suggested_cutoffs))][1],entropy_classification_topography=entropy_topography))
  } else {
    return(list(suggested_cutoff_surface=suggested_cutoffs,entropy_classification_topography=entropy_topography))
  }
}

entropic_analysis <- function(ordered_vector,step_up=1,window_size,bins,verbose=FALSE,export_plots=FALSE,path=NULL)
{
  if (export_plots)
  {
    if (is.null(path))
      path <- paste(getwd(),"Entropic Ranks plots",sep="/")
    if (!file.exists(path))
      dir.create(path)
  }
  
  differences <- ordered_vector[seq(2,length(ordered_vector))]-ordered_vector[seq(1,length(ordered_vector)-1)]
  entropy_plotter <- vector(length=floor((length(differences)-window_size)/step_up))
  for (i in 0:(length(entropy_plotter)-1))
    entropy_plotter[i+1] <- entropy(discretize(differences[(i*step_up+1):(i*step_up+window_size)],numBins=bins),method="Laplace")
    
  entropy_clusters <- eclust(entropy_plotter, "kmeans", k=2, nstart=200, graph=FALSE)
  
  if (verbose)
  {
    cat("Calculating entropies of ",length(entropy_plotter)," ovelapping windows","\n","Suggested cutoff at feature no ",seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1,", at a mean entropy of ",mean(entropy_plotter[1:seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1]),"\n","Last 1/3 minimum entropy: ",min(entropy_plotter[floor(length(entropy_plotter)*2/3):length(entropy_plotter)]),"\n",sep="")
    barplot(entropy_plotter-4/5*min(entropy_plotter),offset=4/5*min(entropy_plotter),border=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],col=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],xlab="Feature index",ylab="Entropy (nats)",main="Entropy value clustering",names.arg=seq(length(entropy_plotter))*step_up)
    legend("bottomright",legend=c(paste("Window size:",window_size,sep=" "),paste("Number of bins:",bins,sep=" "),paste("Granularity: ",step_up,sep=" ")))
    barplot(entropy_clusters$silinfo$widths$sil_width,border=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],col=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
    legend("topright",legend=c("High-entropy cluster","Low-entropy cluster"),pch=15,col=c("gold1","dodgerblue3"))
  }
  
  if (export_plots)
  {
    png(file=file.path(path,paste("Entropy distribution [bins_",bins," - window size_",window_size,"].png", sep = "")),width=1000,height=1000)
    par(mfcol=c(2,1))
    barplot(entropy_plotter-4/5*min(entropy_plotter),offset=4/5*min(entropy_plotter),border=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],col=c("gold1","dodgerblue3")[as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1)],xlab="Feature index",ylab="Entropy (nats)",main="Entropy value clustering",names.arg=seq(length(entropy_plotter))*step_up)
    legend("bottomright",legend=c(paste("Window size:",window_size,sep=" "),paste("Number of bins:",bins,sep=" "),paste("Granularity: ",step_up,sep=" ")))
    barplot(entropy_clusters$silinfo$widths$sil_width,border=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],col=c("gold1","dodgerblue3")[as.integer(entropy_clusters$silinfo$widths$cluster==(as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1))+1],ylim=c(-0.3,1),ylab="Silhouette width",xlab="Clustered elements",main="K-means clustering quality")
    legend("bottomright",legend=c("High-entropy cluster","Low-entropy cluster"),pch=15,col=c("gold1","dodgerblue3"))
    par(mfcol=c(1,1))
    dev.off()
  }
  
  return(list(cutoff_surface=seq(length(entropy_clusters$cluster))[entropy_clusters$cluster!=as.integer(entropy_clusters$centers==max(entropy_clusters$centers))[1]+1][1]-1,entropy_classification=as.vector(as.integer(entropy_clusters$cluster==(as.integer(entropy_clusters$centers==min(entropy_clusters$centers))[1]+1))+1)))
}

entropic_analysis <- entropic_ranks(data_under_analysis,population_vector,data_origin,granularity,supervised,process_log,export_plots,create_output_files,is_logged,logbase,huge_feature_list)

htmlwidgets::saveWidget(entropic_analysis$downregulated_classification_topography, "Downregulated entropy topography.html")
htmlwidgets::saveWidget(entropic_analysis$upregulated_classification_topography, "Upregulated entropy topography.html")


#Code written by Hector-Xavier de Lastic
#Development & testing by Hector-Xavier de Lastic & Irene Liampa
#Contact:
#hector.xavier.de.lastic@gmail.com
#irini.liampa@gmail.com
