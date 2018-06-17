library(ape)
library(ggplot2)
library(ggtree)
library(treeman)
library(dplyr)
library(plyr)
library(reshape)
library(devtools)
#install_github("ggbiplot", "vqv")
library(ggbiplot)
library(Rmisc)

calc_entropy<-function(meta_5, per_plant_merged){
  #the goal of this function is to takee every plant and calculate the tree entropy diversity
  no_id=within(per_plant_merged, rm("plant_id"))
  entropy_append=numeric(length(meta_5))
  names(entropy_append)=meta_5
  for(i in 1:length(meta_5)){
    hm=no_id[which(no_id$meta_genome_id==meta_5[i]),2:dim(no_id)[2]]
    if(dim(hm)[1]>0){
      val=hm[which(hm!=0)]/sum(hm[which(hm!=0)])
      print(val)
      entr=-sum(val*log(val))
      entropy_append[i]=entr
    }
    else{
      entropy_append[i]=NA
      print(meta_5[i])
    }
  }
  return(entropy_append)
}

apply_clade<-function(rec){
  val=rec["sample_tot"]
  if(is.na(match(val,syringae_nodes))==F){
    rec["clade"]="mediumpurple"#brewer.pal(8,"Purples")[4]
  }
  else{
    rec["clade"]="BLACK"}
  return(rec)
}


plot_trees<-function(ent_sort, season){
  postscript(paste(paste("~/Dropbox/germany_pathogen_collections/data_files_rmarkdown/trees", season, sep="_"),".eps", sep=""), height=10, width=200)
  #par(mfrow=c(9,10))
  #now assign phylogeny per plan()
  layout(matrix(1:30, 5, 6))
  par(mar=c(0.3,0.3,0.3,0.3))
  
  #div_ratio=numeric(length(ent_sort))
  names(div_ratio)=names(ent_sort)
  for(val in names(ent_sort)[1:20]){
    print(val)
    keep=meta[c(which(meta$meta_genome_id==val)),]$plate_seq
    color_samp=(meta[c(which(meta$meta_genome_id==val)),])
    hm=assign(paste0("Plant",val), keep)
    match=-match(hm, tree$tip.label)
    pruned.tree<-drop.tip(tree,tree$tip.label[match[is.na(match)==F]])
    matched_color=match( pruned.tree$tip.label,color_samp$plate_seq)
    matched_color=matched_color[matched_color!="NA"]
    reorder_color=as.numeric(color_samp[matched_color,]$clade)
    #black=length(reorder_color[reorder_color=="F"])
    #red=length(reorder_color[reorder_color=="T"])
    #black_ratio=black/(red+black)
    #print(black_ratio)
    #div_ratio[val]=(black_ratio)
    if(length(match[is.na(match)==F])>4){
      plot.phylo(ladderize(pruned.tree), show.tip.label=F, x.lim=1)#, main=val)
      tiplabels(pch=21, col="white", bg=reorder_color, cex=1.2)
      if(val=="S3" || val=="S185"){
        add.scale.bar(.1, 5)
      }#add.scale.bar(.1, 5)
      #
    }}
  dev.off()
  #return(div_ratio)
}