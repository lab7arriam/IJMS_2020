library(ggtree)
library(ggplot2)
library(stringr)
library(ape)
library(dendextend)
library(lattice)
library(cluster)
library(vegan)
library(dplyr)
library(reshape2)
library("ggsci")
library(ggfortify)
library("DECIPHER")

#heatmap function
drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    # panel.text(x, y, e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet,
                   at=seq(0, max(df), length.out=100),
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1)),
                   scales=list(cex=0.7,x=list(rot=45),tck = c(0,0)), xlab=list(label=''),
                   ylab=list(label=''), panel=myPanel_a))
}

setwd("/home/anton/bt_pangenome_project/tree_stuffs/pics/Pics_V4")
setwd("/home/anton/bt_pangenome_project/tree_stuffs/all_trees")

strains_table <- read.csv('/home/anton/bt_pangenome_project/filtered_annotated_strains.csv', header = T, sep = "\t") 
core_snp_tree <-read.tree("ML_core_iter_4.raxml.support")
binary_tree <-read.tree("accessory_binary_genes.fa.newick")


#modify roary-genegated tip labels
core_snp_tree$tip.label <- paste(str_split(core_snp_tree$tip.label, "_", simplify = TRUE)[,1],
       str_split(core_snp_tree$tip.label, "_", simplify = TRUE)[,2],
       sep = "_")

#change labels to strain names
core_snp_tree$tip.label <- make.names(as.character(strains_table[match(core_snp_tree$tip.label, strains_table$GenBank_accession),3]))
ggtree(core_snp_tree)+ geom_tiplab(size=2.5, show.legend=T)+ 
  geom_nodelab(size=2.6,hjust=2)

#apply groping for coloured trees
groupInfo <- split(core_snp_tree$tip.label, as.character(strains_table[match(core_snp_tree$tip.label, strains_table$GenBank_accession),3])) <- as.character(strains_table[match(core_snp_tree$tip.label, strains_table$GenBank_accession),3])
tr_gr <- groupOTU(core_snp_tree, groupInfo)

groupInfo$aizawai

#ggtree(tr_gr, aes(color=group)) +
#  theme(legend.position="right") + geom_tiplab(size=1, , show.legend=T)


binary_tree$tip.label <- paste(str_split(binary_tree$tip.label, "_", simplify = TRUE)[,1],
                                 str_split(binary_tree$tip.label, "_", simplify = TRUE)[,2],
                                 sep = "_")

binary_tree$tip.label <- as.character(strains_table[match(binary_tree$tip.label, strains_table$GenBank_accession),3])
ggtree(binary_tree)+ geom_tiplab(size=2.5, show.legend=T)

#Iteratively drawing trees
#13_9

file_list=list.files(pattern = ".support$", recursive = TRUE)
for (file in file_list){
  tree <- read.tree(file)
  out_name <- paste0(file, ".pdf")
  print(out_name)
  tree$tip.label <- paste(str_split(tree$tip.label, "_", simplify = TRUE)[,1],
                                   str_split(tree$tip.label, "_", simplify = TRUE)[,2],
                                   sep = "_")
  tree$tip.label <- make.names(as.character(strains_table[match(tree$tip.label, 
                                                                strains_table$GenBank_accession),3]))
  p <- ggtree(tree)+ geom_tiplab(size=2.5, show.legend=T)+ 
    geom_nodelab(size=2.6,hjust=2)
  ggsave(filename=out_name,width=13, height=9, dpi=300)
}




#flagellin trees

flag_tree <-read.tree("/home/anton/bt_pangenome_project/tree_stuffs/all_trees/ML_flag_nucl_untr.raxml.support")
flag_tree $tip.label <- paste(str_split(flag_tree $tip.label, "_", simplify = TRUE)[,1],
                               str_split(flag_tree $tip.label, "_", simplify = TRUE)[,2],
                               sep = "_")

flag_tree $tip.label <- as.character(strains_table[match(flag_tree$tip.label, strains_table$GenBank_accession),3])
ggtree(core_snp_tree)+ geom_tiplab(size=2.5, show.legend=T)

core_snp_tree <-ape::read.tree("core_snp_renamed.nwk")
flag_tree <-ape::read.tree("flagellin_renamed.nwk")
flag_tree $tip.label <- paste(str_split(flag_tree $tip.label, "_", simplify = TRUE)[,1],
                              str_split(flag_tree $tip.label, "_", simplify = TRUE)[,2],
                              sep = "_")
core_snp_tree$tip.label <- paste(str_split(core_snp_tree$tip.label, "_", simplify = TRUE)[,1],
                                 str_split(core_snp_tree$tip.label, "_", simplify = TRUE)[,2],
                                 sep = "_")
flag_tree $node.label <- NULL
core_snp_tree $node.label <- NULL
write.tree(flag_tree ,file="flagellin_no_supp.nwk")
write.tree(core_snp_tree ,file="core_no_supp.nwk")


#library('DECIPHER')
dend2 <- ReadDendrogram(file="/home/anton/bt_pangenome_project/tree_stuffs/all_trees/flagellin_no_supp_renamed.nwk") 
dend1 <- ReadDendrogram(file="/home/anton/bt_pangenome_project/tree_stuffs/all_trees/core_no_supp_renamed.nwk") 

## rearrange in ladderized fashion
dnd1 <- ladderize(sort(dend1) )
dnd2 <- ladderize(sort(dend2))
## plot the tanglegram
dndlist <- dendextend::dendlist(dnd1, dnd2)


dndlist %>%set("rank_branches") %>% 
  untangle(method = "DendSer") %>%
  untangle(method = "step1side", k_seq = 2:5) %>% 
  untangle(method = "ladderize",max_n_iterations=60) %>% 
  tanglegram(common_subtrees_color_branches = TRUE,margin_inner = 12.1,
             highlight_distinct_edges = F,
             highlight_branches_lwd = F,
             axes=F,sort=T,lab.cex = 1.75,edge.lwd=3,margin_outer=8,margin_bottom=0.00001)

#dendextend::tanglegram(dndlist, fast = TRUE, margin_inner = 25,common_subtrees_color_branches = TRUE)



gyr_tree <-read.tree("/home/anton/bt_pangenome_project/tree_stuffs/gyrase/ML_gyrAB_nucl.raxml.support")
gyr_tree$tip.label <- paste(str_split(gyr_tree$tip.label, "_", simplify = TRUE)[,1],
                              str_split(gyr_tree$tip.label, "_", simplify = TRUE)[,2],
                              sep = "_")

gyr_tree$tip.label <- as.character(strains_table[match(gyr_tree$tip.label, strains_table$GenBank_accession),3])
ggtree(gyr_tree)+ geom_tiplab(size=2.5, show.legend=T)


nprb_tree <- read.tree("/home/anton/bt_pangenome_project/tree_stuffs/nprB/nucl/ML_nprb.raxml.support.renamed.nwk")
nprb_tree$tip.label <- paste(str_split(nprb_tree$tip.label, "_", simplify = TRUE)[,1],
                            str_split(nprb_tree$tip.label, "_", simplify = TRUE)[,2],
                            sep = "_")

nprb_tree$tip.label <- as.character(strains_table[match(nprb_tree$tip.label, strains_table$GenBank_accession),3])
ggtree(nprb_tree) + geom_tiplab(size=2.5, show.legend=T)



#Similarity heatmaps
#/home/anton/bt_pangenome/assemblies_distance/mash_dist_filt.tsv
#/home/anton/bt_pangenome_project/tree_stuffs/flagellin_id_heatmap.csv
#/home/anton/bt_pangenome_project/full_genome_alignments/full_genomes_corrected.tsv
genome_array_matrix <- read.table('/home/anton/bt_pangenome_project/assemblies_distance/mash_dist_filt.tsv', row.names=1,header = TRUE,  sep = "\t")

genome_array_matrix <- 1-genome_array_matrix


min(genome_array_matrix )
genome_array_matrix  %>% colMeans() %>% 


mena(typeof(genome_array_matrix ))
max(genome_array_matrix )
#perform hieracheal clustering
clusters <- hclust(as.dist(as.matrix(genome_array_matrix)), method = "complete")
#?hclust
plot(clusters, labels=FALSE)


#find optimal number of clusters
s=NULL
for(i in 10:100){
  k=cutree(clusters,i)
  s[i]=summary(silhouette(k,as.dist(genome_array_matrix)))$si.summary[3]
}
which.max(s)
clusterCut <- cutree(clusters, 45)

#reorder genome heatmap with accordance to ordered hc dendrogram
dend <- as.phylo(clusters)
#write.tree(dend,file="hclust_complete_full_genome_corrected.nwk")


names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
dend <- as.dendrogram(clusters)
dend  = reorder(dend, wts = order(match(names_frame$keyName, rownames(genome_array_matrix))))
ord_dend <- rev(reorder(dend, agglo.FUN=sum, 104:1))
plot(dend)
plot(ord_dend)
check <- genome_array_matrix[match(as.vector(labels(ord_dend)),rownames(genome_array_matrix)),match(as.vector(labels(ord_dend)),rownames(genome_array_matrix))]

#plot heatmap
mypal = colorRampPalette(c('#FF4343','#FAFAFA'))
jet = mypal(100)

rownames(check) <- paste(str_split(rownames(check), "_", simplify = TRUE)[,1],
                         str_split(rownames(check), "_", simplify = TRUE)[,2],
                         sep = "_")

colnames(check) <- rownames(check) 

rownames(check) <-make.names(as.character(strains_table[match(rownames(check), strains_table$GenBank_accession),3]), unique = TRUE)
colnames(check) <- rownames(check) 

drawCoolHM(as.matrix(check))

dend <- as.phylo(clusters)
dend$tip.label

dend$tip.label <- paste(str_split(dend$tip.label, "_", simplify = TRUE)[,1],
                         str_split(dend$tip.label, "_", simplify = TRUE)[,2],
                         sep = "_")

dend$tip.label <-make.names(as.character(strains_table[match(dend$tip.label , strains_table$GenBank_accession),3]), unique = TRUE)
ggtree(dend)+ geom_tiplab(size=2.5, show.legend=T)


#PCA phoresis
#/home/anton/bt_pangenome_project/tree_stuffs/Spot_PCA.csv
#/home/anton/bt_pangenome_project/phoresis/dedupped_3/unique_DIGE_names.csv

prot_data = read.table('/home/anton/bt_pangenome_project/phoresis/dedupped_3/unique_DIGE_names.csv',
                       header=T, sep=',', stringsAsFactors = F,  row.names = 1)
prot_data <- prot_data[,-c(6)]
sc_prot<-data.frame(t(prot_data))


wss <- (nrow(sc_prot)-1)*sum(apply(sc_prot,2,var))
for (i in 1:4) wss[i] <- sum(kmeans(sc_prot,
                                       centers=i)$withinss,nstart=25,iter.max=1000)
plot(1:4, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

wss_gf=data.frame(clust=c(1:4), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                    size=16)
  ) 

prot_kmeans <- kmeans(sc_prot,centers=2)
prot_kmeans$cluster
sc_prot$name=rownames(sc_prot)

dist <- vegdist(t(prot_data),  method = "bray")

pcaout = prcomp(dist)
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
pc12$fac=rownames(pc12)
pc12$clust = as.factor(prot_kmeans$cluster)


sc_prot$clust= as.factor(prot_kmeans$cluster)

#autoplot
autoplot(prot_kmeans, 
         data=sc_prot, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 9,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B"))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                    size=24),
         axis.title.y=element_text(color="black", 
                                   size=28),
         panel.background = element_blank(), 
         axis.line = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(color='black', 
                                    size=24),
         axis.title.x = element_text(color="black", 
                                     size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 
scale_y_continuous(expand = c(.1, .1))



ggplot(pc12, aes(x=-PC1, y=PC2,col=fac,shape=clust)) + 
  geom_point(size=5) +
  theme_bw() + xlab('PC1 -59% of variance') +
  guides(fill=guide_legend(title="Sample"))+ 
  ylab('PC2 - 18% of variance')+
  theme( axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(color='black', 
                                    size=12),
         axis.title.x = element_text(face="bold", color="black", 
                                     size=14)
  )

ggplot(df, aes(own, method)) + 
  geom_point(colour="white", shape=21, size = 4, 
             aes(fill = factor(label))) + 
  scale_fill_manual(values=c("blue", "cyan4"))

pc12 = as.data.frame(PCOA$vectors)
pc12$fac=rownames(pc12)
ggplot(pc12, aes(x=-Axis.1, y=Axis.2, fill=fac)) + geom_point(shape=21,size=4) +
  theme_bw() 



#serovas PCA
serovar_data = read.table('/home/anton/bt_pangenome_project/tree_stuffs/serovars_cons_all_trees.tsv', row.names=1,
                       header=T, sep='\t', stringsAsFactors = F)
dist <- vegdist(t(serovar_data),  method = "euclidean")

wss <- (nrow(serovar_data)-1)*sum(apply(serovar_data,2,var))
for (i in 2:13) wss[i] <- sum(kmeans(serovar_data,
                                    centers=i)$withinss,nstart=25,iter.max=1000)

plot(1:13, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
wss_gf=data.frame(clust=c(1:13), wss=wss)

ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =3, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 


pcaout = prcomp(dist)
PCOA = pcoa(dist)
imps = summary(pcaout)$importance

pc12 = as.data.frame(pcaout$x[, 1:2])
pc12$fac=rownames(pc12)
ggplot(pc12, aes(x=-PC1, y=PC2, fill=fac)) + geom_point(shape=21,size=4) +
  theme_bw() 

pc12 = as.data.frame(PCOA$vectors)
pc12$fac=rownames(pc12)
ggplot(pc12, aes(x=-Axis.1, y=Axis.2, fill=fac)) + geom_point(shape=21,size=4) +
  theme_bw() 

sc_ser<-data.frame(t(serovar_data ))
ser_kmeans <- kmeans(sc_ser,centers=3)

#autoplot
autoplot(ser_kmeans, 
         data=sc_ser, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 5,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black'))+
  scale_fill_manual(values = c('#A60B0B','#A9A9A9','#2980B9'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28)
  )+ 
  scale_x_continuous(expand = c(.1, .1)) +
scale_y_continuous(expand = c(.1, .1))


serovar_data <- rbind(serovar_data,colSums(serovar_data))
rownames(serovar_data) <- c(rownames(serovar_data)[-15],'ser_sum')

melt_serovar_data <-melt(as.matrix(serovar_data)) %>%  filter(Var1=='ser_sum')
median(melt_serovar_data$value)

ggplot(data = melt_serovar_data , aes(Var2, value))+
  geom_bar(stat = 'identity', aes(x = reorder(Var2, -value)),
           fill='#A9A9A9',alpha=0.75,color='black')+
  theme_bw() + xlab('Trees') + ylab('Sum of leaves') +
  geom_hline(yintercept =670,linetype="dashed", color = "#2980B9",size =2.5)+
  theme( axis.text.y = element_text(color='black', 
                                    size=24),
         axis.title.y=element_text(color="black", 
                                   size=28),
         panel.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 24, hjust = 1),
         axis.title.x = element_text(color="black", 
                                     size=28)
  )


#Tree topology matrix
#tree_dist_flagellins.tsv.csv tree_dist_matrix.tsv.csv
tree_dist = read.table('/home/anton/bt_pangenome_project/tree_stuffs/tree_dist_matrix.tsv.csv', row.names=1,
                       header=T, sep=',', stringsAsFactors = F)
tree_dist <- round(tree_dist,2)
ord_rows=c("bin","core","mash","minimap","gyr","mmsA" ,"guaB","sucC" ,"dnaK",
           "nprB" , "rph", "atpD","yjlD" , "fusA","groEL","rocA" ,"phbB","inhA","tuf",
           "calY","flagellin")
tree_dist <- tree_dist[ord_rows,ord_rows]

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

tree_dist <- get_upper_tri(tree_dist)


melted_dist <- melt(as.matrix(tree_dist), na.rm = TRUE)


mean_id_for_genes <- melted_dist %>% filter(Var1 %in% c('bin','core','mash','minimap')) %>% group_by(Var2) %>% summarise(m=mean(value))
rownames(tree_dist)



ggplot(data = as.data.frame(mean_id_for_genes)[-c(1:4),], aes(Var2, m))+
  geom_bar(stat = 'identity', aes(x = reorder(Var2, -m)),
           color='black',fill='#A60B0B',alpha=0.75)+
  theme_bw() + ylab('Mean topological similarity') + xlab('Trees') +
  geom_hline(yintercept =0.79,linetype="dashed", color = "#2980B9",size =2.5)+
  theme( axis.text.y = element_text(color='black', 
                                    size=24),
         axis.title.y=element_text( color="black", 
                                   size=28),
         panel.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 24, hjust = 1),
         axis.title.x = element_text(color="black", 
                                     size=28)
  )
median(mean_id_for_genes [-c(1:10),3])

ggplot(data = melted_dist, aes(Var2, Var1, fill = value))+
  geom_tile(color = 'black')+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4.5)+
  scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(0.3,1), space = "Lab", 
                       name="Tree simillarity") +
  theme_bw() + xlab('Trees') + ylab('Trees') +
  theme( axis.text.y = element_text(color='black', 
                                    size=22),
         axis.title.y=element_text(color="black", 
                                   size=24),
         panel.background = element_rect(fill = "#DCDCDC"),
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                         colour = "black"),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 22, hjust = 1),
         axis.title.x = element_text( color="black", 
                                     size=24))+
  coord_fixed()



# Gene presense/absence

binary_dige_dat = read.table('/home/anton/bt_pangenome_project/phoresis/dedupped_3/summarized_proteins_67.tsv', row.names=1,
                       header=T, sep='\t', stringsAsFactors = F)

binary_dige_dat[binary_dige_dat!=0 ]<-1
binary_dige_dat_appl <- data.frame(matrix(unlist(lapply(binary_dige_dat,as.numeric)), nrow=104, byrow=F))
  
rownames(binary_dige_dat_appl) <- rownames(binary_dige_dat) 
colnames(binary_dige_dat_appl) <- colnames(binary_dige_dat) 

binary_dige_dat<-binary_dige_dat_appl

strains_table$GenBank_accession <- as.character(str_split(strains_table$GenBank_accession, '[.]', simplify = TRUE)[,1])
rownames(binary_dige_dat) <- make.names(as.character(strains_table[match(rownames(binary_dige_dat), strains_table$GenBank_accession),3]),unique = T)


binary_dige_dat <-binary_dige_dat[ order(row.names(binary_dige_dat)), ]
binary_dige_dat_filt <- binary_dige_dat[,colSums(binary_dige_dat>0)>0]
binary_dige_dat_filt <- binary_dige_dat_filt[,colSums(binary_dige_dat_filt>0)<103]

melt.dige<-melt(as.matrix(binary_dige_dat_filt ))


binary_dat_isr <- binary_dige_dat[rownames(binary_dige_dat) %in% c("israelensis","israelensis.1",
                                                 'israelensis.2'),] %>% colMeans()

binary_dat_dar <- binary_dige_dat[rownames(binary_dige_dat) %in% c('darmstadiensis','darmstadiensis.1'),] %>% colMeans()

binary_dat_thur <- binary_dige_dat[rownames(binary_dige_dat) %in% c('thuringiensis','thuringiensis.1',
                                                                    'thuringiensis.2'),] %>% colMeans()

genomes_tile <- rbind(binary_dat_thur,binary_dat_dar,binary_dat_isr,binary_dat_isr,binary_dat_isr)
rownames(genomes_tile) <- c('T','D',
                            "I","A",'V') 

prot_data_vis <-prot_data[,-c(6)]
prot_data_vis <-as.data.frame(t(prot_data_vis ))
prot_data_vis <-prot_data_vis[ ,colnames(genomes_tile)] 

merged_dat_vis <- prot_data_vis
for (i in (1:nrow(merged_dat_vis))){
  for (j in (1:ncol(merged_dat_vis))){
    if (prot_data_vis [i,j]==0){
      merged_dat_vis[i,j] <- -genomes_tile[i,j]
    }
    else{
      merged_dat_vis[i,j] <- prot_data_vis [i,j]-genomes_tile[i,j]/5
    }
  }
}

merged_dat_vis[merged_dat_vis==0] <-5
merged_dat_vis[merged_dat_vis==-1] <-0
merged_dat_vis[merged_dat_vis==5] <- -1
merged_dat_vis=round(merged_dat_vis,2)
merged_dat_vis[merged_dat_vis==0.93] <- 0.85
merged_dat_vis[merged_dat_vis==0.8] <- 0.71

melted_merged <- melt(as.matrix(merged_dat_vis ))
genomes_melted <- melt(as.matrix(genomes_tile))
prots_melted <- melt(as.matrix(prot_data_vis ))

ggplot(data =melted_merged   , aes(Var2, Var1, fill = value))+
  geom_tile(color = "black")+
  scale_fill_gradient2(low = "#DCDCDC", high = "yellow", mid = "#A60B0B", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Presence")+
  theme_minimal() + xlab('Proteins') + ylab('Serovars') +
  theme( axis.text.y = element_text(color='black', 
                                    size=24),
         axis.title.y=element_text( color="black", 
                                   size=28),
         panel.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 24, hjust = 1),
         axis.title.x = element_text( color="black", 
                                     size=28))+
  coord_fixed()



ggplot(data =melt.dige , aes(Var1, Var2, fill = value))+
  geom_tile(color = "black", alpha=0.9)+
  scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(0,1), space = "Lab", 
                       name="Presence")+
  coord_fixed()+
  theme_minimal() + xlab('Serovars') + ylab('Proteins') +
  theme( axis.text.y = element_text(color='black', 
                                    size=24),
         axis.title.y=element_text(color="black", 
                                   size=28),
         panel.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 24, hjust = 1),
         axis.title.x = element_text(color="black", 
                                     size=28)
  )

  
