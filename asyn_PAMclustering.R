#install.packages("Rtsne")
library(GGally)
library(ggplot2)
library(cluster)
library(Rtsne)
library(dplyr)
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")
######################################ALPHA-SYNUCLEIN##################################

############################LOGP, FOLDXs, STRUCTURE, P(B)############################
all_data_asyn <- read.csv("asyn_binary_str_plaque.tsv", sep = "\t", header = TRUE)

aggprot_asyn <- all_data_asyn[,1]
APR_asyn <- all_data_asyn[,2]
Uniprot_asyn <- all_data_asyn[,3]
MANGO_asyn <- all_data_asyn[,4]
main_data_asyn <- all_data_asyn[,c(60,61,62,63,64,66)] #get the wanted values
sum(is.na(main_data_asyn))

#set the binary columns as factors for Gower's Distance Calculation
cols <- c("interaction", "structure")
main_data_asyn[,cols] <- lapply(main_data_asyn[,cols], as.factor) #change selected columns to factors for gower distance calculation

#calculate gower_distance
gower_dist_asyn <- as.matrix(daisy(main_data_asyn, metric="gower")) #calculate distance matrix


###CLUSTERING w PAM & tSNE###
#calculate silhouette width with different k using PAM
s_width <- c(NA)
for(i in 2:10){
  pam_fit <- pam(gower_dist_asyn,
                 diss = TRUE,
                 k = i)
  s_width[i] <- pam_fit$silinfo$avg.width
}

#plot the sihouette width
plot(1:10, s_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, s_width)   ### 8 clusters is optimal

#clustering
pam_fit <- pam(gower_dist_asyn, diss = TRUE, k = 8) #apply PAM with 8 clusters

tsne_obj <- Rtsne(gower_dist_asyn, is_distance = TRUE) #apply tSNE for visualization

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         unip = Uniprot_asyn,mango = MANGO_asyn, APRs = APR_asyn, foldX = main_data_asyn$interaction, logP = main_data_asyn$logP, str = main_data_asyn$structure, pb = main_data_asyn$pb, ci = main_data_asyn$cross_int_energy, elong = main_data_asyn$elongation_energy)

#####plots
#with clusters
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster, shape = APRs)) + scale_shape_manual(values = c(16,17,18))

#with foldx interaction data
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = foldX, shape = APRs)) + scale_shape_manual(values = c(16,17,18))

#with logp
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = logP, shape = APRs)) + scale_shape_manual(values = c(16,17,18)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)

#with pb
#P(b)
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = pb, shape = APRs)) + scale_shape_manual(values = c(16,17,18)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 1.2)

#with structure
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = str, shape = APRs)) + scale_shape_manual(values = c(16,17,18))

#with ci_energy
ggplot(aes(x = X, y = Y, shape = APRs, color = ci), data = tsne_data) +
  geom_point(aes(colour = cut(ci, c(-Inf, 0, 1, Inf)))) +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("(-Inf,0]" = "blue",
                                "(0,1]" = "green",
                                "(1, Inf]" = "red"),
                     labels = c("<= 0", "0 < cross-interaction energy <= 1", "> 1")) 

#with elongation energy
ggplot(aes(x = X, y = Y, shape = APRs, color = elong), data = tsne_data) +
  geom_point(aes(colour = cut(elong, c(-Inf, 0, 1, Inf)))) +
  scale_shape_manual(values = c(16,17,18)) +
  scale_color_manual(values = c("(-Inf,0]" = "blue",
                                "(0,1]" = "green",
                                "(1, Inf]" = "red"),
                     labels = c("<= 0", "0 < elongation energy <= 1", "> 1")) 
#with plaque contents
all_data_asyn$X0 <- as.factor(all_data_asyn$X0)
plaque_data <- all_data_asyn$X0
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = plaque_data, shape = APRs)) + scale_shape_manual(values = c(16,17,18)) + scale_color_manual(values = c("plum2", "royalblue1"))


#get clusters
#cluster1
Uniprot1 <- tsne_data[which(tsne_data$cluster == 1),]$unip

#cluster2
Uniprot2 <- tsne_data[which(tsne_data$cluster == 2),]$unip

#cluster3
Uniprot3 <- tsne_data[which(tsne_data$cluster == 3),]$unip

#cluster5
Uniprot4 <- tsne_data[which(tsne_data$cluster == 4),]$unip

#cluster5
Uniprot5 <- tsne_data[which(tsne_data$cluster == 5),]$unip

#cluster6
Uniprot6 <- tsne_data[which(tsne_data$cluster == 6),]$unip

#cluster7
Uniprot7 <- tsne_data[which(tsne_data$cluster == 7),]$unip

#cluster8
Uniprot8 <- tsne_data[which(tsne_data$cluster == 8),]$unip

#save to files
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/GeneOntologyData/Asyn_without_expressions")
write.table(Uniprot1, file = "cluster1.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot2, file = "cluster2.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot3, file = "cluster3.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot4, file = "cluster4.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot5, file = "cluster5.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot6, file = "cluster6.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot7, file = "cluster7.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot8, file = "cluster8.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot_asyn, file = "all.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)



#####Get the Expressions to Color######
expressions <- all_data_asyn[,c(3,9,12,25,52)]
hippo_exp <- expressions$Tissue.RNA...hippocampal.formation..NX.
ggplot(aes(x = X, y = Y, color = hippo_exp, shape = APRs), data = tsne_data) + scale_shape_manual(values = c(16,17,18)) + geom_point(aes(colour = cut(hippo_exp, c(-Inf, 1, 20, 60, Inf)))) + scale_color_manual(values = c("(-Inf,1]" = "blue3",
                                "(1,20]" = "darkturquoise",
                                "(20,60]" = "darksalmon",
                                "(60, Inf]" = "red"),
                     labels = c("<= 1", "1 < hippocampus exp <= 20", "20 < hippocampus exp <= 60", "> 60"))

#color with plaque contents
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = all_data_asyn$X0, shape = APRs)) + scale_shape_manual(values = c(16,17,18))
