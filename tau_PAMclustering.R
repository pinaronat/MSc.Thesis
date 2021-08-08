library(umap)
library(GGally)
library(ggplot2)
library(cluster)
library(Rtsne)
library(dplyr)
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")
all_data_tau <- read.csv("tau_binary_str_plaque.tsv", sep = "\t", header = TRUE)

aggprot_tau <- all_data_tau[,1]
APR_tau <- all_data_tau[,2]
Uniprot_tau <- all_data_tau[,3]
MANGO_tau <- all_data_tau[,4]
main_data_tau <- all_data_tau[,c(60,61,62,78,79,81)] #get the wanted values
main_data_tau$pb <- main_data_tau$pb/6

sum(is.na(main_data_tau))
cols <- c("structure", "interaction")
main_data_tau[,cols] <- lapply(main_data_tau[,cols], as.factor) #change selected columns to factors for gower distance calculation

gower_dist_tau <- as.matrix(daisy(main_data_tau, metric="gower")) #calculate distance matrix

###CLUSTERING w PAM & tSNE###
#calculate silhouette width with different k using PAM
s_width <- c(NA)
for(i in 2:10){
  pam_fit <- pam(gower_dist_tau,
                 diss = TRUE,
                 k = i)
  s_width[i] <- pam_fit$silinfo$avg.width
}

#plot the sihouette width
plot(1:10, s_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, s_width)   #3 clusters is optimal

#clustering
pam_fit <- pam(gower_dist_tau, diss = TRUE, k = 3) #apply PAM with 3 clusters

tsne_obj <- Rtsne(gower_dist_tau, is_distance = TRUE) #apply tSNE for visualization

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         unip = Uniprot_tau,mango = MANGO_tau, APRs = APR_tau, foldX = main_data_tau$interaction, logP = main_data_tau$logP, str = main_data_tau$structure, pb = main_data_tau$pb, ci = main_data_tau$cross_int_energy, elong = main_data_tau$elongation_energy)

#####plots
#with clusters
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = cluster, shape = APRs)) + scale_shape_manual(values = c(16,17))

#with foldx interaction data
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = foldX, shape = APRs)) + scale_shape_manual(values = c(16,17))

#with logp
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = logP, shape = APRs)) + scale_shape_manual(values = c(16,17)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)

#with pb
#P(b)
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = pb, shape = APRs)) + scale_shape_manual(values = c(16,17)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 1.2)

#with expression
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = main_data_tau$Tissue.RNA...amygdala..NX.))

#color with plaque contents
nft_data = as.factor(all_data_tau$X0)
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = nft_data, shape = APRs)) + scale_shape_manual(values = c(16,17)) + scale_color_manual(values = c("plum2", "royalblue1"))

#with structure
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = str, shape = APRs)) + scale_shape_manual(values = c(16,17))

#with ci_energy
ggplot(aes(x = X, y = Y, color = ci, shape = APRs), data = tsne_data) +
  geom_point(aes(colour = cut(ci, alpha = 0.3, c(-Inf, 0, 1, Inf)))) +
  scale_shape_manual(values = c(16,17)) +
  scale_color_manual(values = c("(-Inf,0]" = "blue",
                                "(0,1]" = "green",
                                "(1, Inf]" = "red"),
                     labels = c("<= 0", "0 < cross-interaction energy <= 1", "> 1"))

#with elongation energy
ggplot(aes(x = X, y = Y, color = elong, shape = APRs), data = tsne_data) +
  geom_point(aes(colour = cut(elong, alpha = 0.3, c(-Inf, 0, 1, Inf)))) +
  scale_shape_manual(values = c(16,17)) +
  scale_color_manual(values = c("(-Inf,0]" = "blue",
                                "(0,1]" = "green",
                                "(1, Inf]" = "red"),
                     labels = c("<= 0", "0 < elongation energy <= 1", "> 1"))

#get clusters
#cluster1
Uniprot1 <- tsne_data[which(tsne_data$cluster == 1),]$unip

#cluster2
Uniprot2 <- tsne_data[which(tsne_data$cluster == 2),]$unip

#cluster3
Uniprot3 <- tsne_data[which(tsne_data$cluster == 3),]$unip

#cluster2-left
clust2 <- tsne_data[which(tsne_data$cluster == 2),] #get whole cluster2

color_list <- c()
Uniprot_list <- c()
for (i in 1:length(clust2$Y)) {
  if (clust2$X[i] <= 13 && clust2$X[i] > -19) {
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, clust2$unip[i])
  } 
  else {
    color_list <- append(color_list, "gray")
  }
}

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = color_list))

#cluster2-right +small ones
color_list <- c()
Uniprot_list <- c()
for (i in 1:length(tsne_data$Y)) {
  if (tsne_data$X[i] > 13) {
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, tsne_data$unip[i])
  } else if (tsne_data$Y[i] > 19 && tsne_data$X[i] < -20) {
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, tsne_data$unip[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = color_list))

#save to files
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/GeneOntologyData/tau_without_expressions/ClustersUnip")
write.table(Uniprot1, file = "cluster1.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot2, file = "cluster2.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot3, file = "cluster3.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot_list, file = "cluster2left.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot_tau, file = "all.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)


