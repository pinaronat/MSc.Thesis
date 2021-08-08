library(GGally)
library(ggplot2)
library(cluster)
library(Rtsne)
library(dplyr)
################ LOGP, STRUCTURE, P(B) AND FOLDXs#################
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")

all_data_a_beta <- read.csv("abeta_binary_str_plaque.tsv", sep = "\t", header = TRUE)

aggprot_abeta <- all_data_a_beta[,1]
APR_abeta <- all_data_a_beta[,2]
Uniprot_abeta <- all_data_a_beta[,3]
MANGO_abeta <- all_data_a_beta[,4]
main_data_abeta <- all_data_a_beta[,c(60,61,62,78,79,81)] #get the wanted values
main_data_abeta$pb <- main_data_abeta$pb/6 

cols <- c("structure", "interaction")
main_data_abeta[,cols] <- lapply(main_data_abeta[,cols], as.factor) #change selected columns to factors for gower distance calculation

gower_dist_abeta <- as.matrix(daisy(main_data_abeta, metric="gower")) #calculate distance matrix

###CLUSTERING w PAM & tSNE###
#calculate silhouette width with different k using PAM
s_width <- c(NA)
for(i in 2:10){
  pam_fit <- pam(gower_dist_abeta,
                 diss = TRUE,
                 k = i)
  s_width[i] <- pam_fit$silinfo$avg.width
}

#plot the sihouette width
plot(1:10, s_width,
     xlab = "Number of clusters",
     ylab = "Silhouette Width")
lines(1:10, s_width)   #6 clusters is optimal

#clustering
pam_fit <- pam(gower_dist_abeta, diss = TRUE, k = 6) #apply pam with 6 clusters

tsne_obj <- Rtsne(gower_dist_abeta, is_distance = TRUE) #apply tSNE to visualize

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_fit$clustering),
         unip = Uniprot_abeta,mango = MANGO_abeta, APRs = APR_abeta, foldX = main_data_abeta$interaction, logP = main_data_abeta$logP, str = main_data_abeta$structure, pb = main_data_abeta$pb, ci = main_data_abeta$cross_int_energy, elong = main_data_abeta$elongation_energy)

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
  geom_point(aes(color = pb, shape = APRs)) + scale_shape_manual(values = c(16,17,18)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 1.1)


#with structure
ggplot(aes(x = X, y = Y), data = tsne_data) +
  geom_point(aes(color = str, shape = APRs)) + scale_shape_manual(values = c(16,17,18))

#with ci_energy
ggplot(aes(x = X, y = Y, color = ci, shape = APRs), data = tsne_data) +
  geom_point(aes(colour = cut(ci, alpha = 0.3, c(-Inf, 0, 1, Inf)))) +
  scale_color_manual(values = c("(-Inf,0]" = "blue",
                                "(0,1]" = "green",
                                "(1, Inf]" = "red"),
                     labels = c("<= 0", "0 < cross-interaction energy <= 1", "> 1")) + scale_shape_manual(values = c(16,17,18))

#with elongation energy
ggplot(aes(x = X, y = Y, color = elong, shape = APRs), data = tsne_data) +
  geom_point(aes(colour = cut(elong, alpha = 0.3, c(-Inf, 0, 1, Inf)))) +
  scale_color_manual(values = c("(-Inf,0]" = "blue",
                                "(0,1]" = "green",
                                "(1, Inf]" = "red"),
                     labels = c("<= 0", "0 < elongation energy <= 1", "> 1")) + scale_shape_manual(values = c(16,17,18))


all_data_a_beta$X0 <- as.factor(all_data_a_beta$X0)
plaque_data <- all_data_a_beta$X0
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


#save to files
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/GeneOntologyData/Abeta_without_expressions")
write.table(Uniprot1, file = "cluster1.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot2, file = "cluster2.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot3, file = "cluster3.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot4, file = "cluster4.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot5, file = "cluster5.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot6, file = "cluster6.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)
write.table(Uniprot_abeta, file = "all.txt", quote = FALSE, sep = "\n", row.names = FALSE, col.names = FALSE)


