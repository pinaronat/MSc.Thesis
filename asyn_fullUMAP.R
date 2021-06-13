library(umap)
library(GGally)
library(ggplot2)
library(cluster)

setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")
######################################ALPHA-SYNUCLEIN##################################

############################fullUMAP############################
all_data_asyn <- read.csv("a_syn_binary.tsv", sep = "\t", header = TRUE)
all_data_asyn <- na.omit(all_data_asyn)

aggprot_asyn <- all_data_asyn[,1]
APR_asyn <- all_data_asyn[,2]
Uniprot_asyn <- all_data_asyn[,3]
MANGO_asyn <- all_data_asyn[,4]
main_data_asyn <- all_data_asyn[,-c(1,2,3,4)]

#set the binary columns as factors for Gower's Distance Calculation
cols <- c("alzheimer", "parkinson", "als", "brain_exp", "interaction")
main_data_asyn[,cols] <- lapply(main_data_asyn[,cols], as.factor)

#calculate gower_distance
gower_dist_asyn <- as.matrix(daisy(main_data_asyn, metric="gower"))

#umap
custom_umap_asyn <- umap.defaults
custom_umap_asyn$n_neighbors <- 30
custom_umap_asyn$min_dist <- 0.1

asyn_umap <- umap(gower_dist_asyn, config = custom_umap_asyn)

#umap plot for alpha synuclein
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()
#with expression data
#cerebellum
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...cerebellum..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 70)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 90, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,90]" = "green",
                                "(90, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 90", "> 90")) + scale_alpha(range = c(0.3))
#cerebral cortex
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 170)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 250, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,250]" = "green",
                                "(250, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 250", "> 250"))
#hippocampus
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 90)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 130, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,130]" = "green",
                                "(130, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 130", "> 130"))
#olfactory bulb
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...olfactory.region..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 90)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 120, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,120]" = "green",
                                "(120, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 120", "> 120"))
#basal ganglia
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...basal.ganglia..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 150)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 150, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,150]" = "green",
                                "(150, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 150", "> 150"))
#with foldx interaction data
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, interact = main_data_asyn$interaction)
ggplot(df_asyn, aes(x, y, colour = df_asyn$interact, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with logp
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, logp = main_data_asyn$logP)
ggplot(df_asyn, aes(x, y, colour = df_asyn$logp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#brain exp
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, brain_exp = main_data_asyn$brain_exp)
ggplot(df_asyn, aes(x, y, colour = df_asyn$brain_exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))


#####GET THE INTERESTING CLUSTERS

#cluster1
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,2] > -5 && asyn_umap$layout[i,2] < 2 && asyn_umap$layout[i,1] < -2) {
    MANGO_list <- append(MANGO_list, asyn_umap_matrix$MANGO_asyn[i])
    APR_list <- append(APR_list, asyn_umap_matrix$APR_asyn[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, asyn_umap_matrix$Uniprot_asyn[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = color_list)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_asyn_cluster1.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_asyn_cluster1.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_asyn_cluster1.txt", quote = FALSE, sep = "\n", row.names = FALSE)

#cluster2
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,2] > 0 && asyn_umap$layout[i,1] > -3 && asyn_umap$layout[i,1] < 5) {
    MANGO_list <- append(MANGO_list, asyn_umap_matrix$MANGO_asyn[i])
    APR_list <- append(APR_list, asyn_umap_matrix$APR_asyn[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, asyn_umap_matrix$Uniprot_asyn[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = color_list)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_asyn_cluster2.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_asyn_cluster2.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_asyn_cluster2.txt", quote = FALSE, sep = "\n", row.names = FALSE)

#cluster3
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,2] < 0 && asyn_umap$layout[i,1] > 5 && asyn_umap$layout[i,1] > -5) {
    MANGO_list <- append(MANGO_list, asyn_umap_matrix$MANGO_asyn[i])
    APR_list <- append(APR_list, asyn_umap_matrix$APR_asyn[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, asyn_umap_matrix$Uniprot_asyn[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = color_list)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_asyn_cluster3.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_asyn_cluster3.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_asyn_cluster3.txt", quote = FALSE, sep = "\n", row.names = FALSE)


#cluster4
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,1] > 0 && asyn_umap$layout[i,1] < 5 && asyn_umap$layout[i,2] < -5) {
    MANGO_list <- append(MANGO_list, asyn_umap_matrix$MANGO_asyn[i])
    APR_list <- append(APR_list, asyn_umap_matrix$APR_asyn[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, asyn_umap_matrix$Uniprot_asyn[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = color_list)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_asyn_cluster4.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_asyn_cluster4.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_asyn_cluster4.txt", quote = FALSE, sep = "\n", row.names = FALSE)






