library(umap)
library(GGally)
library(ggplot2)
library(cluster)
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")

all_data_a_beta <- read.csv("a_beta_binary.tsv", sep = "\t", header = TRUE)
all_data_a_beta <- na.omit(all_data_a_beta)

aggprot_abeta <- all_data_a_beta[,1]
APR_abeta <- all_data_a_beta[,2]
Uniprot_abeta <- all_data_a_beta[,3]
MANGO_abeta <- all_data_a_beta[,4]
BBF_data_abeta <- all_data_a_beta[,c(7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78)]

cols <- c("brain_exp", "interaction")
BBF_data_abeta[,cols] <- lapply(BBF_data_abeta[,cols], as.factor)

BBF_dist_abeta <- as.matrix(daisy(BBF_data_abeta, metric="gower"))

custom_umap_abeta <- umap.defaults
custom_umap_abeta$n_neighbors <- 30
custom_umap_abeta$min_dist <- 0.1

abeta_umap <- umap(BBF_dist_abeta, config = custom_umap_abeta)

#umap plot for a-beta
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta)
ggplot(df_abeta, aes(x, y, colour = APRs)) + geom_point()
#with expression data
#olfactory bulb
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta, olfactory_exp = BBF_data_abeta$Tissue.RNA...olfactory.region..NX.)
ggplot(df_abeta, aes(x, y, colour = olfactory_exp, shape = APRs)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = olfactory_exp, shape = APRs)) +
  geom_point(aes(colour = cut(olfactory_exp, c(-Inf, 20, 60, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,60]" = "green",
                                "(60, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 60", "> 60"))
#hippocampus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta, hippocampus_exp = BBF_data_abeta$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_abeta, aes(x, y, colour = hippocampus_exp, shape = APRs)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 60)
ggplot(df_abeta, aes(x, y, colour = hippocampus_exp, shape = APRs)) +
  geom_point(aes(colour = cut(hippocampus_exp, alpha = 0.3, c(-Inf, 20, 90, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,90]" = "green",
                                "(90, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 90", "> 90"))
#cerebral cortex
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta, cerebralcortex_exp = BBF_data_abeta$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_abeta, aes(x, y, colour = cerebralcortex_exp, shape = APRs)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 270)
ggplot(df_abeta, aes(x, y, colour = cerebralcortex_exp, shape = APRs)) +
  geom_point(aes(colour = cut(cerebralcortex_exp, alpha = 0.3, c(-Inf, 30, 300, Inf)))) +
  scale_color_manual(values = c("(-Inf,30]" = "blue",
                                "(30,300]" = "green",
                                "(300, Inf]" = "red"),
                     labels = c("<= 30", "30 < cc_exp <= 300", "> 300"))
#thalamus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta, thalamus_exp = BBF_data_abeta$Tissue.RNA...thalamus..NX.)
ggplot(df_abeta, aes(x, y, colour = thalamus_exp, shape = APRs)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 50)
ggplot(df_abeta, aes(x, y, colour = thalamus_exp, shape = APRs)) +
  geom_point(aes(colour = cut(thalamus_exp, alpha = 0.3, c(-Inf, 25, 80, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,80]" = "green",
                                "(80, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 80", "> 80"))
#foldX interaction data
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta, foldX_interaction = BBF_data_abeta$interaction)
ggplot(df_abeta, aes(x, y, colour = foldX_interaction, shape = APRs)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with logp
logP <- all_data_a_beta[,c()]
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = MANGO_abeta, logp = main_data_abeta$logP)
ggplot(df_abeta, aes(x, y, colour = df_abeta$logp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with brain_dist
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APRs = APR_abeta, braindist = BBF_data_abeta$brain_exp)
ggplot(df_abeta, aes(x, y, colour = braindist, shape = APRs)) + geom_point() + scale_shape_manual(values = c(1,2,5)) 

##GET ABETA CLUSTERS
#big no interaction cluster
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] < 5.3 && abeta_umap$layout[i,2] < 1) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta_nointcluster.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta_nointcluster.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta_nointluster.txt", quote = FALSE, sep = "\n", row.names = FALSE)



#the long line no interaction cluster
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] < -4.9 && abeta_umap$layout[i,2] > 2.5) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta_nointlines.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta_nointlines.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta_nointlines.txt", quote = FALSE, sep = "\n", row.names = FALSE)



#the smallnointcluster
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] < 1 && abeta_umap$layout[i,2] > 10) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta_nointsmallcluster.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta_nointsmallcluster.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta_nointsmallcluster.txt", quote = FALSE, sep = "\n", row.names = FALSE)



#self_cap cluster
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] < 10 && abeta_umap$layout[i,1] > 5 && abeta_umap$layout[i,2] < -2) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else if (abeta_umap$layout[i,2] < 0 && abeta_umap$layout[i,1] > 13) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta_capself.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta_capself.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta_capself.txt", quote = FALSE, sep = "\n", row.names = FALSE)



#tinyclusters
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] > 1.5 && abeta_umap$layout[i,2] > 2.5) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta_tinyclusters.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta_tinyclusters.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta_tinyclusters.txt", quote = FALSE, sep = "\n", row.names = FALSE)
