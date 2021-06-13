library(umap)
library(GGally)
library(ggplot2)
library(cluster)

setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")
############SEPARATE UMAP FOR EACH AGG PROT###############
all_data_a_syn <- read.csv("a_syn_binary.tsv", sep = "\t", header = TRUE)
all_data_tau <- read.csv("tau_binary.tsv", sep = "\t", header = TRUE)
all_data_a_beta <- read.csv("a_beta_binary.tsv", sep = "\t", header = TRUE)
all_data_sod1 <- read.csv("sod_1_binary.tsv", sep = "\t", header = TRUE)

all_data_a_syn <- na.omit(all_data_a_syn)
all_data_a_beta <- na.omit(all_data_a_beta)
all_data_tau <- na.omit(all_data_tau)
all_data_sod1 <- na.omit(all_data_sod1)

aggprot_sod1 <- all_data_sod1[,1]
APR_sod1 <- all_data_sod1[,2]
Uniprot_sod1 <- all_data_sod1[,3]
MANGO_sod1 <- all_data_sod1[,4]
main_data_sod1 <- all_data_sod1[,-c(1,2,3,4)]

aggprot_abeta <- all_data_a_beta[,1]
APR_abeta <- all_data_a_beta[,2]
Uniprot_abeta <- all_data_a_beta[,3]
MANGO_abeta <- all_data_a_beta[,4]
main_data_abeta <- all_data_a_beta[,-c(1,2,3,4)]

aggprot_asyn <- all_data_a_syn[,1]
APR_asyn <- all_data_a_syn[,2]
Uniprot_asyn <- all_data_a_syn[,3]
MANGO_asyn <- all_data_a_syn[,4]
main_data_asyn <- all_data_a_syn[,-c(1,2,3,4)]

aggprot_tau <- all_data_tau[,1]
APR_tau <- all_data_tau[,2]
Uniprot_tau <- all_data_tau[,3]
MANGO_tau <- all_data_tau[,4]
main_data_tau <- all_data_tau[,-c(1,2,3,4)]


#set the binary columns as factor
cols <- c("alzheimer", "parkinson", "als", "brain_exp", "interaction")
main_data_tau[,cols] <- lapply(main_data_tau[,cols], as.factor)
main_data_abeta[,cols] <- lapply(main_data_abeta[,cols], as.factor)
main_data_asyn[,cols] <- lapply(main_data_asyn[,cols], as.factor)
main_data_sod1[,cols] <- lapply(main_data_sod1[,cols], as.factor)

#calculate gower_distance
gower_dist_sod1 <- as.matrix(daisy(main_data_sod1, metric="gower"))
gower_dist_abeta <- as.matrix(daisy(main_data_abeta, metric="gower"))
gower_dist_asyn <- as.matrix(daisy(main_data_asyn, metric="gower"))
gower_dist_tau <- as.matrix(daisy(main_data_tau, metric="gower"))

#umap
custom_umap_sod1 <- umap.defaults
custom_umap_sod1$n_neighbors <- 15
custom_umap_sod1$min_dist <- 0.1

custom_umap_asyn <- umap.defaults
custom_umap_asyn$n_neighbors <- 30
custom_umap_asyn$min_dist <- 0.1

custom_umap_tau <- umap.defaults
custom_umap_tau$n_neighbors <- 30
custom_umap_tau$min_dist <- 0.1

custom_umap_abeta <- umap.defaults
custom_umap_abeta$n_neighbors <- 45
custom_umap_abeta$min_dist <- 0.1

asyn_umap <- umap(gower_dist_asyn, config = custom_umap_asyn)
tau_umap <- umap(gower_dist_tau, config = custom_umap_tau)
sod1_umap <- umap(gower_dist_sod1, config = custom_umap_sod1)
abeta_umap <- umap(gower_dist_abeta, config = custom_umap_abeta)

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

#umap plot for sod1
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1)
ggplot(df_sod1, aes(x, y, colour = df_sod1$APR)) + geom_point()
#with expression value (spinal cord)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, exp = main_data_sod1$Tissue.RNA...spinal.cord..NX.)
ggplot(df_sod1, aes(x, y, colour = df_sod1$exp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 40, breaks = c(0,15, 50))
#with interaction value(foldx)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, inter = main_data_sod1$interaction)
ggplot(df_sod1, aes(x, y, colour = df_sod1$inter, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))
#logP
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = MANGO_sod1, logp = main_data_sod1$logP)
ggplot(df_sod1, aes(x, y, colour = df_sod1$logp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#delta hydrophobicity
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = MANGO_sod1, hyd = main_data_sod1$Delta_hydrophob)
ggplot(df_sod1, aes(x, y, colour = df_sod1$hyd, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)

#umap plot for tau
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()
#with expression data
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, exp = main_data_tau$Tissue.RNA...amygdala..NX.)
ggplot(df_tau, aes(x, y, colour = df_tau$exp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 70)
#with logp
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = MANGO_tau, logp = main_data_tau$logP)
ggplot(df_tau, aes(x, y, colour = df_tau$logp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#with interaction value (foldX)
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, inter = main_data_tau$interaction)
ggplot(df_tau, aes(x, y, colour = inter, shape = df_tau$APR)) + geom_point()
#with hydrophobicity
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = MANGO_tau, hydrob = main_data_tau$Delta_hydrophob)
ggplot(df_tau, aes(x, y, colour = hydrob, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with alzheimer's?
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, disease = main_data_tau$alzheimer)
ggplot(df_tau, aes(x, y, colour = disease, shape = df_tau$APR)) + geom_point()

#umap plot for a-beta
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()
#with expression data
#olfactory bulb
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...olfactory.region..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 60, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,60]" = "green",
                                "(60, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 60", "> 60"))
#hippocampus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 90, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,90]" = "green",
                                "(90, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 90", "> 90"))
#cerebral cortex
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 30, 300, Inf)))) +
  scale_color_manual(values = c("(-Inf,30]" = "blue",
                                "(30,300]" = "green",
                                "(300, Inf]" = "red"),
                     labels = c("<= 30", "30 < cc_exp <= 300", "> 300"))
#thalamus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...thalamus..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 80, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,80]" = "green",
                                "(80, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 80", "> 80"))
#foldX interaction data
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, interact = main_data_abeta$interaction)
ggplot(df_abeta, aes(x, y, colour = df_abeta$interact, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with logp
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, logp = main_data_abeta$logP)
ggplot(df_abeta, aes(x, y, colour = df_abeta$logp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with hydrophobicity
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, hydrop = main_data_abeta$Delta_hydrophob)
ggplot(df_abeta, aes(x, y, colour = df_abeta$hydrop, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)
#with crossint and elongation energy
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, crossint = main_data_abeta$cross_int_energy)
ggplot(df_abeta, aes(x, y, colour = crossint, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -0.5)
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, elong = main_data_abeta$elongation_energy)
ggplot(df_abeta, aes(x, y, colour = elong, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)

###FOR ASYN, GET THE INTERESTING DATAPOINTS
asyn_umap_matrix <- data.frame(asyn_umap$layout)
asyn_umap_matrix <- cbind(asyn_umap_matrix, MANGO_asyn, APR_asyn, Uniprot_asyn)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,1] > -10 && asyn_umap$layout[i,1] < -2.5) {
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

write.table(MANGO_list, file = "cluster_MANGOs_asyn_fullUMAP_capping_elong(bottom).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_asyn_fullUMAP_capping_elong(bottom).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_asyn_fullUMAP_capping_elong(bottom).txt", quote = FALSE, sep = "\n", row.names = FALSE)


####A-BETA INTERESTING CLUSTER
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,2] > -2.2 && abeta_umap$layout[i,2] < 2.5 && abeta_umap$layout[i,1] > -1.25 && abeta_umap$layout[i,1] < 2) {
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
write.table(MANGO_list, file = "cluster_MANGOs_abeta(fullUMAP_onecluster).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_abeta(fullUMAP_onecluster).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_abeta(fullUMAP_onecluster).txt", quote = FALSE, sep = "\n", row.names = FALSE)

#####GET SOD1 CLUSTER#####
sod1_umap_matrix <- data.frame(sod1_umap$layout)
sod1_umap_matrix <- cbind(sod1_umap_matrix, MANGO_sod1, APR_sod1, Uniprot_sod1)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(sod1_umap$layout)) {
  if (sod1_umap$layout[i,1] > -1.8 && sod1_umap$layout[i,2] < 0) {
    MANGO_list <- append(MANGO_list, sod1_umap_matrix$MANGO_sod1[i])
    APR_list <- append(APR_list, sod1_umap_matrix$APR_sod1[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, sod1_umap_matrix$Uniprot_sod1[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], colo = color_list)
ggplot(df_sod1, aes(x, y, colour = df_sod1$colo)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_sod1(fullUMAP).txt", quote = FALSE, sep = "\n", row.names = FALSE)
#write.table(APR_list, file = "cluster_APRs_tau(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_sod1(fullUMAP).txt", quote = FALSE, sep = "\n", row.names = FALSE)



#############UMAP with only expression values, brain_exp & interaction foldX##############
##########################################################################################
all_data_a_syn <- read.csv("a_syn_binary.tsv", sep = "\t", header = TRUE)
all_data_tau <- read.csv("tau_binary.tsv", sep = "\t", header = TRUE)
all_data_a_beta <- read.csv("a_beta_binary.tsv", sep = "\t", header = TRUE)
all_data_sod1 <- read.csv("sod_1_binary.tsv", sep = "\t", header = TRUE)

all_data_a_syn <- na.omit(all_data_a_syn)
all_data_a_beta <- na.omit(all_data_a_beta)
all_data_tau <- na.omit(all_data_tau)
all_data_sod1 <- na.omit(all_data_sod1)

aggprot_sod1 <- all_data_sod1[,1]
APR_sod1 <- all_data_sod1[,2]
Uniprot_sod1 <- all_data_sod1[,3]
MANGO_sod1 <- all_data_sod1[,4]
main_data_sod1 <- all_data_sod1[,c(5:59,61,62,77,78)]

aggprot_asyn <- all_data_a_syn[,1]
APR_asyn <- all_data_a_syn[,2]
Uniprot_asyn <- all_data_a_syn[,3]
MANGO_asyn <- all_data_a_syn[,4]
main_data_asyn <- all_data_a_syn[,c(5:59,61,62,77,78)]

aggprot_abeta <- all_data_a_beta[,1]
APR_abeta <- all_data_a_beta[,2]
Uniprot_abeta <- all_data_a_beta[,3]
MANGO_abeta <- all_data_a_beta[,4]
main_data_abeta <- all_data_a_beta[,c(5:59,61,62,77,78)]

aggprot_tau <- all_data_tau[,1]
APR_tau <- all_data_tau[,2]
Uniprot_tau <- all_data_tau[,3]
MANGO_tau <- all_data_tau[,4]
main_data_tau <- all_data_tau[,c(5:59,61,62,77,78)]

#set the binary columns as factor
cols <- c("brain_exp", "interaction")
main_data_tau[,cols] <- lapply(main_data_tau[,cols], as.factor)
main_data_abeta[,cols] <- lapply(main_data_abeta[,cols], as.factor)
main_data_asyn[,cols] <- lapply(main_data_asyn[,cols], as.factor)
main_data_sod1[,cols] <- lapply(main_data_sod1[,cols], as.factor)

#calculate gower_distance
gower_dist_sod1 <- as.matrix(daisy(main_data_sod1, metric="gower"))
gower_dist_abeta <- as.matrix(daisy(main_data_abeta, metric="gower"))
gower_dist_asyn <- as.matrix(daisy(main_data_asyn, metric="gower"))
gower_dist_tau <- as.matrix(daisy(main_data_tau, metric="gower"))

#umap
custom_umap_sod1 <- umap.defaults
custom_umap_sod1$n_neighbors <- 10
custom_umap_sod1$min_dist <- 0.1

custom_umap_asyn <- umap.defaults
custom_umap_asyn$n_neighbors <- 30
custom_umap_asyn$min_dist <- 0.1

custom_umap_tau <- umap.defaults
custom_umap_tau$n_neighbors <-16
custom_umap_tau$min_dist <- 0.1

custom_umap_abeta <- umap.defaults
custom_umap_abeta$n_neighbors <- 70
custom_umap_abeta$min_dist <- 0.1

asyn_umap <- umap(gower_dist_asyn, config = custom_umap_asyn)
tau_umap <- umap(gower_dist_tau, config = custom_umap_tau)
sod1_umap <- umap(gower_dist_sod1, config = custom_umap_sod1)
abeta_umap <- umap(gower_dist_abeta, config = custom_umap_abeta)

#umap plot for alpha synuclein
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()
#with expression data
#cerebellum
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...cerebellum..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 55, breaks = c(0, 35, 70))
#cerebral cortex
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red",midpoint = 200, breaks = c(0,150,300))

df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 90)
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...olfactory.region..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 90)
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...basal.ganglia..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 150)
#with foldx interaction data
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, interact = main_data_asyn$interaction)
ggplot(df_asyn, aes(x, y, colour = df_asyn$interact, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with brain_exp
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = MANGO_asyn, brain_dist = main_data_asyn$logP)
ggplot(df_asyn, aes(x, y, colour = df_asyn$logp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#brain exp
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, brain_exp = main_data_asyn$brain_exp)
ggplot(df_asyn, aes(x, y, colour = df_asyn$brain_exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))

#umap plot for sod1
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1)
ggplot(df_sod1, aes(x, y, colour = df_sod1$APR)) + geom_point()
#with expression value (spinal cord)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, exp = main_data_sod1$Tissue.RNA...spinal.cord..NX.)
ggplot(df_sod1, aes(x, y, colour = df_sod1$exp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 40, breaks = c(0,15, 50))
#with interaction value(foldx)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, inter = main_data_sod1$interaction)
ggplot(df_sod1, aes(x, y, colour = df_sod1$inter, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))
#logP
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = MANGO_sod1, logp = main_data_sod1$logP)
ggplot(df_sod1, aes(x, y, colour = df_sod1$logp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#delta hydrophobicity
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = MANGO_sod1, hyd = main_data_sod1$Delta_hydrophob)
ggplot(df_sod1, aes(x, y, colour = df_sod1$hyd, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)

#umap plot for tau
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()
#with expression data
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, exp = main_data_tau$Tissue.RNA...amygdala..NX.)
ggplot(df_tau, aes(x, y, colour = df_tau$exp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 55, breaks = c(0,30,70))
#with logp
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, logp = main_data_tau$logP)
ggplot(df_tau, aes(x, y, colour = df_tau$logp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#with interaction value (foldX)
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, inter = main_data_tau$interaction)
ggplot(df_tau, aes(x, y, colour = inter, shape = df_tau$APR)) + geom_point()
#with hydrophobicity
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = MANGO_tau, hydrob = main_data_tau$Delta_hydrophob)
ggplot(df_tau, aes(x, y, colour = hydrob, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)

#umap plot for a-beta
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()
#with expression data
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...olfactory.region..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 35, breaks = c(0,20,50))
#foldX interaction data
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, interact = main_data_abeta$interaction)
ggplot(df_abeta, aes(x, y, colour = df_abeta$interact, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with logp
logP <- all_data_a_beta[,c()]
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = MANGO_abeta, logp = main_data_abeta$logP)
ggplot(df_abeta, aes(x, y, colour = df_abeta$logp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with hydrophobicity
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = MANGO_abeta, hydrop = main_data_abeta$Delta_hydrophob)
ggplot(df_abeta, aes(x, y, colour = df_abeta$hydrop, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)

####A-BETA INTERESTING CLUSTER
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] > 6 && abeta_umap$layout[i,2] >= 3.2) {
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
write.table(MANGO_list, file = "cluster_MANGOs_abeta.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_abeta.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_abeta.txt", quote = FALSE, sep = "\n", row.names = FALSE)


###TAU INTERESTING CLUSTER
tau_umap_matrix <- data.frame(tau_umap$layout)
tau_umap_matrix <- cbind(tau_umap_matrix, MANGO_tau, APR_tau, Uniprot_tau)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(tau_umap$layout)) {
  if (tau_umap$layout[i,2] < 10.6 && tau_umap$layout[i,2] > 9 && tau_umap$layout[i,1] < 2.5 && tau_umap$layout[i,1] > -1.5) {
    MANGO_list <- append(MANGO_list, tau_umap_matrix$MANGO_tau[i])
    APR_list <- append(APR_list, tau_umap_matrix$APR_tau[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, tau_umap_matrix$Uniprot_tau[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = color_list)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "cluster_MANGOs_tau.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_tau.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_tau.txt", quote = FALSE, sep = "\n", row.names = FALSE)


###A-SYN INTERESTING CLUSTER
asyn_umap_matrix <- data.frame(asyn_umap$layout)
asyn_umap_matrix <- cbind(asyn_umap_matrix, MANGO_asyn, APR_asyn, Uniprot_asyn)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,1] < 0 && asyn_umap$layout[i,1] > -5 && asyn_umap$layout[i,2] < -5) {
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

write.table(MANGO_list, file = "cluster_MANGOs_asyn.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_asyn.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_asyn.txt", quote = FALSE, sep = "\n", row.names = FALSE)



#################SEPARATE UMAPS WITH BRAINEXP + BRAINDIST + FOLDX#############
##########################BBF: BRAIN, BRAIN,FOLDX#################### 
####(tau was added with the secondary structure propensity data)###
all_data_a_syn <- read.csv("a_syn_binary.tsv", sep = "\t", header = TRUE)
all_data_tau <- read.csv("tau_binary.tsv", sep = "\t", header = TRUE)
all_data_a_beta <- read.csv("a_beta_binary.tsv", sep = "\t", header = TRUE)
all_data_sod1 <- read.csv("sod1_binary_str.tsv", sep = "\t", header = TRUE)

all_data_a_syn <- na.omit(all_data_a_syn)
all_data_a_beta <- na.omit(all_data_a_beta)
all_data_tau <- na.omit(all_data_tau)
all_data_sod1 <- na.omit(all_data_sod1)

aggprot_sod1 <- all_data_sod1[,1]
APR_sod1 <- all_data_sod1[,2]
Uniprot_sod1 <- all_data_sod1[,3]
MANGO_sod1 <- all_data_sod1[,4]
BBF_data_sod1 <- all_data_sod1[,c(7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78:82)]

aggprot_asyn <- all_data_a_syn[,1]
APR_asyn <- all_data_a_syn[,2]
Uniprot_asyn <- all_data_a_syn[,3]
MANGO_asyn <- all_data_a_syn[,4]
BBF_data_asyn <- all_data_a_syn[,c(7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78)]

aggprot_abeta <- all_data_a_beta[,1]
APR_abeta <- all_data_a_beta[,2]
Uniprot_abeta <- all_data_a_beta[,3]
MANGO_abeta <- all_data_a_beta[,4]
BBF_data_abeta <- all_data_a_beta[,c(7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78)]

aggprot_tau <- all_data_tau[,1]
APR_tau <- all_data_tau[,2]
Uniprot_tau <- all_data_tau[,3]
MANGO_tau <- all_data_tau[,4]
pancreas_exp_tau <- all_data_tau[,34]
BBF_data_tau <- all_data_tau[,c(7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78)]

#set the binary columns as factor
cols <- c("brain_exp", "interaction")
BBF_data_tau[,cols] <- lapply(BBF_data_tau[,cols], as.factor)
BBF_data_abeta[,cols] <- lapply(BBF_data_abeta[,cols], as.factor)
BBF_data_asyn[,cols] <- lapply(BBF_data_asyn[,cols], as.factor)
BBF_data_sod1[,cols] <- lapply(BBF_data_sod1[,cols], as.factor)

#calculate gower_distance
BBF_dist_sod1 <- as.matrix(daisy(BBF_data_sod1, metric="gower"))
BBF_dist_abeta <- as.matrix(daisy(BBF_data_abeta, metric="gower"))
BBF_dist_asyn <- as.matrix(daisy(BBF_data_asyn, metric="gower"))
BBF_dist_tau <- as.matrix(daisy(BBF_data_tau, metric="gower"))

#umap
custom_umap_sod1 <- umap.defaults
custom_umap_sod1$n_neighbors <- 10
custom_umap_sod1$min_dist <- 0.1

custom_umap_asyn <- umap.defaults
custom_umap_asyn$n_neighbors <- 50
custom_umap_asyn$min_dist <- 0.1

custom_umap_tau <- umap.defaults
custom_umap_tau$n_neighbors <-25
custom_umap_tau$min_dist <- 0.1

custom_umap_abeta <- umap.defaults
custom_umap_abeta$n_neighbors <- 60
custom_umap_abeta$min_dist <- 0.1

asyn_umap <- umap(BBF_dist_asyn, config = custom_umap_asyn)
tau_umap <- umap(BBF_dist_tau, config = custom_umap_tau)
sod1_umap <- umap(BBF_dist_sod1, config = custom_umap_sod1)
abeta_umap <- umap(BBF_dist_abeta, config = custom_umap_abeta) #45n right now

#umap plot for alpha synuclein
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn)
ggplot(df_asyn, aes(x, y, colour = df_asyn$APR)) + geom_point()
#with expression data
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = BBF_data_asyn$Tissue.RNA...cerebellum..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 55, breaks = c(0, 35, 70))
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red",midpoint = 200, breaks = c(0,150,300))
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 90)
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...olfactory.region..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 90)
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, exp = main_data_asyn$Tissue.RNA...basal.ganglia..NX.)
ggplot(df_asyn, aes(x, y, colour = df_asyn$exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 150)
#with foldx interaction data
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, interact = BBF_data_asyn$interaction)
ggplot(df_asyn, aes(x, y, colour = df_asyn$interact, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with brain_exp
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = MANGO_asyn, brain_dist = main_data_asyn$logP)
ggplot(df_asyn, aes(x, y, colour = df_asyn$logp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#brain exp
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, brain_exp = main_data_asyn$brain_exp)
ggplot(df_asyn, aes(x, y, colour = df_asyn$brain_exp, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))

#umap plot for sod1
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1)
ggplot(df_sod1, aes(x, y, colour = df_sod1$APR)) + geom_point()
#with expression value (spinal cord)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, exp = BBF_data_sod1$Tissue.RNA...spinal.cord..NX.)
ggplot(df_sod1, aes(x, y, colour = df_sod1$exp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 40, breaks = c(0,15, 50))
ggplot(df_sod1, aes(x, y, colour = df_sod1$exp, shape = df_sod1$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 15, 70, Inf)))) +
  scale_color_manual(values = c("(-Inf,15]" = "blue",
                                "(15,70]" = "green",
                                "(70, Inf]" = "red"),
                     labels = c("<= 15", "15 < pancreas_exp <= 70", "> 70"))
#with interaction value(foldx)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, inter = BBF_data_sod1$interaction)
ggplot(df_sod1, aes(x, y, colour = df_sod1$inter, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))
#logP
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = MANGO_sod1, logp = main_data_sod1$logP)
ggplot(df_sod1, aes(x, y, colour = df_sod1$logp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#delta hydrophobicity
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, str = main_data_sod1$structure)
ggplot(df_sod1, aes(x, y, colour = str, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))
#brain_exp
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, braindist = BBF_data_sod1$brain_exp)
ggplot(df_sod1, aes(x, y, colour = df_sod1$braindist, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))

#umap plot for tau
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()
#with expression data
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, exp = BBF_data_tau$Tissue.RNA...amygdala..NX.)
ggplot(df_tau, aes(x, y, colour = df_tau$exp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 70, breaks = c(0,30,70))
ggplot(df_tau, aes(x, y, colour = df_tau$exp, shape = df_tau$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 90, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,90]" = "green",
                                "(90, Inf]" = "red"),
                     labels = c("<= 20", "20 < amygdala_exp <= 90", "> 90"))
#with logp
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, logp = main_data_tau$logP)
ggplot(df_tau, aes(x, y, colour = df_tau$logp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#with interaction value (foldX)
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, inter = BBF_data_tau$interaction)
ggplot(df_tau, aes(x, y, colour = inter, shape = df_tau$APR)) + geom_point()
#with secondary str propensity
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, str = BBF_data_tau$structure)
ggplot(df_tau, aes(x, y, colour = str, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5))
#with brain dist
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, braindist = BBF_data_tau$brain_exp)
ggplot(df_tau, aes(x, y, colour = braindist, shape = df_tau$APR)) + geom_point()
#with pancreatic expression(external value to umap)
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, exp = pancreas_exp_tau)
ggplot(df_tau, aes(x, y, colour = df_tau$exp, shape = df_tau$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 400, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,400]" = "green",
                                "(400, Inf]" = "red"),
                     labels = c("<= 20", "20 < pancreas_exp <= 400", "> 400"))


#umap plot for a-beta
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()
#with expression data
#olfactory bulb
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...olfactory.region..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, c(-Inf, 20, 60, Inf))), alpha = 0.3) +
  scale_color_manual(values = c("(-Inf,20]" = "purple",
                                "(20,60]" = "green",
                                "(60, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 60", "> 60"))
#hippocampus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 60)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 90, Inf))), alpha = 0.3) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,90]" = "green",
                                "(90, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 90", "> 90"))
#cerebral cortex
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 250)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 30, 300, Inf)))) +
  scale_color_manual(values = c("(-Inf,30]" = "blue",
                                "(30,300]" = "green",
                                "(300, Inf]" = "red"),
                     labels = c("<= 30", "30 < cc_exp <= 300", "> 300"))
#thalamus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...thalamus..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 50)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 80, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,80]" = "green",
                                "(80, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 80", "> 80"))
#foldX interaction data
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, interact = BBF_data_abeta$interaction)
ggplot(df_abeta, aes(x, y, colour = df_abeta$interact, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with logp
logP <- all_data_a_beta[,c()]
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = MANGO_abeta, logp = main_data_abeta$logP)
ggplot(df_abeta, aes(x, y, colour = df_abeta$logp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with brain_dist
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, braindist = BBF_data_abeta$brain_exp)
ggplot(df_abeta, aes(x, y, colour = braindist, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) 


#####GET SOD1 CLUSTER#####
sod1_umap_matrix <- data.frame(sod1_umap$layout)
sod1_umap_matrix <- cbind(sod1_umap_matrix, MANGO_sod1, APR_sod1, Uniprot_sod1)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(sod1_umap$layout)) {
  if (sod1_umap$layout[i,2] < 2.5 && sod1_umap$layout[i,1] < -4 && sod1_umap$layout[i,1] > -8 && sod1_umap$layout[i,2] > -5) {
    MANGO_list <- append(MANGO_list, sod1_umap_matrix$MANGO_sod1[i])
    APR_list <- append(APR_list, sod1_umap_matrix$APR_sod1[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, sod1_umap_matrix$Uniprot_sod1[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], colo = color_list)
ggplot(df_sod1, aes(x, y, colour = df_sod1$colo)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_sod1(BBF+secondarystr).txt", quote = FALSE, sep = "\n", row.names = FALSE)
#write.table(APR_list, file = "cluster_APRs_tau(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_sod1(BBF+secondarystr).txt", quote = FALSE, sep = "\n", row.names = FALSE)


######GET TAU CLUSTERS####
tau_umap_matrix <- data.frame(tau_umap$layout)
tau_umap_matrix <- cbind(tau_umap_matrix, MANGO_tau, APR_tau, Uniprot_tau)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(tau_umap$layout)) {
  if (tau_umap$layout[i,1] > 4 && tau_umap$layout[i,2] < 4.2 && tau_umap$layout[i,2] > -2.5) {
    MANGO_list <- append(MANGO_list, tau_umap_matrix$MANGO_tau[i])
    APR_list <- append(APR_list, tau_umap_matrix$APR_tau[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, tau_umap_matrix$Uniprot_tau[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = color_list)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()

write.table(MANGO_list, file = "cluster_MANGOs_tau(BBF+).txt", quote = FALSE, sep = "\n", row.names = FALSE)
#write.table(APR_list, file = "cluster_APRs_tau(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_tau(BBF+).txt", quote = FALSE, sep = "\n", row.names = FALSE)

##get the smaller tau cluster (just the very high expressions + elongationFoldX)
tau_umap_matrix <- data.frame(tau_umap$layout)
tau_umap_matrix <- cbind(tau_umap_matrix, MANGO_tau, APR_tau, Uniprot_tau)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(tau_umap$layout)) {
  if (tau_umap$layout[i,1] >= 5 && tau_umap$layout[i,1] < 7.5 && tau_umap$layout[i,2] > 0 && tau_umap$layout[i,2] < 5) {
    MANGO_list <- append(MANGO_list, tau_umap_matrix$MANGO_tau[i])
    APR_list <- append(APR_list, tau_umap_matrix$APR_tau[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, tau_umap_matrix$Uniprot_tau[i])
  } else if (tau_umap$layout[i,1] >= 4.9 && tau_umap$layout[i,1] < 7 && tau_umap$layout[i,2] >= 5) {
    MANGO_list <- append(MANGO_list, tau_umap_matrix$MANGO_tau[i])
    APR_list <- append(APR_list, tau_umap_matrix$APR_tau[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, tau_umap_matrix$Uniprot_tau[i])
  } else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = color_list)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()

#write to files
write.table(MANGO_list, file = "cluster_MANGOs_tau(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
#write.table(APR_list, file = "cluster_APRs_tau(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_tau(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)

###GET THE A-SYN CLUSTERS
asyn_umap_matrix <- data.frame(asyn_umap$layout)
asyn_umap_matrix <- cbind(asyn_umap_matrix, MANGO_asyn, APR_asyn, Uniprot_asyn)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,1] > 2.3 && asyn_umap$layout[i,2] < 7 && asyn_umap$layout[i,2] > 2.3) {
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

#write the mango list to a file
write.table(MANGO_list, file = "cluster_MANGOs_asyn(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)


######ABETA CLUSTER (2 small clusters)
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,1] > 7.5 && abeta_umap$layout[i,1] < 12.5 && abeta_umap$layout[i,2] < 7.5 && abeta_umap$layout[i,2] > 0) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } else if (abeta_umap$layout[i,2] > -0.2 && abeta_umap$layout[i,2] < 2.5 && abeta_umap$layout[i,1] < -4 && abeta_umap$layout[i,1] > -7.7) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  } 
  else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta(2smallclusters).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta(2smallclusters).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta(2smallclusters).txt", quote = FALSE, sep = "\n", row.names = FALSE)

#####A-BETA left small cluster
#for the a-beta umap, get the interesting cluster
abeta_umap_matrix <- data.frame(abeta_umap$layout)
abeta_umap_matrix <- cbind(abeta_umap_matrix, MANGO_abeta, APR_abeta, Uniprot_abeta)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(abeta_umap$layout)) {
  if (abeta_umap$layout[i,2] > -0.2 && abeta_umap$layout[i,2] < 2.5 && abeta_umap$layout[i,1] < -4 && abeta_umap$layout[i,1] > -7.7) {
    MANGO_list <- append(MANGO_list, abeta_umap_matrix$MANGO_abeta[i])
    APR_list <- append(APR_list, abeta_umap_matrix$APR_abeta[i])
    color_list <- append(color_list, "pink")
    Uniprot_list <- append(Uniprot_list, abeta_umap_matrix$Uniprot_abeta[i])
  }
  else {
    color_list <- append(color_list, "gray")
  }
}

#color the points in umap to see if we got the right ones
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = color_list)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()

#write the MANGO list to a file
write.table(MANGO_list, file = "MANGOs_BBF_abeta(leftsmallcluster).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "APRs_BBF_abeta(leftsmallcluster).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "Unip_BBF_abeta(leftsmallcluster).txt", quote = FALSE, sep = "\n", row.names = FALSE)

###################UMAP WITH BBF + TPM of AD/PD BRAIN (cerebralcortex)############
TPM_data_a_syn <- read.csv("asyn_TPM_binary.tsv", sep = "\t", header = TRUE)
TPM_data_a_syn <- na.omit(TPM_data_a_syn)

aggprot_asyn <- TPM_data_a_syn[,1]
APR_asyn <- TPM_data_a_syn[,2]
Uniprot_asyn <- TPM_data_a_syn[,3]
MANGO_asyn <- TPM_data_a_syn[,4]
TPM_data_asyn <- TPM_data_a_syn[,c(7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78,79)]

#columns as factors
cols <- c("brain_exp", "interaction")
TPM_data_asyn[,cols] <- lapply(TPM_data_asyn[,cols], as.factor)

#gowers distance
TPM_dist_asyn <- as.matrix(daisy(TPM_data_asyn, metric="gower"))

#umap
TPM_asyn_config <- custom_umap
TPM_asyn_config$n_neighbors <- 50
TPM_asyn_config$min_dist <- 0.1

TPM_asyn_UMAP <- umap(TPM_dist_asyn, config = TPM_asyn_config)

#plot
TPM_df_asyn <- data.frame(x = TPM_asyn_UMAP$layout[,1], y = TPM_asyn_UMAP$layout[,2], APR = APR_asyn)
ggplot(TPM_df_asyn, aes(x, y, colour = TPM_df_asyn$APR)) + geom_point()
#with expression data
TPM_df_asyn <- data.frame(x = TPM_asyn_UMAP$layout[,1], y = TPM_asyn_UMAP$layout[,2], APR = APR_asyn, exp = TPM_data_asyn$Tissue.RNA...cerebellum..NX.)
ggplot(TPM_df_asyn, aes(x, y, colour = TPM_df_asyn$exp, shape = TPM_df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 55, breaks = c(0, 35, 70))
#with TPM data
TPM_df_asyn <- data.frame(x = TPM_asyn_UMAP$layout[,1], y = TPM_asyn_UMAP$layout[,2], APR = APR_asyn, exp = TPM_data_asyn$TPM.cerebralcortex)
ggplot(TPM_df_asyn, aes(x, y, colour = TPM_df_asyn$exp, shape = TPM_df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 2000)

##remove the outlier (RVTTVA)
x <- which.max(TPM_data_asyn$TPM.cerebralcortex)
TPM_data_removed_asyn <- TPM_data_asyn[-c(x),]
TPM_removed_UMAP_asyn <- TPM_asyn_UMAP$layout[-c(x),]
#aggprot_removed <- TPM_data_removed_asyn$AggProt
cc_exp_removed <- TPM_data_removed_asyn$TPM.cerebralcortex
df_TPM_asyn <- data.frame(x = TPM_removed_UMAP_asyn[,1], y = TPM_removed_UMAP_asyn[,2], gradient = cc_exp_removed)
ggplot(df_TPM_asyn, aes(x, y, colour = gradient)) + geom_point() + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 700)



########UMAPS with all + Secondsry Str Propensity#############
setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")
all_data_a_syn <- read.csv("asyn_binary_str.tsv", sep = "\t", header = TRUE)
all_data_tau <- read.csv("tau_binary_str.tsv", sep = "\t", header = TRUE)
all_data_a_beta <- read.csv("abeta_binary_str.tsv", sep = "\t", header = TRUE)
all_data_sod1 <- read.csv("sod1_binary_str.tsv", sep = "\t", header = TRUE)

all_data_a_syn <- na.omit(all_data_a_syn)
all_data_a_beta <- na.omit(all_data_a_beta)
all_data_tau <- na.omit(all_data_tau)
all_data_sod1 <- na.omit(all_data_sod1)

aggprot_sod1 <- all_data_sod1[,1]
APR_sod1 <- all_data_sod1[,2]
Uniprot_sod1 <- all_data_sod1[,3]
MANGO_sod1 <- all_data_sod1[,4]
main_data_sod1 <- all_data_sod1[,-c(1,2,3,4)]

aggprot_abeta <- all_data_a_beta[,1]
APR_abeta <- all_data_a_beta[,2]
Uniprot_abeta <- all_data_a_beta[,3]
MANGO_abeta <- all_data_a_beta[,4]
main_data_abeta <- all_data_a_beta[,-c(1,2,3,4)]

aggprot_asyn <- all_data_a_syn[,1]
APR_asyn <- all_data_a_syn[,2]
Uniprot_asyn <- all_data_a_syn[,3]
MANGO_asyn <- all_data_a_syn[,4]
main_data_asyn <- all_data_a_syn[,-c(1,2,3,4)]

aggprot_tau <- all_data_tau[,1]
APR_tau <- all_data_tau[,2]
Uniprot_tau <- all_data_tau[,3]
MANGO_tau <- all_data_tau[,4]
main_data_tau <- all_data_tau[,-c(1,2,3,4)]


#set the binary columns as factor
cols <- c("alzheimer", "parkinson", "als", "brain_exp", "interaction", "structure")
main_data_tau[,cols] <- lapply(main_data_tau[,cols], as.factor)
main_data_abeta[,cols] <- lapply(main_data_abeta[,cols], as.factor)
main_data_asyn[,cols] <- lapply(main_data_asyn[,cols], as.factor)
main_data_sod1[,cols] <- lapply(main_data_sod1[,cols], as.factor)

#calculate gower_distance
gower_dist_sod1 <- as.matrix(daisy(main_data_sod1, metric="gower"))
gower_dist_abeta <- as.matrix(daisy(main_data_abeta, metric="gower"))
gower_dist_asyn <- as.matrix(daisy(main_data_asyn, metric="gower"))
gower_dist_tau <- as.matrix(daisy(main_data_tau, metric="gower"))

#umap
custom_umap_sod1 <- umap.defaults
custom_umap_sod1$n_neighbors <- 25
custom_umap_sod1$min_dist <- 0.1

custom_umap_asyn <- umap.defaults
custom_umap_asyn$n_neighbors <- 50
custom_umap_asyn$min_dist <- 0.1

custom_umap_tau <- umap.defaults
custom_umap_tau$n_neighbors <- 40
custom_umap_tau$min_dist <- 0.1

custom_umap_abeta <- umap.defaults
custom_umap_abeta$n_neighbors <- 60
custom_umap_abeta$min_dist <- 0.1

asyn_umap <- umap(gower_dist_asyn, config = custom_umap_asyn)
tau_umap <- umap(gower_dist_tau, config = custom_umap_tau)
sod1_umap <- umap(gower_dist_sod1, config = custom_umap_sod1)
abeta_umap <- umap(gower_dist_abeta, config = custom_umap_abeta)

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
#structure
df_asyn <- data.frame(x = asyn_umap$layout[,1], y = asyn_umap$layout[,2], APR = APR_asyn, str = main_data_asyn$structure)
ggplot(df_asyn, aes(x, y, colour = str, shape = df_asyn$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))

#umap plot for sod1
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1)
ggplot(df_sod1, aes(x, y, colour = df_sod1$APR)) + geom_point()
#with expression value (spinal cord)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, exp = main_data_sod1$Tissue.RNA...spinal.cord..NX.)
ggplot(df_sod1, aes(x, y, colour = df_sod1$exp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 40, breaks = c(0,15, 50))
#with interaction value(foldx)
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, inter = main_data_sod1$interaction)
ggplot(df_sod1, aes(x, y, colour = df_sod1$inter, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))
#logP
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, logp = main_data_sod1$logP)
ggplot(df_sod1, aes(x, y, colour = df_sod1$logp, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#delta hydrophobicity
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, hyd = main_data_sod1$Delta_hydrophob)
ggplot(df_sod1, aes(x, y, colour = df_sod1$hyd, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)
#structure
df_sod1 <- data.frame(x = sod1_umap$layout[,1], y = sod1_umap$layout[,2], APR = APR_sod1, str = main_data_sod1$structure)
ggplot(df_sod1, aes(x, y, colour = str, shape = df_sod1$APR)) + geom_point() + scale_shape_manual(values = c(1,5))

#umap plot for tau
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau)
ggplot(df_tau, aes(x, y, colour = df_tau$APR)) + geom_point()
#with expression data
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, exp = main_data_tau$Tissue.RNA...amygdala..NX.)
ggplot(df_tau, aes(x, y, colour = df_tau$exp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 70)
#with logp
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, logp = main_data_tau$logP)
ggplot(df_tau, aes(x, y, colour = df_tau$logp, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -1)
#with interaction value (foldX)
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, inter = main_data_tau$interaction)
ggplot(df_tau, aes(x, y, colour = inter, shape = df_tau$APR)) + geom_point()
#with hydrophobicity
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, hydrob = main_data_tau$Delta_hydrophob)
ggplot(df_tau, aes(x, y, colour = hydrob, shape = df_tau$APR)) + geom_point() + scale_shape_manual(values = c(1,5)) + scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with structure
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, str = main_data_tau$structure)
ggplot(df_tau, aes(x, y, colour = str, shape = df_tau$APR)) + geom_point()


#umap plot for a-beta
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta)
ggplot(df_abeta, aes(x, y, colour = df_abeta$APR)) + geom_point()
#with expression data
#olfactory bulb
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...olfactory.region..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 60, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,60]" = "green",
                                "(60, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 60", "> 60"))
#hippocampus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...hippocampal.formation..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 20, 90, Inf)))) +
  scale_color_manual(values = c("(-Inf,20]" = "blue",
                                "(20,90]" = "green",
                                "(90, Inf]" = "red"),
                     labels = c("<= 20", "20 < cc_exp <= 90", "> 90"))
#cerebral cortex
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...cerebral.cortex..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 30, 300, Inf)))) +
  scale_color_manual(values = c("(-Inf,30]" = "blue",
                                "(30,300]" = "green",
                                "(300, Inf]" = "red"),
                     labels = c("<= 30", "30 < cc_exp <= 300", "> 300"))
#thalamus
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, exp = main_data_abeta$Tissue.RNA...thalamus..NX.)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 45)
ggplot(df_abeta, aes(x, y, colour = df_abeta$exp, shape = df_abeta$APR)) +
  geom_point(aes(colour = cut(exp, alpha = 0.3, c(-Inf, 25, 80, Inf)))) +
  scale_color_manual(values = c("(-Inf,25]" = "blue",
                                "(25,80]" = "green",
                                "(80, Inf]" = "red"),
                     labels = c("<= 25", "25 < cc_exp <= 80", "> 80"))
#foldX interaction data
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, interact = main_data_abeta$interaction)
ggplot(df_abeta, aes(x, y, colour = df_abeta$interact, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))
#with logp
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, logp = main_data_abeta$logP)
ggplot(df_abeta, aes(x, y, colour = df_abeta$logp, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0)
#with hydrophobicity
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, hydrop = main_data_abeta$Delta_hydrophob)
ggplot(df_abeta, aes(x, y, colour = df_abeta$hydrop, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)
#with crossint and elongation energy
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, crossint = main_data_abeta$cross_int_energy)
ggplot(df_abeta, aes(x, y, colour = crossint, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = -0.5)
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, elong = main_data_abeta$elongation_energy)
ggplot(df_abeta, aes(x, y, colour = elong, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5)) +scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0.5)
#with structure
df_abeta <- data.frame(x = abeta_umap$layout[,1], y = abeta_umap$layout[,2], APR = APR_abeta, str = main_data_abeta$structure)
ggplot(df_abeta, aes(x, y, colour = str, shape = df_abeta$APR)) + geom_point() + scale_shape_manual(values = c(1,2,5))


########A-SYN CLUSTER
asyn_umap_matrix <- data.frame(asyn_umap$layout)
asyn_umap_matrix <- cbind(asyn_umap_matrix, MANGO_asyn, APR_asyn, Uniprot_asyn)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(asyn_umap$layout)) {
  if (asyn_umap$layout[i,1] < -0.1 && asyn_umap$layout[i,2] < 2 && asyn_umap$layout[i,2] > -3.5) {
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

#write to a file
write.table(APR_list, file = "cluster_APRs_asyn_propensity.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(MANGO_list, file = "cluster_MANGOs_asyn_propensity.txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Uniprot_asyn_propensity.txt", quote = FALSE, sep = "\n", row.names = FALSE)
