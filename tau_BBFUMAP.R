setwd("/Users/pinaronat/Dropbox (SwitchLab)/Nikos_Pinar/agawrapper_all/agadir_wrapper/build/Release/Scripts/Separate_UMAPs")
all_data_tau <- read.csv("tau_binary_str.tsv", sep = "\t", header = TRUE)

main_data_tau <- all_data_tau[,c(1,2,3,4,7,9,12,13,16,25,26,31,32,36,38,48,52,61,62,77,78)]
main_data_tau <- na.omit(main_data_tau)
aggprot_tau <- main_data_tau[,1]
APR_tau <- main_data_tau[,2]
Uniprot_tau <- main_data_tau[,3]
MANGO_tau <- main_data_tau[,4]
main_data_tau <- main_data_tau[,-c(1,2,3,4)]

cols <- c("brain_exp", "interaction")
main_data_tau[,cols] <- lapply(main_data_tau[,cols], as.factor)

gower_dist_tau <- as.matrix(daisy(main_data_tau, metric="gower"))

custom_umap_tau <- umap.defaults
custom_umap_tau$n_neighbors <- 25
custom_umap_tau$min_dist <- 0.1

tau_umap <- umap(gower_dist_tau, config = custom_umap_tau)


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
#with braindist
df_tau <- data.frame(x = tau_umap$layout[,1], y = tau_umap$layout[,2], APR = APR_tau, braindist = main_data_tau$brain_exp)
ggplot(df_tau, aes(x, y, colour = braindist, shape = df_tau$APR)) + geom_point()


######GET TAU CLUSTERS####
#main cluster, cap_elong_self
tau_umap_matrix <- data.frame(tau_umap$layout)
tau_umap_matrix <- cbind(tau_umap_matrix, MANGO_tau, APR_tau, Uniprot_tau)

#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(tau_umap$layout)) {
  if (tau_umap$layout[i,1] > 5.8 && tau_umap$layout[i,2] > -4) {
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

write.table(MANGO_list, file = "cluster_MANGOs_taucluster1(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_taucluster1(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_taucluster1(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)


#cluster2
#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(tau_umap$layout)) {
  if (tau_umap$layout[i,1] <= 5.8 && tau_umap$layout[i,2] > 0) {
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

write.table(MANGO_list, file = "cluster_MANGOs_taucluster2(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_taucluster2(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_taucluster2(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)



#cluster4
#get the MANGOs that are in the right coordinates
MANGO_list <- vector(mode="character")
color_list <- c()
APR_list <- c()
Uniprot_list <- c()
for (i in 1:nrow(tau_umap$layout)) {
  if (tau_umap$layout[i,1] > 0 && tau_umap$layout[i,2] < -8) {
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

write.table(MANGO_list, file = "cluster_MANGOs_taucluster4(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(APR_list, file = "cluster_APRs_taucluster4(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
write.table(Uniprot_list, file = "cluster_Unip_taucluster4(BBF).txt", quote = FALSE, sep = "\n", row.names = FALSE)
