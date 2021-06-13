install.packages("data.table")
install.packages("plotly")
library(plotly)
library(data.table)

setwd("/Users/pinaronat/Desktop/protein_frq")
#get the frequencies for aminoacids
aa_freq <- read.table("aa_freq.txt", sep = "-")
aa_freq <- transpose(aa_freq)
colnms <- aa_freq[1,]
colnames(aa_freq) <- colnms
aa_freq <- aa_freq[-1,]
#for (i in ncol(aa_freq)) {
#  aa_freq[1,i] = as.numeric(aa_freq[1,i])
#}
aa_freq[1,] <- as.data.frame(sapply(aa_freq[1,], as.numeric))

#get the random sequences
rand_1 <- read.table("rand_sequences_1.txt")
rand_2 <- read.table("rand_sequences_2.txt")
rand_3 <- read.table("rand_sequences_3.txt")
rand_4 <- read.table("rand_sequences_4.txt")
rand_5 <- read.table("rand_sequences_5.txt")
rand_6 <- read.table("rand_sequences_6.txt")

rand_1 <- cbind(rand_1, data.frame(matrix(data = rep(1,20))))
colnames(rand_1) <- c("seq", "complexity")
rand_2 <- cbind(rand_2, data.frame(matrix(data = rep(2,100))))
colnames(rand_2) <- c("seq", "complexity")
rand_3 <- cbind(rand_3, data.frame(matrix(data = rep(3,100))))
colnames(rand_3) <- c("seq", "complexity")
rand_4 <- cbind(rand_4, data.frame(matrix(data = rep(4,100))))
colnames(rand_4) <- c("seq", "complexity")
rand_5 <- cbind(rand_5, data.frame(matrix(data = rep(5,100))))
colnames(rand_5) <- c("seq", "complexity")
rand_6 <- cbind(rand_6, data.frame(matrix(data = rep(6,100))))
colnames(rand_6) <- c("seq", "complexity")
random <- rbind(rand_1, rand_2, rand_3, rand_4, rand_5, rand_6)
#add an empty column to random for the color
random <- cbind(random, data.frame(matrix(nrow = 520, ncol = 1)))
#add another empty column to random for the frequency
random <- cbind(random, data.frame(matrix(nrow = 520, ncol = 1)))
#add another amoty column for nr of matches
random <- cbind(random, data.frame(matrix(nrow = 520, ncol = 1)))
#add column for legend labels
random <- cbind(random, data.frame(matrix(nrow = 520, ncol = 1)))
colnames(random) <- c("sequence", "complexity", "color", "frequency", "matches", "label")

#add the APRs to the dataframe
random <- rbind(random, c("EGVLYV", 5, NA, NA, NA, NA))
random <- rbind(random, c("VQIVYK", 5, NA, NA, NA, NA))
random <- rbind(random, c("VQIINK", 5, NA, NA, NA, NA))
random <- rbind(random, c("GAIIGL", 4, NA, NA, NA, NA))
random <- rbind(random, c("KLVFFA", 5, NA, NA, NA, NA))
random <- rbind(random, c("MVGGVV", 3, NA, NA, NA, NA))
random <- rbind(random, c("GVATVA", 4, NA, NA, NA, NA))
random <- rbind(random, c("DSVISL", 5, NA, NA, NA, NA))
random <- rbind(random, c("SVISLS", 4, NA, NA, NA, NA))
random <- rbind(random, c("GVIGIA", 4, NA, NA, NA, NA))
random <- rbind(random, c("VIGIAQ", 5, NA, NA, NA, NA))
random <- rbind(random, c("GAVVTG", 4, NA, NA, NA, NA))
random <- rbind(random, c("AVVTGV", 4, NA, NA, NA, NA))
random <- rbind(random, c("VVTGVT", 3, NA, NA, NA, NA))
random <- rbind(random, c("VTGVTA", 4, NA, NA, NA, NA))
random <- rbind(random, c("TGVTAV", 4, NA, NA, NA, NA))
random <- rbind(random, c("GVTAVA", 4, NA, NA, NA, NA))










#get the pepsearch matches
matches <- read.table("GSIAAA_nr_of_matches.txt", sep = "\t", header = TRUE)
nr_matches <- matches[,c(2,4)]

for (i in 1:nrow(random)) {
  #add the colors
  if (random[i,]$complexity == 1) {
    random[i,]$color = "firebrick"
    random[i,]$label = "1-aa-type"
  } else if (random[i,]$complexity == 2) {
    random[i,]$color = "deepskyblue"
    random[i,]$label = "2-aa-types"
  } else if (random[i,]$complexity == 3) {
    random[i,]$color = "darkolivegreen1"
    random[i,]$label = "3-aa-types"
  } else if (random[i,]$complexity == 4) {
    random[i,]$color = "darkgoldenrod1"
    random[i,]$label = "4-aa-types"
  } else if (random[i,]$complexity == 5) {
    random[i,]$color = "darksalmon"
    random[i,]$label = "5-aa-types"
  } else if (random[i,]$complexity == 6) {
    random[i,]$color = "blueviolet"
    random[i,]$label = "6-aa-types"
  }
  if (random[i,]$sequence == "EGVLYV" || random[i,]$sequence == "VQIVYK"
  || random[i,]$sequence == "VQIINK" || random[i,]$sequence == "GAIIGL"
  || random[i,]$sequence == "GVATVA" || random[i,]$sequence == "MVGGVV"
  || random[i,]$sequence == "KLVFFA" || random[i,]$sequence == "GAVVTG"
  || random[i,]$sequence == "AVVTGV" || random[i,]$sequence == "VVTGVT" 
  || random[i,]$sequence == "VTGVTA" || random[i,]$sequence == "TGVTAV"
  || random[i,]$sequence == "GVTAVA" || random[i,]$sequence == "GVIGIA"
  || random[i,]$sequence == "VIGIAQ" || random[i,]$sequence == "DSVISL"
  || random[i,]$sequence == "SVISLS") {
    random[i,]$color = "black"
    random[i,]$label = "hexa-APR"
  }
  #add the frequencies
  sequence = random[i,]$sequence
  freq = 0
  for (j in 1:nchar(sequence)) {
    residue = substr(sequence, j, j)
    freq = freq + as.numeric(aa_freq[,which(colnames(aa_freq) == residue)])
  }
  random[i,]$frequency = freq
    
  #nr of matches
  if (sequence == nr_matches[i,]$sequence.window) {
    random$matches = nr_matches$total.matches
  }
}

#plot the outcome
fig <- plot_ly(x = random$complexity, y = random$frequency, z = random$matches, 
        color = random$label, colors = random$color, type = "scatter3d", size = 3)
fig <- fig %>% layout(scene = list(xaxis = list(title = "Complexity"), yaxis = list(title = "Frequency"), zaxis = list(title = "# of Matches"))) 
fig
