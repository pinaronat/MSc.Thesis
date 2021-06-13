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

setwd("/Users/pinaronat/Desktop/protein_frq/with_nonamyloid")
#get the random sequences
rand_2 <- read.table("2complexity_nonamyloid.txt")
rand_3 <- read.table("3complexity_nonamyloid.txt")
rand_4 <- read.table("4complexity_nonamyloid.txt")
rand_5 <- read.table("5complexity_nonamyloid.txt")
rand_6 <- read.table("6complexity_nonamyloid.txt")

rand_2 <- cbind(rand_2, data.frame(matrix(data = rep(2,24))))
colnames(rand_2) <- c("seq", "complexity")
rand_3 <- cbind(rand_3, data.frame(matrix(data = rep(3,111))))
colnames(rand_3) <- c("seq", "complexity")
rand_4 <- cbind(rand_4, data.frame(matrix(data = rep(4,238))))
colnames(rand_4) <- c("seq", "complexity")
rand_5 <- cbind(rand_5, data.frame(matrix(data = rep(5,341))))
colnames(rand_5) <- c("seq", "complexity")
rand_6 <- cbind(rand_6, data.frame(matrix(data = rep(6,187))))
colnames(rand_6) <- c("seq", "complexity")
random <- rbind(rand_2, rand_3, rand_4, rand_5, rand_6)
#add an empty column to random for the color
random <- cbind(random, data.frame(matrix(nrow = 901, ncol = 1)))
#add another empty column to random for the frequency
random <- cbind(random, data.frame(matrix(nrow = 901, ncol = 1)))
#add another amoty column for nr of matches
random <- cbind(random, data.frame(matrix(nrow = 901, ncol = 1)))
#add column for legend labels
random <- cbind(random, data.frame(matrix(nrow = 901, ncol = 1)))
colnames(random) <- c("sequence", "complexity", "color", "frequency", "matches", "label")

uniq_random <- unique(random)

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

uniq_matches <- unique(nr_matches)

for (i in 1:nrow(uniq_random)) {
  #add the colors
  if (uniq_random[i,]$complexity == 2) {
    uniq_random[i,]$color = "deepskyblue"
    uniq_random[i,]$label = "2-aa-types"
  } else if (uniq_random[i,]$complexity == 3) {
    uniq_random[i,]$color = "darkolivegreen1"
    uniq_random[i,]$label = "3-aa-types"
  } else if (uniq_random[i,]$complexity == 4) {
    uniq_random[i,]$color = "deeppink"
    uniq_random[i,]$label = "4-aa-types"
  } else if (uniq_random[i,]$complexity == 5) {
    uniq_random[i,]$color = "darksalmon"
    uniq_random[i,]$label = "5-aa-types"
  } else if (uniq_random[i,]$complexity == 6) {
    uniq_random[i,]$color = "blueviolet"
    uniq_random[i,]$label = "6-aa-types"
  }
  if (uniq_random[i,]$sequence == "EGVLYV" || uniq_random[i,]$sequence == "VQIVYK"
  || uniq_random[i,]$sequence == "VQIINK" || uniq_random[i,]$sequence == "GAIIGL"
  || uniq_random[i,]$sequence == "GVATVA" || uniq_random[i,]$sequence == "MVGGVV"
  || uniq_random[i,]$sequence == "KLVFFA" || uniq_random[i,]$sequence == "GAVVTG"
  || uniq_random[i,]$sequence == "AVVTGV" || uniq_random[i,]$sequence == "VVTGVT" 
  || uniq_random[i,]$sequence == "VTGVTA" || uniq_random[i,]$sequence == "TGVTAV"
  || uniq_random[i,]$sequence == "GVTAVA" || uniq_random[i,]$sequence == "GVIGIA"
  || uniq_random[i,]$sequence == "VIGIAQ" || uniq_random[i,]$sequence == "DSVISL"
  || uniq_random[i,]$sequence == "SVISLS") {
    uniq_random[i,]$color = "black"
    uniq_random[i,]$label = "hexa-APR"
  }
  #add the frequencies
  sequence = uniq_random[i,]$sequence
  freq = 0
  for (j in 1:nchar(sequence)) {
    residue = substr(sequence, j, j)
    freq = freq + as.numeric(aa_freq[,which(colnames(aa_freq) == residue)])
  }
  uniq_random[i,]$frequency = freq
    
  #nr of matches
  if (sequence == uniq_matches[i,]$sequence.window) {
    uniq_random$matches = uniq_matches$total.matches
  }
}

#plot the outcome
fig <- plot_ly(uniq_random, x = ~complexity, y = ~frequency, z = ~matches, 
               colors = c("deepskyblue", "darkolivegreen3", "darkgoldenrod1", "darksalmon", "blueviolet", "black") ,color = ~label, type = "scatter3d", size = 3)
fig <- fig %>% layout(scene = list(xaxis = list(title = "Complexity"), yaxis = list(title = "Frequency"), zaxis = list(title = "# of Matches"))) 
fig
