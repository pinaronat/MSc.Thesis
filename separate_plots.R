setwd("/Users/pinaronat/Desktop/Thesis Work/protein_frq/snippets/repeats")
i = 5 # select one of the repeat files
filename = capture.output(cat("PSX_Global_", i, ".txt", sep = "")) #get the matches filename
matches = read.table(filename, sep = "\t", header = TRUE)
total = matches$X1mut + matches$X2mut + matches$identical
nr_matches = cbind(matches$baitsequence, total) #create sequence-nr of matches dataframe
#add 3 columns to nr_matches for complexity, color and frequency
nr_matches <- cbind(nr_matches, data.frame(matrix(nrow = 1500, ncol = 1)))
nr_matches <- cbind(nr_matches, data.frame(matrix(nrow = 1500, ncol = 1)))
nr_matches <- cbind(nr_matches, data.frame(matrix(nrow = 1500, ncol = 1)))
colnames(nr_matches) <- c("sequence", "matches", "complexity", "frequency", "color")
#for each sequence in the ith file
for (j in 1:nrow(nr_matches)) {
  sequence = nr_matches[j,1]
  #get complexity
  complexity = length(unique(strsplit(sequence, "")[[1]])) 
  
  #get frequency
  res_list = strsplit(sequence, "")
  freq1 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][1])]
  freq2 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][2])]
  freq3 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][3])]
  freq4 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][4])]
  freq5 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][5])]
  freq6 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][6])]
  frequency = as.numeric(freq1) + as.numeric(freq2) + as.numeric(freq3) + as.numeric(freq4) + as.numeric(freq5) + as.numeric(freq6)
  #add complexity and frequency to the dataframe
  nr_matches$complexity[j] = complexity
  nr_matches$frequency[j] = frequency
  if (complexity == 3) {
    nr_matches$color[j] = "coral1"
  } else if (complexity == 4) {
    nr_matches$color[j] = "springgreen"
  } else if (complexity == 5) {
    nr_matches$color[j] = "cornflowerblue"
  }
}


#get the aprs
setwd("/Users/pinaronat/Desktop/Thesis Work/protein_frq/snippets")
APR_matches = read.table("APR_matches_new.txt", sep = "\t", header = TRUE)
APR_matches = APR_matches[,c(2,4)]
APR_matches <- cbind(APR_matches, data.frame(matrix(nrow = 13, ncol = 1)))
APR_matches <- cbind(APR_matches, data.frame(matrix(nrow = 13, ncol = 1)))
APR_matches <- cbind(APR_matches, data.frame(matrix(nrow = 13, ncol = 1)))
colnames(APR_matches) <- c("sequence", "matches", "complexity", "frequency", "color")
for (i in 1:nrow(APR_matches)) {
  sequence = APR_matches[i,1]
  #get complexity
  complexity = length(unique(strsplit(sequence, "")[[1]]))
  #  if (complexity == 3) {
  #    color = "coral1"
  #  } else if (complexity == 4) {
  #    color = "springgreen"
  #  } else if (complexity == 5) {
  #    color = "cornflowerblue"
  #  }
  #get frequency
  res_list = strsplit(sequence, "")
  freq1 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][1])]
  freq2 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][2])]
  freq3 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][3])]
  freq4 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][4])]
  freq5 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][5])]
  freq6 = aa_freq[,which(colnames(aa_freq) == res_list[[1]][6])]
  frequency = as.numeric(freq1) + as.numeric(freq2) + as.numeric(freq3) + as.numeric(freq4) + as.numeric(freq5) + as.numeric(freq6)
  APR_matches$complexity[i] = complexity
  APR_matches$frequency[i] = frequency
  APR_matches$color[i] = "black"
}
apr_separate <- split(APR_matches, APR_matches$complexity)
APR_3 <- apr_separate$`3`
APR_4 <- apr_separate$`4`
APR_5 <- apr_separate$`5`


####fit regressions
#separate dataframe based on complexity
separate <- split(nr_matches, nr_matches$complexity)
complex_3 <- separate$`3`
complex_4 <- separate$`4`
complex_5 <- separate$`5`
complex_3$frequency <- as.numeric(complex_3$frequency)
complex_4$frequency <- as.numeric(complex_4$frequency)
complex_5$frequency <- as.numeric(complex_5$frequency)
complex_3$matches <- as.numeric(complex_3$matches)
complex_4$matches <- as.numeric(complex_4$matches)
complex_5$matches <- as.numeric(complex_5$matches)
reg3 <- lm(matches ~ poly(frequency,2), data = complex_3)
reg4 <- lm(matches ~ poly(frequency,2), data = complex_4)
reg5 <- lm(matches ~ poly(frequency,2), data = complex_5)
temp_var3 <- predict(reg3, interval="prediction", se.fit = TRUE)
temp_var4 <- predict(reg4, interval="prediction", se.fit = TRUE)
temp_var5 <- predict(reg5, interval="prediction", se.fit = TRUE)

##3-complexity
#prediction interval
total_data3 <- rbind(APR_3, complex_3)
total_data3$frequency <- as.numeric(total_data3$frequency)
total_data3$matches <- as.numeric(total_data3$mathces)
new_df3 <- cbind(complex_3, temp_var3)

c3_plot <- ggplot(total_data3, aes(as.numeric(frequency), matches))+
  geom_point(color = total_data3$color) +
  xlab("Frequency") + ylab("# of Matches")+
  geom_line(data = new_df3, aes(y=fit.lwr), color = "red", linetype = "dashed")+
  geom_line(data = new_df3, aes(y=fit.upr), color = "red", linetype = "dashed")+
  geom_smooth(data = complex_3, method = "lm", se = TRUE, formula = y ~ poly(x, 2), color = "darkgrey")

print(c3_plot + ggtitle("3-Complexity Hexapeptides"))


##4-complexity
total_data4 <- rbind(APR_4, complex_4)
total_data4$frequency <- as.numeric(total_data4$frequency)
total_data4$matches <- as.numeric(total_data4$mathces)
new_df4 <- cbind(complex_4, temp_var4)
ggplot(total_data4, aes(as.numeric(frequency), matches))+
  geom_point(color = total_data4$color) +
  xlab("Frequency") + ylab("# of Matches")+
  geom_line(data = new_df4, aes(y=fit.lwr), color = "red", linetype = "dashed")+
  geom_line(data = new_df4, aes(y=fit.upr), color = "red", linetype = "dashed")+
  geom_smooth(data = complex_4, method="lm", formula = y ~ poly(x, 2), se=TRUE, color = "darkgray")+
  ggtitle("4-complexity Hexapeptides")

##5-complexity
total_data5 <- rbind(APR_5, complex_5)
total_data5$frequency <- as.numeric(total_data5$frequency)
total_data5$matches <- as.numeric(total_data5$mathces)
new_df5 <- cbind(complex_5, temp_var5)
ggplot(total_data5, aes(as.numeric(frequency), matches))+
  geom_point(color = total_data5$color) +
  xlab("Frequency") + ylab("# of Matches")+
  geom_line(data = new_df5, aes(y=fit.lwr), color = "red", linetype = "dashed")+
  geom_line(data = new_df5, aes(y=fit.upr), color = "red", linetype = "dashed")+
  geom_smooth(data = complex_5, method="lm", formula = y ~ poly(x, 2), se=TRUE, color = "darkgray")+
  ggtitle("5-complexity Hexapeptides")
