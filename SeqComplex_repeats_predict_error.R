library(plotly)
library(data.table)
library(stringr)
library(dplyr)

setwd("/Users/pinaronat/Desktop/Thesis Work/protein_frq")
#get the frequencies for aminoacids
aa_freq <- read.table("aa_freq.txt", sep = "-")
aa_freq <- transpose(aa_freq)
colnms <- aa_freq[1,]
colnames(aa_freq) <- colnms
aa_freq <- aa_freq[-1,]

setwd("/Users/pinaronat/Desktop/Thesis Work/protein_frq/snippets")

#get points to predict (APRs)
APR_matches = read.table("APR_matches_new.txt", sep = "\t", header = TRUE)
APR_matches = APR_matches[,c(2,4)]


#add complexities, frequencies and colors for APRs
APR_matches <- cbind(APR_matches, data.frame(matrix(nrow = 13, ncol = 1)))
APR_matches <- cbind(APR_matches, data.frame(matrix(nrow = 13, ncol = 1)))
colnames(APR_matches) <- c("sequence", "matches", "complexity", "frequency")
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
  APR_matches$complexity[i] = complexity
  #  APR_matches$color[i] = color
  
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
  #add frequency to the dataframe
  APR_matches$frequency[i] = frequency
}

#create dataframe to put standard errors in
se_data <- data.frame()
apr_se <- data.frame()

setwd("/Users/pinaronat/Desktop/Thesis Work/protein_frq/snippets/repeats")
#for each repeat (generated sequence set)
for (i in 0:99) {
  #get the matches file
  filename = capture.output(cat("PSX_Global_", i, ".txt", sep = ""))
  matches = read.table(filename, sep = "\t", header = TRUE)
  total = matches$X1mut + matches$X2mut + matches$identical
  nr_matches = cbind(matches$baitsequence, total) #create sequence-nr of matches dataframe
  #add 2 columns to nr_matches for complexity and frequency
  nr_matches <- cbind(nr_matches, data.frame(matrix(nrow = 1500, ncol = 1)))
  nr_matches <- cbind(nr_matches, data.frame(matrix(nrow = 1500, ncol = 1)))
  colnames(nr_matches) <- c("sequence", "matches", "complexity", "frequency")
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
  }
  ####fit regressions
  #separate dataframe based on complexity
  separate <- split(nr_matches, nr_matches$complexity)
  complex_3 <- separate$`3`
  complex_4 <- separate$`4`
  complex_5 <- separate$`5`
  reg3 <- lm(matches ~ poly(as.numeric(frequency),2), data = complex_3)
  reg4 <- lm(matches ~ poly(as.numeric(frequency),2), data = complex_4)
  reg5 <- lm(matches ~ poly(as.numeric(frequency),2), data = complex_5)
  temp_var3 <- predict(reg3, interval="prediction", se.fit = TRUE)
  temp_var4 <- predict(reg4, interval="prediction", se.fit = TRUE)
  temp_var5 <- predict(reg5, interval="prediction", se.fit = TRUE)
  #get the predicted values
  fits_3 <- temp_var3$fit[,1]
  fits_4 <- temp_var4$fit[,1]
  fits_5 <- temp_var5$fit[,1]
  #get the difference between observed-predicted
  complex_3$matches <- as.numeric(complex_3$matches)
  complex_4$matches <- as.numeric(complex_4$matches)
  complex_5$matches <- as.numeric(complex_5$matches)
  diff3 <- complex_3$matches - fits_3
  diff4 <- complex_4$matches - fits_4
  diff5 <- complex_5$matches - fits_5
  label3 <- capture.output(cat(i+1, "-", 3, sep = ""))
  label4 <- capture.output(cat(i+1, "-", 4, sep = ""))
  label5 <- capture.output(cat(i+1, "-", 5, sep = ""))
  se_data3 <- cbind(c(diff3), rep(label3, length(temp_var3$se.fit)), rep("3-complexity", length(temp_var3$se.fit)), temp_var3$fit[,2], temp_var3$fit[,3], rep("black", length(temp_var3$se.fit)))
  se_data4 <- cbind(c(diff4), rep(label4, length(temp_var4$se.fit)), rep("4-complexity", length(temp_var4$se.fit)), temp_var4$fit[,2], temp_var4$fit[,3], rep("black", length(temp_var3$se.fit)))
  se_data5 <- cbind(c(diff5), rep(label5, length(temp_var5$se.fit)), rep("5-complexity", length(temp_var5$se.fit)), temp_var5$fit[,2], temp_var5$fit[,3], rep("black", length(temp_var3$se.fit)))
  se_data <- rbind(se_data, se_data3, se_data4, se_data5)
  #for APRs
  separate_apr <- split(APR_matches, APR_matches$complexity)
  apr_3 <- separate_apr$`3`
  apr_4 <- separate_apr$`4`
  apr_5 <- separate_apr$`5`
  apr_var3 <- predict(reg3, newdata = apr_3, interval="prediction", se.fit = TRUE)
  apr_var4 <- predict(reg4, newdata = apr_4, interval="prediction", se.fit = TRUE)
  apr_var5 <- predict(reg5, newdata = apr_5, interval="prediction", se.fit = TRUE)
  #get the predicted values
  apr_fits_3 <- apr_var3$fit[,1]
  apr_fits_4 <- apr_var4$fit[,1]
  apr_fits_5 <- apr_var5$fit[,1]
  #get the difference between observed-predicted
  apr_3$matches <- as.numeric(apr_3$matches)
  apr_4$matches <- as.numeric(apr_4$matches)
  apr_5$matches <- as.numeric(apr_5$matches)
  apr_diff3 <- apr_3$matches - apr_fits_3
  apr_diff4 <- apr_4$matches - apr_fits_4
  apr_diff5 <- apr_5$matches - apr_fits_5
  apr_se3 <- cbind(c(apr_diff3), rep(label3, length(apr_var3$se.fit)), rep("3-complexity", length(apr_var3$se.fit)), apr_var3$fit[,2], apr_var3$fit[,3], rep("red2", length(apr_var3$se.fit)))
  apr_se4 <- cbind(c(apr_diff4), rep(label4, length(apr_var4$se.fit)), rep("4-complexity", length(apr_var4$se.fit)), apr_var4$fit[,2], apr_var4$fit[,3], rep("green4", length(apr_var4$se.fit)))
  apr_se5 <- cbind(c(apr_diff5), rep(label5, length(apr_var5$se.fit)), rep("5-complexity", length(apr_var5$se.fit)), apr_var5$fit[,2], apr_var5$fit[,3], rep("mediumblue", length(apr_var5$se.fit)))
  apr_se <- rbind(apr_se, apr_se3, apr_se4, apr_se5)
}
colnames(se_data) <- c("prediction errors", "repeat-complexity", "complexity", "lwr", "upper", "color")
colnames(apr_se) <- c("prediction errors", "repeat-complexity", "complexity", "lwr", "upper", "color")

se_data <- as.data.frame(se_data)
se_data$`prediction errors` <- as.numeric(se_data$`prediction errors`)
se_data$`repeat-complexity` <- as.character(se_data$`repeat-complexity`)
se_data$`complexity` <- as.factor(se_data$complexity)
se_data$`repeat-complexity` <- factor(se_data$`repeat-complexity`, levels = unique(se_data$`repeat-complexity`))

total_se <- rbind(se_data, apr_se)
total_se$`repeat-complexity` <- factor(total_se$`repeat-complexity`, levels = unique(total_se$`repeat-complexity`))
total_se$`prediction errors` <- as.numeric(total_se$`prediction errors`)
total_se <- total_se[order(total_se$`repeat-complexity`),]

#separate se_data into four for better visualization
first_end <- as.numeric(max(which(se_data$`repeat-complexity` == "25-5")))
second_end <- as.numeric(max(which(se_data$`repeat-complexity` == "50-5")))
third_end <- as.numeric(max(which(se_data$`repeat-complexity` == "75-5")))
fourth_end <- as.numeric(max(which(se_data$`repeat-complexity` == "100-5")))


#for errorbars - standard error
errbar_lims <- group_by(se_data, `repeat-complexity`) %>% 
  summarize(mean=mean(`prediction errors`), se=sd(`prediction errors`)/sqrt(n()), 
            upper=mean+(2*se), lower=mean-(2*se))
errbar_lims$`repeat-complexity` <- factor(errbar_lims$`repeat-complexity`, levels = unique(errbar_lims$`repeat-complexity`))
errorbar1 <- as.numeric(max(which(errbar_lims$`repeat-complexity` == "25-5")))
errorbar2 <-as.numeric(max(which(errbar_lims$`repeat-complexity` == "50-5")))
errorbar3 <-as.numeric(max(which(errbar_lims$`repeat-complexity` == "75-5")))
errorbar4 <-as.numeric(max(which(errbar_lims$`repeat-complexity` == "100-5")))
errbar_lims1 <- errbar_lims[c(1:75),]
errbar_lims2 <- errbar_lims[c(76:150),]
errbar_lims3 <- errbar_lims[c(151:225),]
errbar_lims4 <- errbar_lims[c(226:300),]


#with jitter - just aprs
apr_se$`repeat-complexity` <- factor(apr_se$`repeat-complexity`, levels = unique(apr_se$`repeat-complexity`))
apr_se <- apr_se[order(apr_se$`repeat-complexity`),]
aprs_limit1 <- as.numeric(max(which(apr_se$`repeat-complexity` == "25-5")))
apr_se$`prediction errors` <- as.numeric(apr_se$`prediction errors`)
#plot 1
ggplot(se_data[c(1:37501),], aes(x = `repeat-complexity`, y = `prediction errors`, fill = complexity)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + geom_jitter(data = apr_se[c(1:325),], size = 0.3, color = apr_se$color[c(1:325)]) + geom_errorbar(data = errbar_lims1, aes(x=`repeat-complexity`, ymax=upper, 
                                                                                                                                                                                                                                                                                     ymin=lower), color = "purple", stat='identity', width=1, inherit.aes = FALSE)
aprs_limit2 <- as.numeric(max(which(apr_se$`repeat-complexity` == "50-5")))
#plot 2
ggplot(se_data[c(37502:75001),], aes(x = `repeat-complexity`, y = `prediction errors`, fill = complexity)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + geom_jitter(data = apr_se[c(326:650),], size = 0.3, color = apr_se$color[c(326:650)]) + geom_errorbar(data = errbar_lims2, aes(x=`repeat-complexity`, ymax=upper, 
                                                                                                                                                                                                                                                                                           ymin=lower), color = "purple", stat='identity', width=1, inherit.aes = FALSE)
aprs_limit3 <- as.numeric(max(which(apr_se$`repeat-complexity` == "75-5")))
#plot 3
ggplot(se_data[c(75002:112501),], aes(x = `repeat-complexity`, y = `prediction errors`, fill = complexity)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + geom_jitter(data = apr_se[c(651:975),], size = 0.3, color = apr_se$color[c(651:975)]) + geom_errorbar(data = errbar_lims3, aes(x=`repeat-complexity`, ymax=upper, 
                                                                                                                                                                                                                                                                                                    ymin=lower), color = "purple", stat='identity', width=1, inherit.aes = FALSE)

aprs_limit4 <- as.numeric(max(which(apr_se$`repeat-complexity` == "100-5")))
#plot 4
ggplot(se_data[c(112502:150001),], aes(x = `repeat-complexity`, y = `prediction errors`, fill = complexity)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, hjust=1)) + geom_jitter(data = apr_se[c(976:1300),], size = 0.3, color = apr_se$color[c(976:1300)]) + geom_errorbar(data = errbar_lims4, aes(x=`repeat-complexity`, ymax=upper, 
                                                                                                                                                                                                                                                                                                                        ymin=lower), color = "purple", stat='identity', width=1, inherit.aes = FALSE)


#####plot last 3 complexity as a test
apr_var3 <- as.data.frame(apr_var3)
temp_var3 <- as.data.frame(temp_var3)
temp_var <- rbind(temp_var3, apr_var3)
dummy_complex3 <- rbind(complex_3, apr_3)
new_df <- cbind(dummy_complex3, temp_var)
new_df$matches <- as.numeric(new_df$matches)
ggplot(new_df, aes(as.numeric(frequency), matches))+
  geom_point() +
  xlab("Frequency") + ylab("# of Matches")+
  geom_line(aes(y=fit.lwr), color = "red", linetype = "dashed")+
  geom_line(aes(y=fit.upr), color = "red", linetype = "dashed")+
  geom_smooth(method="lm", se=TRUE, formula = y ~ poly(x, 2))
