rm(list=ls())

#library(DMwR)
library(GGally)
library(ggplot2)
library(pROC)
library(caret)
library(RColorBrewer)
library(pheatmap)
library(factoextra)
library(reshape2)
library(tidyverse)
library(pls)
library(ggplot2)
library(Boruta)
library(VIM)
library(ggpubr)
library(mice)
library(visdat)
library(naniar)

## Set data and results dirs
data.dir <- "/home/user/project/data/"
setwd(data.dir)

results.dir <- "/home/user/project/results/"
results.tag <- paste(tag, ".res_", resampling, ".boruta_", do.boruta, sep="")

## PARAMETERS TO CUSTOMIZE
args <- commandArgs(trailingOnly = TRUE)
tag <- args[1]
resampling <- args[2]
do.boruta <- args[3]
k <- as.numeric(args[4]) ## k-cross fold validation
subsampling.times <- as.numeric(args[5])

#k <- 10 ## k-cross fold validation
times.the.size <- 1 # number of times that should multiply the cases to make the train/test datasets
seed.sum <- 1
lucky.number <- 101
percentage.test <- 0.8
black.and.white <- "no"

tag.classes <- c(rep("numeric", 10))

tag.mice.methods <- "pmm"


tag.table <- read.table(paste(tag, ".csv", sep=""), sep="\t", header=T, na.strings = "NA", 
                              stringsAsFactors = F, colClasses = tag.classes)
matrix.frame <- tag.table

## set patients names
rownames(tag.table) <- c(1:dim(tag.table)[1])

#treatment.table <- read.table("treatment.csv", sep="\t", header=T, na.strings = "NA", stringsAsFactors = F, colClasses = "numeric")
#relapse.table <- read.table("relapse.csv", sep="\t", header=T, na.strings = "NA", stringsAsFactors = F, colClasses = "numeric")

# Find patients that have many missing data
how_many_missing <- 25
patients <- rownames(matrix.frame)
patients.to.remove <- c()

for(i in 1:length(patients)){
  missing <- sum(is.na(matrix.frame[i,]))
  total <- length(matrix.frame[i,])
  perc_missing <- (missing*100)/total
  if(perc_missing > how_many_missing){
    patients.to.remove <- c(patients.to.remove, i)
  }
}
if(length(patients.to.remove)>0){
  matrix.frame <- matrix.frame[-patients.to.remove,]
}

## Find columns that have >how_many_missing% of missing data
variables <- colnames(matrix.frame)
to.remove <- c()
for(i in 1:length(variables)){
  variables[i]
  missing <- sum(is.na(matrix.frame[,i]))
  missing
  total <- length(matrix.frame[,i])
  perc_missing <- (missing*100)/total
  perc_missing
  if(perc_missing > how_many_missing){
    to.remove <- c(to.remove, i)
  }
}
## Remove - only if there is anything to remove...
if(length(to.remove)>0){
  matrix.frame <- matrix.frame[,-to.remove]
  if(length(tag.mice.methods)>1){
    tag.mice.methods.filtered <- tag.mice.methods[-to.remove]
  }
  else{
    tag.mice.methods.filtered <- tag.mice.methods 
  }
} else{
  tag.mice.methods.filtered <- tag.mice.methods
}

## Write removed variables
variables.removed.file <- paste(data.dir, results.tag, ".variables_removed.txt", sep="")
write.table(variables[to.remove], file=variables.removed.file, quote=F)

# IMPUTATION
# get plot with missing percentages
missing.plot.file <- paste(results.dir, "missing.", results.tag, ".png", sep="")
png(file=missing.plot.file, width=500, height = 800)
gg_miss_var(matrix.frame)
dev.off()

# now impute data
tempData <- mice(matrix.frame,m=1,maxit=100,meth=tag.mice.methods.filtered,seed=666)
summary(tempData)

#plot(tempData)

imputation.plot.file <- paste(results.dir, "imputation.", results.tag, ".png", sep="")
png(file=imputation.plot.file, width=600, height = 800)
stripplot(tempData)
dev.off()

imputed.matrix <- complete(tempData)

## We remove the columns that have no variability
variables <- colnames(imputed.matrix)
no.variables <- c()
for(i in 1:length(variables)){
  if(length(unique(imputed.matrix[,i])) < 2){
    no.variables <- c(no.variables, i)
  }
}
## Remove - only if there is anything to remove...
if(length(no.variables)>0){
  imputed.matrix <- imputed.matrix[,-no.variables]
}

## NOW WE FIX THE DATA AND THE MATRIX SORTING OUTPUT 1 AND 0 AS THE PROGRAM REQUIRES
fix.imputed.matrix <- rbind(imputed.matrix[imputed.matrix$OUTPUT==1,], 
                            imputed.matrix[imputed.matrix$OUTPUT==0,])

fix.matrix.file <- paste(data.dir, results.tag, ".fix.matrix.txt", sep="")
write.table(fix.imputed.matrix, file=fix.matrix.file, sep="\t", row.names=F, quote=F)

### PCA ###

## Leave OUTPUT out
 ## Same matrix but with imputed values
pca.imputed.matrix <- complete(tempData)[,-length(names(fix.imputed.matrix))] ## without last column

## converted to numeric
pca.imputed.matrix <- sapply(pca.imputed.matrix, as.numeric)

## We remove the columns that have no variability, NOW in the PCA matrix
variables <- colnames(pca.imputed.matrix)
no.variables <- c()
for(i in 1:length(variables)){
  if(length(unique(pca.imputed.matrix[,i])) < 2){
    no.variables <- c(no.variables, i)
  }
}
## Remove - only if there is anything to remove...
if(length(no.variables)>0){
  pca.imputed.matrix <- pca.imputed.matrix[,-no.variables]
}


tag.pca <- prcomp(pca.imputed.matrix, scale = TRUE) ## and the pca

groups <- as.factor(fix.imputed.matrix$OUTPUT) ## maybe change the name of the groups




pca.plot.file <- paste(results.dir, "pca.", results.tag, ".png", sep="")
png(file=pca.plot.file, width=500, height = 500)

fviz_pca_ind(tag.pca,
             axes = c(1,2),
             col.ind = groups, # color by groups
             geom="point",
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "OUTPUT",
             repel = TRUE) + scale_shape_manual(values=c(19,19,19,19))
dev.off()

# fviz_pca_var(tag.pca,
#              axes = c(1,2),
#              repel = TRUE)


## Heatmap ##
dist.samples <- dist(pca.imputed.matrix, method="euclidean")
hclust.samples <- hclust(dist.samples, method="average")

heatmap.plot.file <- paste(results.dir, "heatmap.", results.tag, ".png", sep="")
png(file=heatmap.plot.file, width=900, height = 900)
my_heatmap <- pheatmap(pca.imputed.matrix, cluster_rows=hclust.samples)
dev.off()

## Boruta ##
if(do.boruta == "yes"){
  set.seed(111)
  boruta.bank_train <- Boruta(OUTPUT~., data = fix.imputed.matrix, doTrace = 2, maxRuns=100)
  print(boruta.bank_train)

  boruta.bank <- TentativeRoughFix(boruta.bank_train)
  print(boruta.bank)

  ## Plot Boruta features selection
  boruta.plot.file <- paste(results.dir, "boruta.features.importance.", results.tag, ".png", sep="")
  png(file=boruta.plot.file, width=900, height = 800)
  plot(boruta.bank, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(boruta.bank$ImpHistory),function(i)
            boruta.bank$ImpHistory[is.finite(boruta.bank$ImpHistory[,i]),i])
  names(lz) <- colnames(boruta.bank$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),
        at = 1:ncol(boruta.bank$ImpHistory), cex.axis = 0.7)
  dev.off()

  getSelectedAttributes(boruta.bank, withTentative = F)
  bank_df <- attStats(boruta.bank)
  print(bank_df)


  ## Get the important features
  selected.matrix.file <- paste(data.dir, "features.selected.", results.tag ,".txt", sep="")
  last.selected.frame <- fix.imputed.matrix[,c(boruta.bank$finalDecision == "Confirmed", TRUE)]
  write.table(last.selected.frame, file=selected.matrix.file, sep="\t", quote=F)
} else{
  selected.matrix.file <- paste(data.dir, "features.selected.", results.tag ,".txt", sep="")
  last.selected.frame <- fix.imputed.matrix
  write.table(last.selected.frame, file=selected.matrix.file, sep="\t", quote=F)
}

## Subsampling j times
number.to.keep <- sum(last.selected.frame$OUTPUT==1)*times.the.size
if(resampling=="no"){
  number.to.keep <- sum(last.selected.frame$OUTPUT==0)
  subsampling.times <- 1
}

## Vectors that save each performance to calculate mean afterwards
auc.vect <- rep(NA,k*subsampling.times) ## store the error in this vector (AUC)
accuracy.vect <- rep(NA, k*subsampling.times) ## store the accuracy in this vector
sensitivity.vector <- rep(NA, k*subsampling.times)
specificity.vector <- rep(NA, k*subsampling.times)

## Importance
importance.frame <- data.frame(Parameter=character(0),value=numeric(0))

count <- 0

## Draw the first part of the roc curve
if(black.and.white == "yes"){
  my.colors <- rep("grey", 10)
}else{
  my.colors <- c("red", "deepskyblue", "yellow", "green4", "hotpink", "orange", "blueviolet", "navy", "sienna", "darksalmon")
}

## an object to save the AUC plot
png(filename=paste(results.dir, "auc.", results.tag, ".results.png", sep=""), width=500, height=480)

# We repeat the analysis if we go for subsampling
for (j in 1:subsampling.times){
  matrix.frame <- last.selected.frame
  
  first.value <- sum(matrix.frame$OUTPUT==1)+1
  last.value <- dim(matrix.frame)[1]
  rows <- c(1:sum(matrix.frame$OUTPUT==1), sample(first.value:last.value, number.to.keep, replace=F))
  matrix.frame <- matrix.frame[rows,]
  
  ## We record vectors to draw a ROC curve per each resampling
  classes.local.vector <- vector(mode="logical")
  predictions.local.vector <- vector(mode="logical")
  
  ## 10-fold cross validation loop
  for(i in 1:k){
    set.seed(i + seed.sum + (j*100))
    count <- count + 1
    
    index <- sample(1:nrow(matrix.frame),round(percentage.test*nrow(matrix.frame)))
    
    ## creating training and test sets
    train.nn = matrix.frame[index, ]
    test.nn = matrix.frame[-index, ]
    
    ## To avoid 
    while(length(unique(test.nn$OUTPUT)) == 1 || length(unique(train.nn$OUTPUT)) == 1){
      set.seed(i + seed.sum + (j*100) + lucky.number)
      index <- sample(1:nrow(matrix.frame),round(0.9*nrow(matrix.frame)))
      
      ## creating training and test sets
      train.nn = matrix.frame[index, ]
      test.nn = matrix.frame[-index, ]
      
      lucky.number <- lucky.number + 1
    }
    
    ## Convert response variable to factor (MORTALITY) to avoid NN does regression
    train.nn$OUTPUT <- as.character(train.nn$OUTPUT)
    train.nn$OUTPUT <- as.factor(train.nn$OUTPUT)
    test.nn$OUTPUT <- as.character(test.nn$OUTPUT)
    test.nn$OUTPUT <- as.factor(test.nn$OUTPUT)
    
    ## Get the model and prediciton
    crtl <- trainControl(method="none",classProbs = T, summaryFunction=twoClassSummary)
    my.grid <- expand.grid(.decay = seq(from = 0.1, to = 0.5, by = 0.1), .size = c(1,2,3,4,5))
    nn.caret <- train(OUTPUT~., data=train.nn, method="nnet", trainControl=crtl, tuneGrid=my.grid)
    
    #plotnet(nn.caret)
    predict_test.nn.prob <- predict(nn.caret, test.nn, type = "prob")
    predict_test.nn.class <- predict(nn.caret, test.nn)
    
    t <- table(predictions=predict_test.nn.class, real=test.nn$OUTPUT) ## predictions vs real
    t
    
    results <- data.frame(actual = test.nn$OUTPUT, prediction = predict_test.nn.prob$"0")
    
    ## Accuracy
    accuracy.vect[count] <- sum(diag(t))/sum(t) ## percentage accuracy
    
    ## ROC element
    my.roc <- roc(predictor=predict_test.nn.prob$"0", response=results$actual)
    ## AUC calculated using all classes, other methods only allow two classes
    auc.vect[count] <- my.roc$auc[1]
    
    ## Calculate sensitivity & specificity
    specificity.vector[count] <- mean(mean(my.roc$specificities))
    sensitivity.vector[count] <- mean(mean(my.roc$sensitivities))
    
    classes.local.vector <- c(classes.local.vector, as.numeric(levels(results$actual))[results$actual])
    predictions.local.vector <- c(predictions.local.vector, results$prediction)
    
    ## Importance 
    for(h in 1:length(colnames(matrix.frame)[1:length(colnames(matrix.frame))-1])){
      measurement <- colnames(matrix.frame)[h]
      importance.frame <- rbind(importance.frame, data.frame(Parameter=measurement, 
                                                             value=as.vector(varImp(nn.caret)$importance[h,1])))
    }
  }
  
  ## Generate the ROC for a re-sampled matrix
  intermediate.roc <- roc(factor(classes.local.vector), predictions.local.vector, ci=T)
  
  if(j==1){
    plot.roc(intermediate.roc, col=my.colors[j], add=F, lwd=3, cex=1.2)
  }else{plot.roc(intermediate.roc, col=my.colors[j], add=T, lwd=3, cex=1.2)}
}

#dev.off()
while (!is.null(dev.list())){dev.off()}

## Now deal with other paramenters (as before adding the ROC curve)
print(paste("Average AUC:", mean(auc.vect)))
print(paste("Average Accuracy:", mean(accuracy.vect)))
print(paste("Average Sensitivity:", mean(sensitivity.vector)))
print(paste("Average Specificity:", mean(specificity.vector)))

## Print results to file
results.frame <- data.frame(AUC=numeric(0),
                            ACCURACY=numeric(0), 
                            SENSITIVITY=numeric(0), 
                            SPECIFICITY=numeric(0))

results.frame <- rbind(results.frame, data.frame(AUC=round(mean(auc.vect), digits=2), 
                                                 ACCURACY=round(mean(accuracy.vect), digits=2),
                                                 SENSITIVITY=round(mean(sensitivity.vector), digits=2),
                                                 SPECIFICITY=round(mean(specificity.vector), digits=2)))

results.file <- paste(results.dir, results.tag, ".results.txt", sep="")
write.table(t(results.frame), file=results.file,
            sep="\t", quote=F, col.names = F)


## Plot Importance parameters
levels(importance.frame$Parameter)

## This reorders the legend
importance.reorder <- reorder(importance.frame$Parameter, importance.frame$value, FUN=mean)
levels.reordered <- rev(attributes(importance.reorder)$levels)

importance.frame$Parameter <- factor(importance.frame$Parameter, levels = levels.reordered)

## Calculate significance of differences and plot the importance

fix.imputed.matrix <- read.table(selected.matrix.file, sep="\t", header=T, stringsAsFactors = F)

matrix.good <- sapply(fix.imputed.matrix[fix.imputed.matrix$OUTPUT == 0,], as.numeric)
matrix.bad <- sapply(fix.imputed.matrix[fix.imputed.matrix$OUTPUT == 1,], as.numeric)

p.vector <- c()
direction.vector <- c()
median.direction.vector <- c()

for(i in 1:length(levels.reordered)){
  #colors <- c()
  levels.reordered[i]
  my.index <- match(levels.reordered[i], colnames(fix.imputed.matrix))
  
  if(length(unique(fix.imputed.matrix[,my.index]))>2){
    ## Non-binary variable
    good.vector <- matrix.good[,my.index]
    bad.vector <- matrix.bad[,my.index]
    
    w <- wilcox.test(good.vector, bad.vector)
    p.value <- w$p.value
    my.ratio <- mean(good.vector)/mean(bad.vector)
    my.median <- median(good.vector)/median(bad.vector)
    p.vector <- c(p.vector, p.value)
    direction.vector <- c(direction.vector, my.ratio)
    median.direction.vector <- c(median.direction.vector, my.median)
  }
  else{
    ## Binary variable
    good.feature <- sum(matrix.good[,my.index])
    good.nofeature <- abs(length(matrix.good[,my.index]) - good.feature)
    
    bad.feature <- sum(matrix.bad[,my.index])
    bad.nofeature <- abs(length(matrix.bad[,my.index]) - bad.feature)
    
    contingency.matrix <- matrix(c(good.feature, bad.feature, good.nofeature, bad.nofeature), nrow = 2,
                                 dimnames = list(c("GOOD", "BAD"), c("feature","NO_feature")))
    f <- fisher.test(contingency.matrix)
    p.value <- f$p.value
    odds.ratio <- f$estimate
    p.vector <- c(p.vector, p.value)
    direction.vector <- c(direction.vector, odds.ratio)
    median.direction.vector <- c(median.direction.vector, odds.ratio)
  }
}

adjusted <- p.adjust(p.vector, method = "fdr")
significance.frame <- data.frame(direction.vector, median.direction.vector, p.vector, adjusted)
rownames(significance.frame) <- levels.reordered

## Print significance results
results.file <- paste(results.dir, ".", results.tag, ".variables.significance.txt", sep="")
write.table(significance.frame, file=results.file,
            sep="\t", quote=F, col.names = T)

## And choose the colors
colors.vector <- c()

for(i in 1:length(levels.reordered)){
  new.color <- "#999999"
  if(significance.frame$direction.vector[i]>1){
    new.color <- "#9999FF"
    if(significance.frame$adjusted[i]<0.05){
      new.color <- "#0000FF"
    }
  }
  else if(significance.frame$direction.vector[i]<1){
    new.color <- "#FF9999"
    if(significance.frame$adjusted[i]<0.05){
      new.color <- "#FF0000"
    }
  }
  colors.vector <- c(colors.vector,new.color)
}

## And print the parameters
p1<-ggplot(importance.frame, aes(x=reorder(Parameter, value, FUN=mean), y=value)) +
  geom_boxplot(fill=rev(colors.vector), color=rev(colors.vector)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
               panel.background = element_rect(fill="white", size=1, color="grey", linetype = "solid"),
               panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey90")) + 
  scale_y_continuous("Importance", limits=c(min(importance.frame$value), 
                      max(importance.frame$value))) + 
  scale_x_discrete("Parameters")

dat <- ggplot_build(p1)$data[[1]]
p1 + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                               y=middle, yend=middle), colour="grey90", size=0.8)
p1.file <- paste(results.dir, "importance.", results.tag, ".results.png", sep="")
ggsave(p1.file, width=8, height = 4)
#dev.off()
while (!is.null(dev.list())){dev.off()}

## Now with the median
## And choose the colors
median.colors.vector <- c()

for(i in 1:length(levels.reordered)){
  new.color <- "#999999"
  if(significance.frame$adjusted[i]<0.05){
    new.color <- "#666666"
  }
  
  if(is.na(significance.frame$median.direction.vector[i])){
    new.color <- "#ABABAB"
  }
  else if(significance.frame$median.direction.vector[i]>1){
    new.color <- "#9999FF"
    if(significance.frame$adjusted[i]<0.05){
      new.color <- "#0000FF"
    }
  }
  else if(significance.frame$median.direction.vector[i]<1){
    new.color <- "#FF9999"
    if(significance.frame$adjusted[i]<0.05){
      new.color <- "#FF0000"
    }
  }
  median.colors.vector <- c(median.colors.vector, new.color)
}

p2<-ggplot(importance.frame, aes(x=reorder(Parameter, value, FUN=mean), y=value)) +
  geom_boxplot(fill=rev(median.colors.vector), color=rev(median.colors.vector)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8),
        panel.background = element_rect(fill="white", size=1, color="grey", linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, linetype = "solid", colour = "grey90")) + 
  scale_y_continuous("Importance", limits=c(min(importance.frame$value), 
                                            max(importance.frame$value))) + 
  scale_x_discrete("Parameters")

dat <- ggplot_build(p2)$data[[1]]
p2 + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                y=middle, yend=middle), colour="grey90", size=0.8)

p2.file <- paste(results.dir, "importance.", results.tag, ".median.results.png", sep="")
ggsave(p2.file, width=8, height = 4)
#dev.off()
while (!is.null(dev.list())){dev.off()}


## Calculate significance of differences and plot the importance
matrix.good <- fix.imputed.matrix[fix.imputed.matrix$OUTPUT == 0,]
matrix.bad <- fix.imputed.matrix[fix.imputed.matrix$OUTPUT == 1,]
p.vector <- c()
direction.vector <- c()
alternative.direction.vector <- c()

for(i in 1:length(colnames(fix.imputed.matrix)[-length(colnames(fix.imputed.matrix))])){
  if(length(unique(fix.imputed.matrix[,i]))>2){
    ## Non-binary variable
    good.vector <- matrix.good[,i]
    bad.vector <- matrix.bad[,i]
    
    w <- wilcox.test(good.vector, bad.vector)
    p.value <- w$p.value
    my.ratio <- mean(good.vector)/mean(bad.vector)
    my.median <- median(good.vector)/median(bad.vector)
    p.vector <- c(p.vector, p.value)
    direction.vector <- c(direction.vector, my.ratio)
    alternative.direction.vector <- c(alternative.direction.vector, my.median)
  }
  else{
    ## Binary variable
    good.feature <- sum(matrix.good[,i])
    good.nofeature <- abs(length(matrix.good[,i]) - good.feature)
    
    bad.feature <- sum(matrix.bad[,i])
    bad.nofeature <- abs(length(matrix.bad[,i]) - bad.feature)
    
    contingency.matrix <- matrix(c(good.feature, bad.feature, good.nofeature, bad.nofeature), nrow = 2,
                                 dimnames = list(c("GOOD", "BAD"), c("feature","NO_feature")))
    f <- fisher.test(contingency.matrix)
    p.value <- f$p.value
    odds.ratio <- f$estimate
    
    relative.odds.ratio <- odds.ratio
    if(odds.ratio<1){
      relative.odds.ratio <- 1/odds.ratio
    }
    
    p.vector <- c(p.vector, p.value)
    direction.vector <- c(direction.vector, odds.ratio)
    alternative.direction.vector <- c(alternative.direction.vector, relative.odds.ratio)
  }
}

adjusted <- p.adjust(p.vector, method = "fdr")
significance.frame <- data.frame(direction.vector, alternative.direction.vector, p.vector, adjusted)
rownames(significance.frame) <- colnames(fix.imputed.matrix)[-length(colnames(fix.imputed.matrix))]


## Print enrichment results
significance.frame$features <- row.names(significance.frame)
significance.frame$log_pvals <- -log(significance.frame$adjusted)

features_sorted <- arrange(significance.frame, alternative.direction.vector)
significance.frame$features <- factor(significance.frame$features, levels=features_sorted$features)

features.colors <- c()
for(i in 1:length(rownames(significance.frame))){
  new.color <- "#999999"
  if(features_sorted$direction.vector[i]>1){
    new.color <- "#9999FF"
    if(features_sorted$adjusted[i]<0.05){
      new.color <- "#0000FF"
    }
  }
  else if(features_sorted$direction.vector[i]<1){
    new.color <- "#FF9999"
    if(features_sorted$adjusted[i]<0.05){
      new.color <- "#FF0000"
    }
  }
  features.colors <- c(features.colors,new.color)
}

enrichment.tag <- "wilcoxon"
enrichment.tag2 <- "relative mean"

if(length(unique(matrix.frame[,1]))==2){
  enrichment.tag <- "fisher"
  enrichment.tag2 <- "relative odds ratio"
  
}

features.histogram <- ggplot(significance.frame, aes(features, alternative.direction.vector, fill=features)) + 
  geom_bar(stat="identity", width=0.5) + 
  theme(panel.background = element_rect(fill='white', colour="black"), legend.position="none", 
        axis.text=element_text(size=10)) + 
  scale_fill_manual(values=features.colors) +
  labs(x=results.tag, y=enrichment.tag2, fill="odds ratio") + coord_flip() #+
#geom_hline(yintercept=2, linetype="dashed", color = "grey50")

features.histogram
fisher.file <- paste(results.dir, enrichment.tag, ".", results.tag, ".results.png", sep="")
ggsave(fisher.file, width=4, height = 8)
#dev.off()
while (!is.null(dev.list())){dev.off()}
#}


