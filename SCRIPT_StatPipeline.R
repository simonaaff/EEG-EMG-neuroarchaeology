############# ------> CHANNEL SELECTION <------ #############

# Load the datasets ("Hold" condition from task1 and task2)
Hhold<-"EEG+EMG_Hhold.csv"
Fhold<-"EEG+EMG_Fhold.csv"

task1_t<-read.table(Hhold,sep=";",dec=",",header = T)
task1<-task1_t[,-1]
task1

task2_t<-read.table(Fhold,sep=";",dec=",",header = T)
task2<-task2_t[,-1]
task2

# Extract the first column to store the participant IDs for later use
nameID_1<-task1_t[,1]
nameID_2<-task2_t[,1]
combined_ID<-c(nameID_1,nameID_2)

# Add a new column to indicate the group each observation belongs to
task1$group <- "Hhold"
task2$group <- "Fhold"
task1
task2

# Remove specific channels that are not required for analysis
task1<-task1[,-which(colnames(task1)%in%c("FP1","TP9","TP10","FT9","FT10","FP2","FCz"))]
task2<-task2[,-which(colnames(task2)%in%c("FP1","TP9","TP10","FT9","FT10","FP2","FCz"))]
task1
task2

# Combine the two datasets into a single dataframe  
combined_data <- rbind(task1, task2)
str(combined_data)
combined_data





### OPTION 1. WILCOXON SIGNED-RANK TEST 

task1_clean <- task1[, !(colnames(task1) %in% "group")]
task2_clean <- task2[, !(colnames(task2) %in% "group")]
WilcPvalue <- numeric(ncol(task1_clean))
WilcZscore <- numeric(ncol(task1_clean))

# Perform Wilcoxon Signed-Rank Test for each channel
for(i in 1:ncol(task1_clean)){
  testi<-wilcox.test(task1_clean[,i],task2_clean[,i],paired = T,exact = F)
  WilcZscore[i]<- qnorm(testi$p.value/2)
  WilcPvalue[i]<-testi[[3]]
}

# Create a table with Wilcoxon results
names(WilcPvalue)<-colnames(task1_clean)
names(WilcZscore)<-colnames(task1_clean)
Wilc_table<-data.frame("Z-score"=WilcZscore,"p-value"=WilcPvalue)
write.table(Wilc_table,"Wilc_table.csv",sep=";",dec=",",col.names=NA) #Save results to a CSV file

# Extract channels with significant p-values (??? 0.05)
sel<-which(Wilc_table$p.value<=0.05)
colnames(task1)[sel]

####Select channels for PCA
# Select the top N channels based on absolute Z.score values
howmany<-11
selected_vars<-order(Wilc_table$Z.score)[1:howmany]
selected_vars
rownames(Wilc_table)[selected_vars]

# Manually select specific channels for analysis based on your predefined criteria (e.g., 5 EEG+5 EMG based on absolute Z.scores)
colnames(task1)
selected_vars<-c("FCU","FCR","FPL","HTE","DI1","O2","CP1","FC1","Cz","F3") 
selected_vars





### OPTION 2. LDA-BASED PERMUTATION TEST

combined_data_numeric <- combined_data[, sapply(combined_data, is.numeric)] #Select only numeric columns
howmanychannels<-5 #Define the number of channels to select in each permutation
combinchannels<-combn(colnames(combined_data_numeric),5,simplify = T)
dim(combn(c(1:ncol(combined_data_numeric)),5,simplify = T)) #Max number of possible permutations

permutations<-201376 #Define the number of permutations to perform

Combinchannels<-combinchannels[,sample(ncol(combinchannels),permutations)]
accmatrix<-matrix(NA,ncol=length(colnames(combined_data_numeric)),nrow=permutations)
accs<-NULL

# Perform Linear Discriminant Analysis (LDA) on each permutation
library(MASS)

for(i in 1:permutations){
  print(i)
  selcha<-Combinchannels[,i]
  permi<-data.frame(combined_data_numeric[,selcha],"group"=combined_data$group)
  lda_analisi <- lda(group ~ ., data = permi)
  lda_predictions<-predict(lda_analisi)
  actual_classes <- combined_data$group
  predicted_classes <- lda_predictions$class
  accuracy <- mean(predicted_classes == actual_classes) * 100
  accmatrix[i,which(colnames(combined_data_numeric)%in%selcha)]<-accuracy
}

# Compute the mean accuracy for each channel
table_accs<-data.frame("accuracy"=colMeans(accmatrix,na.rm = T),"channels"=colnames(combined_data_numeric))

# Rank channels based on accuracy
channels_name<-table_accs$channels[order(table_accs$accuracy,decreasing = T)]
acc_values<-table_accs$accuracy[order(table_accs$accuracy,decreasing = T)]
table<-data.frame(channels_name,acc_values)
table
write.table(table,file = paste("LDA_permut",".csv",sep="_"), row.names=TRUE, col.names = NA, sep=";", dec=",") #Save results to a CSV file

###Select channels for PCA
# Select the top N channels based on accuracy
howmany<-10
selected_vars<-table_accs$channels[order(table_accs$accuracy,decreasing = TRUE)][1:howmany]
selected_vars

# Manually select specific channels for analysis based on your predefined criteria (e.g., 5 EEG+5 EMG with highest accuracy values)
selected_vars<-c("FCR","HTE","FCU","FPL","DI1","O2","Oz","CP6","F3","CP5")
selected_vars






############ ---------> PCA <----------- ##################

dataPCA<-rbind(task1[,selected_vars],task2[,selected_vars])
meanPCA<-(task1[,selected_vars]+task2[,selected_vars])/2
dataPCA<-rbind(task1[,selected_vars]-meanPCA,task2[,selected_vars]-meanPCA) #Center the data 
taskvar<-c(rep("orange2",dim(task1)[1]),rep("red4",dim(task2)[1])) #Assign colors to different groups for visualization

# Perform PCA
pcaCor<-princomp(dataPCA,cor = T,scores=T)
eigv <- pcaCor$sdev^2
Variance <- cbind(sqrt(eigv), eigv/sum(eigv), cumsum(eigv)/sum(eigv)) * 100


# Define which principal components to plot
PCx<-1
PCy<-2

# Plot PCA results
plot(pcaCor$scores[,c(PCx,PCy)],pch=19,col=taskvar,asp=1,cex=1,
     xlab=paste("PC",PCx, " (",round(Variance[PCx,2],2),"%)",sep=""),
     ylab=paste("PC",PCy, " (",round(Variance[PCy,2],2),"%)",sep=""),
     cex.axis=0.75,cex.lab=0.8)

conv1<-chull(pcaCor$scores[taskvar=="orange2",c(PCx,PCy)])
conv1<-c(conv1,conv1[1])
polygon(pcaCor$scores[taskvar=="orange2",c(PCx,PCy)][conv1,],border=FALSE,col=adjustcolor("yellow",alpha=0.2),lwd=0.5)
conv2<-chull(pcaCor$scores[taskvar=="red4",c(PCx,PCy)])
conv2<-c(conv2,conv2[1])
polygon(pcaCor$scores[taskvar=="red4",c(PCx,PCy)][conv2,],border=FALSE,col=adjustcolor("red",alpha=0.2),lwd=0.5)
points(pcaCor$scores[,c(PCx,PCy)],pch=19,col=taskvar)
abline(v=0,h=0,lty=1,lwd=0.5,col="lightgray")

# Add loadings (arrows) to the plot to indicate variable contributions
loadfact<-2
for(i in 1:dim(pcaCor$loadings)[1]){
  arrows(0,0,pcaCor$loadings[i,c(PCx)]*loadfact,pcaCor$loadings[i,c(PCy)]*loadfact,length = 0.07)
  text(t(rbind(c(0,0),pcaCor$loadings[i,c(PCx,PCy)]*loadfact)[2,]),labels =rownames(pcaCor$loadings)[i],cex=1.2,pos=1)
}

#Print variance values
colnames(Variance)<-c("eigenvalues","% Variance", "Cumulative %")
Variance
write.table(Variance, "Variance.csv",sep=";",dec=",",col.names=NA,row.names = T) #Save to a CSV file

# Perform statistical tests on the selected principal components (scores)
wilcox.test(pcaCor$scores[,c(PCx)]~as.numeric(as.factor(taskvar)),paired=T)
wilcox.test(pcaCor$scores[,c(PCy)]~as.numeric(as.factor(taskvar)),paired=T)

# Compute Wilcoxon test statistics for all principal components (scores)
pvalues<-NULL
statistics<-NULL
for (i in c(1:dim(pcaCor$scores)[2])) {
  testp<-wilcox.test(pcaCor$scores[,i]~as.numeric(as.factor(taskvar)),paired=TRUE)
  pvalues[i]<-testp$p.value
  statistics[i]<-testp$statistic
}

table<-cbind(statistics,pvalues)
write.table(table, "Wiltest.csv",sep=";",dec=",",col.names=NA,row.names = T) #Save to a CSV file

# Extract PCA loadings values
ab<-print(pcaCor$loadings, digits=8, cutoff=0.01)
write.table(t(ab), "pca_loadings.csv",sep=";",dec=",",col.names=NA,row.names = T) #Save to a CSV file

# Identify the presence of outliers
boxplot(pcaCor$scores[1:23,c(PCx)],pcaCor$scores[24:46,c(PCx)],names=c("H hold","F hold"))
text(rep(1,23),pcaCor$scores[1:23,c(PCx)],combined_ID[1:23])
text(rep(2,23),pcaCor$scores[24:46,c(PCx)],combined_ID[24:46])

boxplot(pcaCor$scores[1:23,c(PCy)],pcaCor$scores[24:46,c(PCy)],names=c("H hold","F hold"))
text(rep(1,23),pcaCor$scores[1:23,c(PCy)],combined_ID[1:23])
text(rep(2,23),pcaCor$scores[24:46,c(PCy)],combined_ID[24:46])





############ ---------> PLS <----------- ##################

# Define the selected EEG and EMG channels for analysis
selEEG_Block1<-c("O2","Oz","CP6","F3","CP5")
selEMG_Block2<-c("FCR","HTE","FCU","FPL","DI1","TE")
selChannels<-c(selEEG_Block1,selEMG_Block2)


dataPLS<-rbind(task1[,selChannels],task2[,selChannels])
meanPLS<-(task1[,selChannels]+task2[,selChannels])/2            
dataPLS<-scale(rbind(task1[,selChannels]-meanPLS,task2[,selChannels]-meanPLS) ) #Center and scale the data
taskvar<-c(rep("orange2",dim(task1)[1]),rep("red4",dim(task2)[1])) #Assign colors to different groups for visualization

library(Morpho)

# Separate EEG and EMG data blocks for PLS analysis
block1<-dataPLS[,which(colnames(task1[,selChannels])%in%selEEG_Block1)]
block2<-dataPLS[,which(colnames(task1[,selChannels])%in%selEMG_Block2)]

# Perform PLS analysis
pls<-pls2B(block1,block2,rounds = 1000)
plot(pls$Xscores[,1],pls$Yscores[,1],col=taskvar,pch=19,xlab="EEG",ylab="EMG") #Plot PLS results

# Covariation and correlation coefficients from PLS results
Stat_pls <-pls$CoVar
write.table(Stat_pls,sep=";",dec=",",file = paste("Covariance","Hhold","Fhold",".csv",sep="_")) #Save to a CSV file

# Compute the PLS effects (loadings)
plsEffects <- plsCoVar(pls,i=1) # select PLS axis
colnames(plsEffects$x)<-selEEG_Block1
colnames(plsEffects$y)<-selEMG_Block2

write.table(plsEffects$x,sep=";",dec=",",file = paste("PLS_EEG_effect","Hhold","Fhold",".csv",sep="_")) #Save to a CSV file
write.table(plsEffects$y,sep=";",dec=",",file = paste("PLS_EMG_effect","Hhold","Fhold",".csv",sep="_")) #Save to a CSV file
