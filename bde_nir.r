#install.packages("SoyNAM")
#install.packages("NAM")
#install.packages("rrBLUP")
#install.packages('FSelector')
#install.packages('RWeka')
#install.packages('Biocomb')
#install.packages('arules')
#install.packages('lsr')
#install.packages('e1071')
#install.packages('pROC')
#install.packages('caret')

#Load packages
library(openxlsx) 
library(SoyNAM)
library(NAM)
library(rrBLUP)
library(FSelector)
library(stats)
library(lsr)
library(e1071)
library(caret)
#library(pROC)
#library(arules)
#library(Biocomb)
#library('RWeka')

rm(list=ls())

#Load data
data(soybase)
dim(gen.qa)
dim(data.line.qa)
Markers.src = gen.qa
any(duplicated(rownames(Markers.src)))
any(duplicated(colnames(Markers.src)))
idx_arr <- data.line.qa$strain %in% rownames(Markers.src)
Pheno.src <- data.line.qa[idx_arr,]
any(duplicated(Pheno.src$strain))
idx_vec <- match(Pheno.src$strain, rownames(Markers.src))
Markers.src <- Markers.src[idx_vec,]
dim(Markers.src)

#Data visualisation
hist(Pheno.src$yield, col="gold")
boxplot(Pheno.src$yield ~ Pheno.src$location, xlab = "Location", ylab = "Yield", main = "Yield by Location", col = "pink")
boxplot(Pheno.src$yield ~ Pheno.src$environ, xlab = "Evironment", ylab = "Yield", main = "Yield by Environment", col = "cyan")

#Let's take a part of data to speed up computations
USED_YEAR = c(2012)
podmnogestvo <- Pheno.src$year %in% USED_YEAR
Pheno <- Pheno.src[podmnogestvo,]
Markers <- Markers.src[podmnogestvo,]
USED_LOCATION = c('NE')
podmnogestvo <- Pheno$location %in% USED_LOCATION
Pheno <- Pheno[podmnogestvo,]
Markers <- Markers[podmnogestvo,]
podmnogestvo <- which(is.na(Pheno$yield)) #delete 'NA' yield-data
Pheno <- Pheno[-podmnogestvo,]
Markers <- Markers[-podmnogestvo,]

Markers.2012.NE <- Markers
Pheno.2012.NE <- Pheno

USED_SAMPLES = 500
set.seed(12)
podmnogestvo <- sample(1:nrow(Markers), USED_SAMPLES)
Markers <- Markers[podmnogestvo,]
Pheno <- Pheno[podmnogestvo,]

# Imputation
sum(is.na(Markers))/nrow(Markers)/ncol(Markers)
Markers[Markers == 0] <- -1
Markers[Markers == 1] <- 0
Markers[Markers == 2] <- 1
summary(Markers[,1:10])
impute <- A.mat (Markers, max.missing = 0.4, impute.method = "mean", return.imputed = T, n.core = 12)
Markers_impute <- impute$imputed
dim(Markers_impute)
summary(Markers_impute[,1:10])

USED_MARKERS = 300
podmnogestvo <- sample(1:ncol(Markers_impute), USED_MARKERS)
Markers_impute <- Markers_impute[,podmnogestvo]

# Form training & test samples
dim(Markers_impute)
n <- nrow(Markers_impute)
train_set.idx <- as.matrix(sample(1:n, 0.75 * n))
test_set.idx <- setdiff(1:n, train_set.idx)

Pheno_train <- Pheno[train_set.idx,]
m_train <- Markers_impute[train_set.idx,]
num_of_samples <- nrow(m_train)

############### Count feature scores by 3 groups ##################

#Create a list of Features
feature_names <- colnames(m_train)
Features <- list(feature_names)

#Create feature_score_dimnames list to name feature_score matrices
feature_score_dimnames <- list()
feature_score_dimnames[['score']] <- c('average_value', 'variance', 'standard_deviation', 'number_of_samples', 'Sum', 'sum_of_squares')
feature_score_dimnames[['groups']] <- c('negative', 'null', 'positive')

for (j in 1:length(feature_names)) {
  Features[[feature_names[j]]][['yield']] <- list(feature_score_dimnames[['groups']]) #Form 3 groups of feature_yield data
  Features[[feature_names[j]]][['score']] <- matrix(data = NA, nrow = 6, ncol = 3, dimnames = feature_score_dimnames) #Form feature_score matrix
  #Write feature_yield data by groups
  for (i in 1:nrow(m_train)) {
    if (m_train[i,j] < -0.1) {
      Features[[feature_names[j]]][['yield']][['negative']] <- c(Features[[feature_names[j]]][['yield']][['negative']], Pheno_train$yield[i])
      m_train[i,j] <- -1
    } else if (m_train[i,j] > 0.1) {
      Features[[feature_names[j]]][['yield']][['positive']] <- c(Features[[feature_names[j]]][['yield']][['positive']], Pheno_train$yield[i])
      m_train[i,j] <- 1
    } else {
      Features[[feature_names[j]]][['yield']][['null']] <- c(Features[[feature_names[j]]][['yield']][['null']], Pheno_train$yield[i])
      m_train[i,j] <- 0
    }
  }
  #Count and write scores ('average_value', 'variance', 'standard_deviation', 'number_of_samples', 'Sum', 'sum_of_squares') in feature_score matrix:
  #for 'negative' group
  if (is.null(Features[[feature_names[j]]][['yield']][['negative']]) == F) {
    Features[[feature_names[j]]][['score']][1,1] <- mean(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][2,1] <- var(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][3,1] <- sd(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][4,1] <- length(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][5,1] <- sum(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][6,1] <- sum((Features[[feature_names[j]]][['yield']][['negative']])^2)
  }
  #for 'null' group
  if (is.null(Features[[feature_names[j]]][['yield']][['null']]) == F) {
    Features[[feature_names[j]]][['score']][1,2] <- mean(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][2,2] <- var(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][3,2] <- sd(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][4,2] <- length(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][5,2] <- sum(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][6,2] <- sum((Features[[feature_names[j]]][['yield']][['null']])^2)
  }
  #for 'positive' group
  if (is.null(Features[[feature_names[j]]][['yield']][['positive']]) == F) {
    Features[[feature_names[j]]][['score']][1,3] <- mean(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][2,3] <- var(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][3,3] <- sd(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][4,3] <- length(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][5,3] <- sum(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][6,3] <- sum((Features[[feature_names[j]]][['yield']][['positive']])^2)
  }
  
  #Count average value of yield by feature
  Features[[feature_names[j]]][['average_yield']] <- mean(c(Features[[feature_names[j]]][['score']][1,1], Features[[feature_names[j]]][['score']][1,2],
                                                            Features[[feature_names[j]]][['score']][1,3]))
  
  #Count One_Way_ANOVA scores
  Features[[feature_names[j]]][['One_Way_ANOVA']][['N']] <- sum(Features[[feature_names[j]]][['score']][4,1], Features[[feature_names[j]]][['score']][4,2],
                                                                Features[[feature_names[j]]][['score']][4,3]) #Total number of samples used in feature_score counts
  Features[[feature_names[j]]][['One_Way_ANOVA']][['SSt']] <- sum(c((Features[[feature_names[j]]][['yield']][['negative']] - Features[[feature_names[j]]][['average_yield']])^2,
                                                                    (Features[[feature_names[j]]][['yield']][['null']] - Features[[feature_names[j]]][['average_yield']])^2,
                                                                    (Features[[feature_names[j]]][['yield']][['positive']] - Features[[feature_names[j]]][['average_yield']])^2)) #Total sum of squares (SSt = SSbw + SSwg)
  Features[[feature_names[j]]][['One_Way_ANOVA']][['SSwg']] <- sum(c((Features[[feature_names[j]]][['yield']][['negative']] - Features[[feature_names[j]]][['score']][1,1])^2,
                                                                     (Features[[feature_names[j]]][['yield']][['null']] - Features[[feature_names[j]]][['score']][1,2])^2,
                                                                     (Features[[feature_names[j]]][['yield']][['positive']] - Features[[feature_names[j]]][['score']][1,3])^2)) #Sum of squares within groups
  Features[[feature_names[j]]][['One_Way_ANOVA']][['SSbg']] <- sum(c(Features[[feature_names[j]]][['score']][4,1]*(Features[[feature_names[j]]][['score']][1,1] - Features[[feature_names[j]]][['average_yield']])^2,
                                                                     Features[[feature_names[j]]][['score']][4,2]*(Features[[feature_names[j]]][['score']][1,2] - Features[[feature_names[j]]][['average_yield']])^2,
                                                                     Features[[feature_names[j]]][['score']][4,3]*(Features[[feature_names[j]]][['score']][1,3] - Features[[feature_names[j]]][['average_yield']])^2)) #Sum of squares between groups
  Features[[feature_names[j]]][['One_Way_ANOVA']][['MSwg']] <- Features[[feature_names[j]]][['One_Way_ANOVA']][['SSwg']] / (Features[[feature_names[j]]][['One_Way_ANOVA']][['N']] - 3) #Mean of squares within groups
  Features[[feature_names[j]]][['One_Way_ANOVA']][['MSbg']] <- Features[[feature_names[j]]][['One_Way_ANOVA']][['SSbg']] / 2 #Mean of squares between groups
  Features[[feature_names[j]]][['One_Way_ANOVA']][['F_value']] <- Features[[feature_names[j]]][['One_Way_ANOVA']][['MSbg']] / Features[[feature_names[j]]][['One_Way_ANOVA']][['MSwg']] #Result score of the test
  
  #Count Fisher_Score
  Features[[feature_names[j]]][['Fisher_Score']] <- sum(c(Features[[feature_names[j]]][['score']][4,1]*(Features[[feature_names[j]]][['score']][1,1] - Features[[feature_names[j]]][['average_yield']])^2,
                                                          Features[[feature_names[j]]][['score']][4,2]*(Features[[feature_names[j]]][['score']][1,2] - Features[[feature_names[j]]][['average_yield']])^2,
                                                          Features[[feature_names[j]]][['score']][4,3]*(Features[[feature_names[j]]][['score']][1,3] - Features[[feature_names[j]]][['average_yield']])^2)
  ) / sum(c(Features[[feature_names[j]]][['score']][4,1]*Features[[feature_names[j]]][['score']][2,1]^2,
            Features[[feature_names[j]]][['score']][4,2]*Features[[feature_names[j]]][['score']][2,2]^2,
            Features[[feature_names[j]]][['score']][4,3]*Features[[feature_names[j]]][['score']][2,3]^2))
}

#Write scores in vectors
One_Way_ANOVA <- vector(length = length(feature_names))
Fisher_Score <- vector(length = length(feature_names))
for (j in 1:length(feature_names)) {
  One_Way_ANOVA[j] <- Features[[feature_names[j]]][['One_Way_ANOVA']][['F_value']]
  Fisher_Score[j] <- Features[[feature_names[j]]][['Fisher_Score']]
}
names(Fisher_Score) <- feature_names
names(One_Way_ANOVA) <- feature_names

#Compare One_Way_ANOVA & Fisher_Score results
One_Way_ANOVA.sorted <- sort(One_Way_ANOVA, decreasing = T)
Fisher_Score.sorted <- sort(Fisher_Score, decreasing = T)
which(rank(-One_Way_ANOVA.sorted) <= 10)
which(rank(-Fisher_Score.sorted) <= 10)

#Split 'yield' data by 16 classes & form a dataframe to count IG
yield.max <- max(Pheno_train$yield)
yield.min <- min(Pheno_train$yield) - 1
yield.seq_16 <- (yield.max - yield.min) / 16
hist(Pheno_train$yield, breaks = seq(yield.min, yield.max, yield.seq_16), col = 'orange')
yield_classes_16 <- cut(Pheno_train$yield, breaks = seq(yield.min, yield.max, yield.seq_16), labels = letters[1:16])
any(is.na(yield_classes_16))
m_train_classes_16 <- cbind(m_train, as.vector(yield_classes_16))
colnames(m_train_classes_16) <- c(colnames(m_train), 'classes')
m_train_classes_16 <- as.data.frame(m_train_classes_16)
is.data.frame(m_train_classes_16)

#Count IG & check results
IG.FSelector_16 <- information.gain(classes ~ ., data=m_train_classes_16)
IG.FSelector.vector_16 <- as.vector(IG.FSelector_16[,1])
names(IG.FSelector.vector_16) <- rownames(IG.FSelector_16)
IG.FSelector.sorted_16 <- sort(IG.FSelector.vector_16, decreasing = T)
which(rank(-IG.FSelector.sorted_16) <= 10)

#Normalization of the scores
normalize <- function(x) {(x - min(x)) / (max(x) - min(x))} #Create a function

One_Way_ANOVA.norm <- normalize(One_Way_ANOVA)
Fisher_Score.norm <- normalize(Fisher_Score)
IG.FSelector_16.norm <- normalize(IG.FSelector.vector_16)

#Leave 100 best from each filter method
One_Way_ANOVA.best_100 <- One_Way_ANOVA.sorted[1:100]
Fisher_Score.best_100 <- Fisher_Score.sorted[1:100]
IG.FSelector_16.best_100 <- IG.FSelector.sorted_16[1:100]

#Count overlapping percentage between three filter methods
One_Way_ANOVA.best_100_VS_Fisher_Score.best_100 <- length(which(names(One_Way_ANOVA.best_100) %in% names(Fisher_Score.best_100)))
One_Way_ANOVA.best_100_VS_IG.FSelector_16.best_100 <- length(which(names(One_Way_ANOVA.best_100) %in% names(IG.FSelector_16.best_100)))
Fisher_Score.best_100_VS_IG.FSelector_16.best_100 <- length(which(names(Fisher_Score.best_100) %in% names(IG.FSelector_16.best_100)))

#Form feature pool & count total score
feature_pool.idx <- unique(c(match(names(One_Way_ANOVA.best_100), names(One_Way_ANOVA)), match(names(Fisher_Score.best_100), names(Fisher_Score)),
                             match(names(IG.FSelector_16.best_100), names(IG.FSelector.vector_16))))

Total_Score <- One_Way_ANOVA.norm[feature_pool.idx] + Fisher_Score.norm[feature_pool.idx] + IG.FSelector_16.norm[feature_pool.idx]

#Count probability of features
p_of_feature <- (max(Total_Score) - Total_Score) / (max(Total_Score) - min(Total_Score))
D <- length(p_of_feature) #number of features in feature pool
for (i in 1:D) {
  p_of_feature[i] <- 1 - min(0.9, max(0.1, p_of_feature[i])) 
}

#Count Correlation based Feature Selection (CFS) in feature pool
podmnogestvo <- match(names(p_of_feature), colnames(m_train))
m_train_tot.sc <- as.data.frame(m_train[,podmnogestvo])
dim(m_train_tot.sc)

k <- length(p_of_feature) #cardinality of the feature subset (k = D)
CFS_score <- vector(length = k)
names(CFS_score) <- names(p_of_feature)
class_corr <- vector(length = k)
names(class_corr) <- names(p_of_feature)
inter_corr <- matrix(nrow = k, ncol = k)
colnames(inter_corr) <- names(p_of_feature)
rownames(inter_corr) <- names(p_of_feature)
inter_corr.av <- vector(length = k) #mean value of each feature from inter_corr
names(inter_corr.av) <- names(p_of_feature)

for (i in 1:k) {
  tbl <- table(m_train_tot.sc[,i], Pheno_train$yield)
  class_corr[i] <- cramersV(tbl)
}
class_corr.av <- mean(class_corr)
for (i in 1:k) {
  for (j in 1:k) {
    tbl <- table(m_train_tot.sc[,i], m_train_tot.sc[,j])
    inter_corr[i,j] <- cramersV(tbl)
  }
  inter_corr.av[i] <- mean(inter_corr[i,])
  CFS_score[i] <- k*class_corr.av / sqrt(k + k*(k-1)*inter_corr.av[i])
}
 
CFS_score.df <- as.data.frame(CFS_score)
CFS_score.best_names_50 <- cutoff.k(CFS_score.df, 50)

####### BDE ######
NP <- 20
Maxiter <- 15
mutation_factor <- 0.3
CR <- 0.5
#set.seed(12)

#Create initial population - Population[[1]]
Population <- list()
Population[['X']][[1]] <- matrix(nrow = NP, ncol = D)
colnames(Population[['X']][[1]]) <- names(p_of_feature)
x_init_CFS.idx_best_50 <- match(CFS_score.best_names_50, names(p_of_feature))
Population[['X']][[1]][1,] <- rep.int(1, D)
Population[['X']][[1]][1, x_init_CFS.idx_best_50] <- rep.int(0, length(x_init_CFS.idx_best_50)) #first individual (x_1) is best_50 from CFS

for (i in 2:NP) {
  for (j in 1:k) {
    if (p_of_feature[j] < runif(1)) {
      Population[['X']][[1]][i,j] <- 1
    } else {
      Population[['X']][[1]][i,j] <- 0
    }
  }
}

### Fitness evaluation with SVM ###
Population[['Fitness']][[1]] <- vector(length = NP) #procentage of coincidences with default parameters of svm

for (i in 1:NP) {
  #Features acting in fitness evaluation
  svm_feature.idx <- which(Population[['X']][[1]][i,] == 0)
  m_train_individual <- m_train_tot.sc[,svm_feature.idx]
  #Creating dummies from (-1,0,1) encoding for svm calculations
  m_train_svm <- matrix(data = 0, nrow = num_of_samples, ncol = 3*length(svm_feature.idx))
  m_train_svm.names <- vector()
  for (j in 1:length(svm_feature.idx)) {
    m_train_svm.names <- c(m_train_svm.names, paste(names(m_train_individual)[j],"homo(-1)",sep = "_"), paste(names(m_train_individual)[j],"hetero(0)",sep = "_"),
                           paste(names(m_train_individual)[j],"homo(1)",sep = "_"))
  }
  colnames(m_train_svm) <- m_train_svm.names
  rownames(m_train_svm) <- rownames(m_train)
  
  for (l in 1:length(svm_feature.idx)) {
    for (j in 1:num_of_samples) {
      if (m_train_individual[j,l] == -1) {
        m_train_svm[j,(3*l-2)] <- 1
      } else if (m_train_individual[j,l] == 1) {
        m_train_svm[j,3*l] <- 1
      } else {
        m_train_svm[j,(3*l-1)] <- 1
      }
    }
  }
  
  svm_model <- svm(x=m_train_svm,y=Pheno_train$yield,type='eps-regression',cross = 10, cost = 1, gamma = 0.01)
  
  #Write Fitness
  Population[['Fitness']][[1]][i] <- sqrt(svm_model$tot.MSE)
}

#Find the best individual from the population
Population$x_best[[1]] <- which.min(Population[['Fitness']][[1]])

#Count hamming based population diversity for initial population
gtype <- vector()
l <- 2
for (i in 1:(NP-1)) {
  for (j in l:NP) {
    ham <- sum(abs(Population[['X']][[1]][i,] - Population[['X']][[1]][j,]))
    gtype <- sum(gtype, ham)
  }
  l <- l + 1
}
Population[['gtype']][[1]] <- (gtype/D) / ((NP*(NP-1))/2)


### Go on with BDE ###

for (G in 1:Maxiter) {
  Population[['X']][[G+1]] <- matrix(nrow = NP, ncol = D)
  colnames(Population[['X']][[G+1]]) <- names(p_of_feature)
  Population[['Fitness']][[G+1]] <- vector(length = NP)
  
  #Count C_min
  C_min <- ceiling(13*(1 - G/Maxiter)) + 4
  
  for (i in 1:NP) {
    target <- Population[['X']][[G]][i,]
    #Generation of the mutation operator
    r1 <- sample(c(1:NP)[-i], 1)
    r2 <- sample(c(1:NP)[-c(i,r1)], 1)
    dif <- vector(length = D)
    for (j in 1:D) {
      if (Population[['X']][[G]][r1,j] == Population[['X']][[G]][r2,j]) {
        dif[j] <- 0
      } else {
        dif[j] <- Population[['X']][[1]][r1,j]
      }
    }
    #Generation of the mutant vector
    mutant <- vector(length = D)
    for (j in 1:D) {
      if (dif[j] == 1 && runif(1) < mutation_factor) {
        mutant[j] <- 1
      } else {
        mutant[j] <- Population[['X']][[G]][Population$x_best[[G]],j]
      }
    }
    #Generation of the trial vector
    trial <- vector(length = D)
    for (j in 1:D) {
      if (runif(1) < CR) {
        trial[j] <- mutant[j]
      } else {
        trial[j] <- target[j]
      }
    }
    #Modification of the trial vector with C_min
    modif <- sample(1:D, C_min)
    trial[modif[j]] <- abs(trial[modif[j]] - 1)
    
    ### Fitness evaluation with SVM ###
    
    #Features acting in fitness evaluation
    svm_feature.idx <- which(Population[['X']][[G]][i,] == 0)
    m_train_individual <- m_train_tot.sc[,svm_feature.idx]
    #Creating dummies from (-1,0,1) encoding for svm calculations
    m_train_svm <- matrix(data = 0, nrow = num_of_samples, ncol = 3*length(svm_feature.idx))
    m_train_svm.names <- vector()
    for (j in 1:length(svm_feature.idx)) {
      m_train_svm.names <- c(m_train_svm.names, paste(names(m_train_individual)[j],"homo(-1)",sep = "_"), paste(names(m_train_individual)[j],"hetero(0)",sep = "_"),
                             paste(names(m_train_individual)[j],"homo(1)",sep = "_"))
    }
    colnames(m_train_svm) <- m_train_svm.names
    rownames(m_train_svm) <- rownames(m_train)
    
    for (l in 1:length(svm_feature.idx)) {
      for (j in 1:num_of_samples) {
        if (m_train_individual[j,l] == -1) {
          m_train_svm[j,(3*l-2)] <- 1
        } else if (m_train_individual[j,l] == 1) {
          m_train_svm[j,3*l] <- 1
        } else {
          m_train_svm[j,(3*l-1)] <- 1
        }
      }
    }
    
    svm_model <- svm(x=m_train_svm,y=Pheno_train$yield,type='eps-regression',cross = 10, cost = 1, gamma = 0.01)
    
    #Write Fitness
    trial.fit <- sqrt(svm_model$tot.MSE)
    
    # Write best (target or trial) in the next generation
    if (trial.fit < Population[['Fitness']][[G]][i]) {
      Population[['X']][[G+1]][i,] <- trial
      Population[['Fitness']][[G+1]][i] <- trial.fit
    } else {
      Population[['X']][[G+1]][i,] <- target
      Population[['Fitness']][[G+1]][i] <- Population[['Fitness']][[G]][i]
    }
  }
  
  #Find the best individual from the population
  Population[['x_best']][[G+1]] <- which.min(Population[['Fitness']][[G+1]])
  
  #Count hamming based population diversity & CR
  gtype <- vector()
  l <- 2
  for (i in 1:(NP-1)) {
    for (j in l:NP) {
      ham <- sum(abs(Population[['X']][[G+1]][i,] - Population[['X']][[G+1]][j,]))
      gtype <- sum(gtype, ham)
    }
    l <- l + 1
  }
  Population[['gtype']][[G+1]] <- (gtype/D) / ((NP*(NP-1))/2)
  CR <- Population[['gtype']][[G+1]] / Population[['gtype']][[1]]
  CR <- 1 - min(0.9, max(0.1, CR))
}

#Best Fitness
Population$Fitness[[G+1]][Population$x_best[[G+1]]]
#Final features
final_features.names <- colnames(Population$X[[Maxiter+1]][,which(Population$X[[Maxiter+1]][Population$x_best[[Maxiter+1]],] == 0)])
final_features.names
length(final_features.names)
#Best features in first generation
best_features_in_G1.names <- colnames(Population$X[[1]][,which(Population$X[[1]][Population$x_best[[1]],] == 0)])
best_features_in_G1.names
length(best_features_in_G1.names)

#Number of coincidences of feautures in first and last Generartion
length(which(final_features.names %in% best_features_in_G1.names))

write.csv(final_features.names, file = "final_features_BDE.csv")
write.csv(best_features_in_G1.names, file = "best_features_in_G1_BDE.csv")

### Predicting test data ###
m_test <- Markers_impute[test_set.idx,]
Pheno_test <- Pheno[test_set.idx,]

#Train svm
m_train_individual <- m_train_tot.sc[,final_features.names]
m_train_svm <- matrix(data = 0, nrow = num_of_samples, ncol = 3*length(final_features.names))
m_train_svm.names <- vector()
for (j in 1:length(final_features.names)) {
  m_train_svm.names <- c(m_train_svm.names, paste(names(m_train_individual)[j],"homo(-1)",sep = "_"), paste(names(m_train_individual)[j],"hetero(0)",sep = "_"),
                         paste(names(m_train_individual)[j],"homo(1)",sep = "_"))
}
colnames(m_train_svm) <- m_train_svm.names
rownames(m_train_svm) <- rownames(m_train)

for (l in 1:length(final_features.names)) {
  for (j in 1:num_of_samples) {
    if (m_train_individual[j,l] == -1) {
      m_train_svm[j,(3*l-2)] <- 1
    } else if (m_train_individual[j,l] == 1) {
      m_train_svm[j,3*l] <- 1
    } else {
      m_train_svm[j,(3*l-1)] <- 1
    }
  }
}
svm_model <- svm(x=m_train_svm,y=Pheno_train$yield,type='eps-regression', cost = 1, gamma = 0.01)

#Predict 
m_test_svm <- matrix(data = 0, nrow = nrow(m_test), ncol = 3*length(final_features.names))
m_test_individual <- m_test[,final_features.names]
colnames(m_test_svm) <- m_train_svm.names
rownames(m_test_svm) <- rownames(m_test)
for (l in 1:length(final_features.names)) {
  for (j in 1:nrow(m_test_svm)) {
    if (m_test_individual[j,l] == -1) {
      m_test_svm[j,(3*l-2)] <- 1
    } else if (m_test_individual[j,l] == 1) {
      m_test_svm[j,3*l] <- 1
    } else {
      m_test_svm[j,(3*l-1)] <- 1
    }
  }
}
prediction <- predict(svm_model, newdata = m_test_svm)
sqrt(mean((Pheno_test$yield - prediction)^2))

#Plot best Fitness by Generations
Fitness_best <- vector(length = Maxiter+1)
for (i in 1:(Maxiter+1)) {
  Fitness_best[i] <- Population$Fitness[[i]][Population$x_best[[i]]]
}
plot(Fitness_best, type = 'b', col = 'blue', xlab = 'Number of Population', ylab = 'RMSE',
     main = 'Improvement of prediction \n from population to population')

