rm(list=ls())

# Загрузка пакетов
library(SoyNAM)
# Загрузка данных
data(soybase)
dim(gen.qa)
dim(data.line.qa)
ls()
Markers.src = gen.qa
idx_arr <- data.line.qa$strain %in% rownames(Markers.src)
Pheno.src <- data.line.qa[idx_arr,]
idx_vec <- match(Pheno.src$strain, rownames(Markers.src))
Markers.src <- Markers.src[idx_vec,]


# Визуализация данных
hist(Pheno.src$yield, col="gold")
boxplot(Pheno.src$yield ~ Pheno.src$location, xlab = "Location", ylab = "Yield", main = "Yield
by Location", col = "pink")
boxplot(Pheno.src$yield ~ Pheno.src$environ, xlab = "Evironment", ylab = "Yield", main = "Yield by Environment", col = "cyan")

# Возьмем подмножество данных для ускорения вычислений
USED_SAMPLES = 10000
set.seed(12)
podmnogestvo <- sample(1:nrow(Markers.src), USED_SAMPLES)
Markers <- Markers.src[podmnogestvo,]
Pheno <- Pheno.src[podmnogestvo,]
USED_MARKERS = 300
podmnogestvo <- sample(1:ncol(Markers.src), USED_MARKERS)
Markers <- Markers[,podmnogestvo]
Markers.sv <- Markers
Pheno.sv <- Pheno
USED_YEAR = c(2012)
podmnogestvo <- Pheno$year %in% USED_YEAR
Pheno <- Pheno[podmnogestvo,]
Markers <- Markers[podmnogestvo,]

# Heritability Наследуемость
library(NAM)
FIT = reml ( y = Pheno$yield, X = ~ as.factor(Pheno$location), Z = ~ as.factor(Pheno$fam) )
FIT$VC
#Общая дисперсия (VT) может быть представлена в виде суммы дисперсии,
#связанной с различиями в генотипе (VG), и дисперсии, связанной с влиянием среды (VE). 
#Наследуемость в широком смысле понимается как коэффициент генетической детерминации (H??)

# Загрузка пакета
library(rrBLUP)
USED_LOCATION = c('NE')
podmnogestvo <- Pheno$location %in% USED_LOCATION
Pheno <- Pheno[podmnogestvo,]
Markers <- Markers[podmnogestvo,]

# Импутация
sum(is.na(Markers))/nrow(Markers)/ncol(Markers)
Markers[Markers == 2] <- -1
summary(Markers[,1:10])
impute <- A.mat (Markers, max.missing = 0.5, impute.method = "mean", return.imputed = T)
Markers_impute <- impute$imputed
summary(Markers_impute[,1:10])
Markers_impute2 <- Markers_impute

# Выборки
n <- nrow(Markers_impute2)
p <- ncol(Markers_impute2)
n
p
train <- as.matrix(sample(1:n, 0.75 * n))
test <- setdiff(1:n, train)


#При каждом запуске, т.к. используется генератор псевдослучайных чисел,
#будут разные разбиения и разные результаты. Для работающей модели
#результаты должны отличаться не сильно.

Pheno_train <- Pheno[train,]
m_train <- Markers_impute2[train,] # feature для каждого столбца считаем score


write.csv(Pheno_train, file = "C://Users//minee//Documents//MyData.csv")



############### Count feature scores by classes ##################
library('stats')

#Form list of features
feature_names <- colnames(m_train)
Features <- list(feature_names)

#Form dimnames for feature_score matrices
feature_score_dimnames <- list()
feature_score_dimnames[['score']] <- c('average_value', 'dispersion', 'standard_deviation', 'number_of_samples')
feature_score_dimnames[['classes']] <- c('negative', 'null', 'positive')

for (j in 1:ncol(m_train)) {
  Features[[feature_names[j]]][['yield']] <- list(c('negative', 'null', 'positive')) #Form 3 class_vectors of feature_yield data
  Features[[feature_names[j]]][['score']] <- matrix(data = NA, nrow = 4, ncol = 3, dimnames = feature_score_dimnames) #Form feature_score matrix
  #Write feature_yield data in class_vectors
  for (i in 1:nrow(m_train)) {
    if (is.na(Pheno_train$yield[i]) == F) {
      if (m_train[i,j] < 0) {
        Features[[feature_names[j]]][['yield']][['negative']] <- c(Features[[feature_names[j]]][['yield']][['negative']], Pheno_train$yield[i])
      } else if (m_train[i,j] > 0) {
        Features[[feature_names[j]]][['yield']][['positive']] <- c(Features[[feature_names[j]]][['yield']][['positive']], Pheno_train$yield[i])
      } else {
        Features[[feature_names[j]]][['yield']][['null']] <- c(Features[[feature_names[j]]][['yield']][['null']], Pheno_train$yield[i])
      }
    }
  }
  #Count and write scores ('average_value', 'dispersion', 'standard_deviation', 'number_of_samples') in feature_score matrix:
  #for 'negative' feature_class
  if (is.null(Features[[feature_names[j]]][['yield']][['negative']]) == F) {
    Features[[feature_names[j]]][['score']][1,1] <- mean(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][2,1] <- var(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][3,1] <- sd(Features[[feature_names[j]]][['yield']][['negative']])
    Features[[feature_names[j]]][['score']][4,1] <- length(Features[[feature_names[j]]][['yield']][['negative']])
  }
  #for 'null' feature_class
  if (is.null(Features[[feature_names[j]]][['yield']][['null']]) == F) {
    Features[[feature_names[j]]][['score']][1,2] <- mean(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][2,2] <- var(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][3,2] <- sd(Features[[feature_names[j]]][['yield']][['null']])
    Features[[feature_names[j]]][['score']][4,2] <- length(Features[[feature_names[j]]][['yield']][['null']])
  }
  #for 'positive' feature_class
  if (is.null(Features[[feature_names[j]]][['yield']][['positive']]) == F) {
    Features[[feature_names[j]]][['score']][1,3] <- mean(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][2,3] <- var(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][3,3] <- sd(Features[[feature_names[j]]][['yield']][['positive']])
    Features[[feature_names[j]]][['score']][4,3] <- length(Features[[feature_names[j]]][['yield']][['positive']])
  }
  #Count T-Statistics
  Features[[feature_names[j]]][['T-Statistics']] <- (Features[[feature_names[j]]][['score']][1,3] - Features[[feature_names[j]]][['score']][1,2] - Features[[feature_names[j]]][['score']][1,1]
  ) / sqrt((Features[[feature_names[j]]][['score']][3,3])^2/Features[[feature_names[j]]][['score']][4,3]
           + (Features[[feature_names[j]]][['score']][3,2])^2/Features[[feature_names[j]]][['score']][4,2]
           + (Features[[feature_names[j]]][['score']][3,1])^2/Features[[feature_names[j]]][['score']][4,1])
  #Count Average_Value by feature
  Features[[feature_names[j]]][['Average_Value']] <- mean(c(Features[[feature_names[j]]][['yield']][['negative']],Features[[feature_names[j]]][['yield']][['null']],
                                                            Features[[feature_names[j]]][['yield']][['positive']]))
  #Count Fisher_Score
  Features[[feature_names[j]]][['Fisher_Score']] <- (Features[[feature_names[j]]][['score']][4,1]*(Features[[feature_names[j]]][['score']][1,1] - Features[[feature_names[j]]][['Average_Value']])^2
                                                     + Features[[feature_names[j]]][['score']][4,2]*(Features[[feature_names[j]]][['score']][1,2] - Features[[feature_names[j]]][['Average_Value']])^2
                                                     + Features[[feature_names[j]]][['score']][4,3]*(Features[[feature_names[j]]][['score']][1,3] - Features[[feature_names[j]]][['Average_Value']])^2
  ) / (Features[[feature_names[j]]][['score']][4,1]*Features[[feature_names[j]]][['score']][2,1]
       + Features[[feature_names[j]]][['score']][4,2]*Features[[feature_names[j]]][['score']][2,2]
       + Features[[feature_names[j]]][['score']][4,3]*Features[[feature_names[j]]][['score']][2,3])
}








Pheno_test <- Pheno[test,]
m_test <- Markers_impute2[test,]
yield_train <- Pheno_train[,14]

#Вызов функции
yield_model <- mixed.solve(yield_train, Z = m_train, K = NULL, SE = FALSE, return.Hinv = FALSE)
g <- yield_model$u
head(g)
mu <- yield_model$beta[1]

#Оценка точности
yield_predicted <- mu + m_test %*% g
yield_test <- Pheno_test[,14]
yield_accuracy <- cor.test(yield_predicted, yield_test)
yield_accuracy

