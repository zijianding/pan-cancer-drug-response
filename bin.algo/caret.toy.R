library(caret)
set.seed(1)
dat <- twoClassSim(1000, intercept = -50,linearVars=1000,noiseVars=50,corrVars=10)
datpart = createDataPartition(y=dat$Class, times=1, p=0.5, list=F)
training = dat[datpart[,1],]
testing = dat[-datpart[,1],]

table(training$Class)

#nmin <- sum(training$Class == "Class2")
nmin <- nrow(training)/2
mintraining = training[training$Class=="Class2",]
training = rbind(training,mintraining,mintraining, mintraining)

ctrl <- trainControl(method = "cv",
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)
library(doParallel)
library(foreach)
no_cores = detectCores()
cl = makeCluster(no_cores)
registerDoParallel(cl)
set.seed(2)
library(randomForest);library(pROC)
rfDownsampled <- train(Class ~ ., data = training,
                       method = "rf",
                       ntree = 1500,
                       tuneLength = 5,
                       metric = "ROC",
                      trControl = ctrl,
                       ## Tell randomForest to sample by strata. Here, 
                         ## that means within each class
                         strata = training$Class,
                       ## Now specify that the number of samples selected
                         ## within each class should be the same
                         sampsize = rep(nmin, 2),
                      importance=T)
stopImplicitCluster()
stopCluster(cl)
downProbs <- predict(rfDownsampled, testing, type = "prob")
downsampledROC <- roc(response = testing$Class, 
                        predictor = downProbs,
                        levels = rev(levels(testing$Class)))
plot(downsampledROC, col = rgb(1, 0, 0, .5), lwd = 2)
