rm(list=ls())
gc()


###load libraries
library(randomForest)
library(e1071)
library(ggthemes)
library(grid)
library(ggplot2)
library(tidyr)
library(viridis)

### set working directory
setwd("/data/cecchgu/LHG/plotdbfz36/")

#####load dataset
DATA_BIOMASS <- read.csv('DataBiomass.csv')

##### reproducibility number
set.seed(2511) 
DATA_BIOMASS <- DATA_BIOMASS[,3:29]

###################################################################################################

####     BIOMASS model with training/calibration (85% training, 15% validation) START

###################################################################################################
###########train

# randomly pick 85% of the number of observations
index <- sample(1:nrow(DATA_BIOMASS),size = 0.85*nrow(DATA_BIOMASS))

# subset weather to include only the elements in the index
train <- DATA_BIOMASS[index,]
# train_cat <- DATA_BIOMASS_cat [index,]
# subset weather to include all but the elements in the index
test <- DATA_BIOMASS [-index,]

#print number of data for train and validation
nrow(train)
nrow(test)


##################
# 1 Exponenetial Model
##################

# Create a multiple (log)linear regression model using the training data
log.reg_H_with <- glm(log(BIOMASS+1) ~., data = train)
summary(log.reg_H_with )




# Apply the model to the testing data (i.e., make predictions) ...
# (Don't forget to exponentiate the results to revert the log transformation)
test.pred.lin <- exp(predict(log.reg_H_with,train))-1

# ...and evaluate the accuracy
RMSE.log.reg_H_with <- sqrt(mean((test.pred.lin-train$BIOMASS)^2))
RMSE.log.reg_H_with


MAE.log.reg_H_with <- mean(abs(test.pred.lin-train$BIOMASS))
MAE.log.reg_H_with

Rel.RMSEexp <- sqrt(mean((test.pred.lin-train$BIOMASS)^2)) / diff(range(train$BIOMASS))
Rel.RMSEexp


R2exp <- 1 - sum((train$BIOMASS-test.pred.lin)^2)/sum((train$BIOMASS-mean(train$BIOMASS))^2)
R2exp
##################
# 2 Linear Model
##################

# Create a multiple linear regression model using the training data
lin.reg_H_with <- glm(BIOMASS ~ ., data = train)

summary(lin.reg_H_with )

# Apply the model to the testing data (i.e., make predictions) ...
test.pred.lin_ <- (predict(lin.reg_H_with,train))-1

# ...and evaluate the accuracy
RMSE.lin.reg_H_with <- sqrt(mean((test.pred.lin_-train$BIOMASS)^2))
RMSE.lin.reg_H_with


MAE.lin.reg_H_with <- mean(abs(test.pred.lin_-train$BIOMASS))
MAE.lin.reg_H_with


Rel.RMSElin <- sqrt(mean((test.pred.lin_-train$BIOMASS)^2)) / diff(range(train$BIOMASS))
Rel.RMSElin

R2lin <- 1 - sum((train$BIOMASS-test.pred.lin_)^2)/sum((train$BIOMASS-mean(train$BIOMASS))^2)

R2lin
# R2RF <- 1 - sum((test$BIOMASS-test.pred.forest)^2)/sum((test$BIOMASS-mean(test$BIOMASS))^2)

##################
# 3 RANDOM FOREST
##################

rf_H_with <- randomForest(BIOMASS ~., data = train,importance = TRUE, proximity = T,ntree=500)
print(rf_H_with)
varImpPlot(rf_H_with)
# How many trees are needed to reach the minimum error estimate? 
# This is a simple problem; it appears that about 100 trees would be enough. 
which.min(rf_H_with$mse)

# Plot rf to see the estimated error as a function of the number of trees
# (not running it here)
# plot(rf) 

# Using the importance()  function to calculate the importance of each variable
imp_H_with <- as.data.frame(sort(importance(rf_H_with)[,1],decreasing = TRUE),optional = T,col.names=T)
names(imp_H_with) <- "% Inc MSE"
imp_H_with



test2 <- rbind(train, train)

# As usual, predict and evaluate on the test set
test.pred.forest <- predict(rf_H_with,test2)

test.pred.forest <- test.pred.forest[1:nrow(train)]

RMSE.forest_H_with <- sqrt(mean((test.pred.forest-train$BIOMASS)^2))
RMSE.forest_H_with


MAE.forest_H_with <- mean(abs(test.pred.forest-train$BIOMASS))
MAE.forest_H_with


# Rel.RMSE <- sqrt(mean((test.pred.forest-test$BIOMASS)^2)) / sd(test$BIOMASS)
Rel.RMSErf <- sqrt(mean((test.pred.forest-train$BIOMASS)^2)) / diff(range(train$BIOMASS))
Rel.RMSErf
### ref Remote Sens. 2016, 8, 583; doi:10.3390/rs8070583 
R2RF <- 1 - sum((train$BIOMASS-test.pred.forest)^2)/sum((train$BIOMASS-mean(train$BIOMASS))^2)
R2RF

##################
# 4 SVM
##################


## create rmse function
rmse <- function(error)
{
  sqrt(mean(error^2))
}

##load a model SVM
SVM_H <- svm(BIOMASS ~., data = train)
predictedY <- predict(SVM_H, train)
predictedY <- predictedY[1:nrow(train)]
plot(train$BIOMASS, predictedY, col = "red", pch=4)

error <- train$BIOMASS - predictedY
svrPredictionRMSE <- rmse(error)
svrPredictionRMSE

tuneResult <- tune(svm, BIOMASS ~., data = train, kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

print(tuneResult)
# Draw the tuning graph
plot(tuneResult)

##tune the SVM model
tunedModel <- tuneResult$best.model
test.pred.svm <- predict(tunedModel, train) 
test.pred.svm <- predictedY[1:nrow(train)]
error <- train$BIOMASS - test.pred.svm  


##compute the performance metrics
tunedModelRMSE <- rmse(error)
tunedModelRMSE


RMSE.svm_H_with <- sqrt(mean((test.pred.svm-train$BIOMASS)^2))
RMSE.svm_H_with


MAE.svm_H_with <- mean(abs(test.pred.svm-train$BIOMASS))
MAE.svm_H_with

Rel.RMSEsvm <- sqrt(mean((test.pred.svm-train$BIOMASS)^2)) / diff(range(train$BIOMASS))
Rel.RMSEsvm
# SVM
R2svm <- 1 - sum((train$BIOMASS-test.pred.svm)^2)/sum((train$BIOMASS-mean(train$BIOMASS))^2)
R2svm



# Create a data frame with the error metrics for each method
accuracy_H_with <- data.frame(Method = c("Exp Regression","Linear Regression","Random forest","SVM"),
                              RMSE   = c(RMSE.log.reg_H_with,RMSE.lin.reg_H_with,RMSE.forest_H_with,RMSE.svm_H_with),
                              MAE    = c(MAE.log.reg_H_with,MAE.lin.reg_H_with,MAE.forest_H_with,MAE.svm_H_with),
                              relRMSE    = c(Rel.RMSEexp,Rel.RMSElin,Rel.RMSErf,Rel.RMSEsvm),
                              r2    = c(R2exp,R2lin,R2RF,R2svm)) 
# Round the values and print the table
accuracy_H_with$RMSE <- round(accuracy_H_with$RMSE,2)
accuracy_H_with$MAE <- round(accuracy_H_with$MAE,2) 
accuracy_H_with$r2 <- round(accuracy_H_with$r2,2) 

#print results
accuracy_H_with



# Create a data frame with the predictions for each method
all.predictions <- data.frame(actual = train$BIOMASS,
                              exponential.regression = test.pred.lin,
                              linear.regression_ = test.pred.lin_,
                              random.forest = test.pred.forest,
                              support.vector.machine = test.pred.svm)

index_0df <- which(all.predictions$actual !=0)
all.predictions <- all.predictions[index_0df,]



# Gather the prediction variables (columns) into a single row (i.e., wide to long)
# Recall the ggplot2 prefers the long data format
all.predictions <- gather(all.predictions,key = model,value = predictions,2:5)
all.predictions[,4] <- rep(test[,3],4)

###plot for publication
theme_Publication <- function(base_size=14, base_family="helvetica") {
  
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(face = "bold",size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -0.2),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           # legend.key = element_rect(colour = NA),
           # legend.position = "bottom",
           # legend.direction = "horizontal",
           # legend.key.size= unit(0.2, "cm"),
           # legend.margin = unit(0, "cm"),
           # legend.title = element_text(face="italic"),
           plot.margin=unit(c(10,5,5,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(face="bold")
   ))
  
}




##### create a dataframe with the summary of the stats for the four models
all.predictions$model[all.predictions$model== "linear.regression_"] = "Linear Regression"
all.predictions$model[all.predictions$model== "exponential.regression"] = "Exponential Regression"
all.predictions$model[all.predictions$model== "random.forest"] = "Random Forest"
all.predictions$model[all.predictions$model== "support.vector.machine"] = "Support Vector Machine"

all.predictions$RMSE <- rep(accuracy_H_with$RMSE, each=435)
all.predictions$MAE <- rep(accuracy_H_with$MAE, each=435)
all.predictions$relRMSE <- rep(round(accuracy_H_with$relRMSE, digit=3), each=435)
all.predictions$r2 <- rep(round(accuracy_H_with$r2, digit=2), each=435)



####first plot(base theme)

Map_1 <- ggplot(data = all.predictions,aes(x = actual, y = predictions)) + ###,color= V4
  geom_point(colour = "blue") + #colour = "blue"
  geom_abline(intercept = 0, slope = 1, colour = "black") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
  facet_wrap(~ model,ncol = 2)  +
  xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
  ggtitle("Biomass, Predicted vs. Actual, by model")+ coord_fixed(ratio = 1)+theme_Publication()+ 
  xlab("Actual Biomass [t/ha]") +
  ylab("Predicted Biomass [t/ha]")+
  geom_text(aes(x=12, y=200, label= paste("RMSE= ", RMSE ), hjust=0))+
  geom_text(aes(x=12, y=185, label= paste("MAE= ", MAE ), hjust=0))+
  geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0)) # data= accuracy_H_with,


Map_1

# 
# ##save as png
# namepng <- paste('Plot1.png', sep = '')
# png(namepng, width=2800, height=2800, units='px', res=300)
# print(Map_1)
# dev.off()


# 
# ######alternative plot with viridis library (i.e. fancy colors)
# Map_2 <- ggplot(data = all.predictions,aes(x = actual, y = predictions)) + ###,color= V4
#   geom_abline(intercept = 0, slope = 1, colour = "black",linetype="dashed") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
#   geom_point(aes(size=actual, fill=actual), colour="black",pch=21,show.legend = FALSE) + #colour = "blue"
#   facet_wrap(~ model,ncol = 2)  +xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
#   ggtitle("Biomass, Predicted vs. Actual, by model")+ coord_fixed(ratio = 1)+theme_Publication()+ 
#   xlab("Actual Biomass [t/ha]") +
#   ylab("Predicted Biomass [t/ha]")+
#   geom_text(aes(x=12, y=200, label= paste("RMSE= ", RMSE ), hjust=0))+
#   geom_text(aes(x=12, y=185, label= paste("MAE= ", MAE ), hjust=0))+
#   geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0))+scale_fill_viridis() # data= accuracy_H_with,
# 
# 
# Map_2
# ##save as png
# 
# namepng <- paste('Plot2.png', sep = '')
# png(namepng, width=2800, height=2800, units='px', res=300)
# print(Map_2)
# dev.off()




#another one
ggplot(data = all.predictions,aes(x = actual, y = predictions)) +
  geom_smooth(method=lm, se=FALSE)+# data= accuracy_H_with, ###,color= V4
  geom_point(aes(size=actual, fill=actual), colour="black",pch=21,show.legend = FALSE) + #colour = "blue"
  geom_abline(intercept = 0, slope = 1, colour = "black") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
  facet_wrap(~ model,ncol = 2)  +xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
  ggtitle("Biomass, Predicted vs. Actual, by model")+ coord_fixed(ratio = 1)+theme_Publication()+ 
  xlab("Actual Biomass [t/ha]") +
  ylab("Predicted Biomass [t/ha]")+
  geom_text(aes(x=12, y=200, label= paste("RMSE= ", RMSE ), hjust=0))+
  geom_text(aes(x=12, y=185, label= paste("MAE= ", MAE ), hjust=0))+
  geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0))+scale_fill_viridis()

###...and another one


#another one

#another one
Map_4 <- ggplot(data = all.predictions,aes(x = actual, y = predictions)) +
  geom_smooth(method=lm)+# data= accuracy_H_with, ###,color= V4
  #colour = "blue"
  geom_abline(intercept = 0, slope = 1, colour = "black",linetype="dashed") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
  geom_point(aes(size=actual, fill=actual), colour="black",pch=21,show.legend = T) +
  facet_wrap(~ model,ncol = 2)  +xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
  ggtitle("AboveGround Biomass, Predicted vs. Actual- by model - Training Dataset")+ coord_fixed(ratio = 1)+theme_Publication()+ 
  ylab(expression("Predicted AboveGround Biomass ["~t~ha^-1~"]"))+
  xlab(expression("Actual AboveGround Biomass ["~t~ha^-1~"]"))+
  # 
  # xlab("Actual AboveGround Biomass [t haa^-1]") +
  # ylab("Predicted AboveGround Biomass [t ha^-1]")+
  geom_text(aes(x=12, y=200, label= paste("RMSE== ", RMSE, "~t~ha^-1"), hjust=0),parse=TRUE)+
  geom_text(aes(x=12, y=185, label= paste("MAE== ", MAE , "~t~ha^-1"), hjust=0),parse=TRUE)+
  geom_text(aes(x=12, y=170, label= paste("R^2== ", r2 ), hjust=0),parse=TRUE)+
  
  # geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0))+
  scale_fill_viridis()


Map_4
namepng <- paste('Plot4Train.png', sep = '')
png(namepng, width=2800, height=2800, units='px', res=300)
print(Map_4)
dev.off()
