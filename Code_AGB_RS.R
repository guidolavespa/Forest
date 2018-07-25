## *********************************************************************************************************************************************************
##              Repository for technical data dealing with Above Ground Biomass (AGB) Modelling in Tanzania (Version 1.0)
## *********************************************************************************************************************************************************
## 
##  * JRC of the European Commission
##  *
##  * Purpose:  A R script allows predicting Above Ground Biomass.
##  *
##  *           We tested four predictive models to relate the remotely sensed parameters to AGB and evaluated their accuracies.
##  *           Two models using inferential statistics (i.e. a generalised linear model and a generalised exponential model) 
##  *           and two models machine learning (i.e. a Random Forest model and a Support Vector Machine (SVM) model). 
##  *           The predictors of the models are image bands, their textures and spectral indices, while the response variable 
##  *           is the AGB.
##  *           
##  *           - Updates to the Forest Article:           https://doi.org/TOBEDEFINED
##  *
##  *
##  *
##  * Author:   Guido Ceccherini
##  * Email:    guido.ceccherini@gmail.com, guido.ceccherini@ec.europa.eu
##
##



###load libraries
library(randomForest)
library(e1071)
library(ggthemes)
library(grid)
library(ggplot2)
library(tidyr)
library(viridis)

#####load dataset
DATA_BIOMASS <- read.csv('DataBiomassRS.csv')

##### reproducibility number
set.seed(2511) 


#### remove the first column (sequence from 1 to 512)
DATA_BIOMASS <- DATA_BIOMASS[,2:27]

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
test.pred.lin <- exp(predict(log.reg_H_with,test))-1

# ...and evaluate the accuracy
RMSE.log.reg_H_with <- sqrt(mean((test.pred.lin-test$BIOMASS)^2))
RMSE.log.reg_H_with


MAE.log.reg_H_with <- mean(abs(test.pred.lin-test$BIOMASS))
MAE.log.reg_H_with

Rel.RMSEexp <- sqrt(mean((test.pred.lin-test$BIOMASS)^2)) / diff(range(test$BIOMASS))
Rel.RMSEexp


R2exp <- 1 - sum((test$BIOMASS-test.pred.lin)^2)/sum((test$BIOMASS-mean(test$BIOMASS))^2)

##################
# 2 Linear Model
##################

# Create a multiple linear regression model using the training data
lin.reg_H_with <- glm(BIOMASS ~ ., data = train)

summary(lin.reg_H_with )

# Apply the model to the testing data (i.e., make predictions) ...
test.pred.lin_ <- (predict(lin.reg_H_with,test))-1

# ...and evaluate the accuracy
RMSE.lin.reg_H_with <- sqrt(mean((test.pred.lin_-test$BIOMASS)^2))
RMSE.lin.reg_H_with


MAE.lin.reg_H_with <- mean(abs(test.pred.lin_-test$BIOMASS))
MAE.lin.reg_H_with


Rel.RMSElin <- sqrt(mean((test.pred.lin_-test$BIOMASS)^2)) / diff(range(test$BIOMASS))
Rel.RMSElin

R2lin <- 1 - sum((test$BIOMASS-test.pred.lin_)^2)/sum((test$BIOMASS-mean(test$BIOMASS))^2)


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



test2 <- rbind(test, train)

# As usual, predict and evaluate on the test set
test.pred.forest <- predict(rf_H_with,test2)

test.pred.forest <- test.pred.forest[1:nrow(test)]

RMSE.forest_H_with <- sqrt(mean((test.pred.forest-test$BIOMASS)^2))
RMSE.forest_H_with


MAE.forest_H_with <- mean(abs(test.pred.forest-test$BIOMASS))
MAE.forest_H_with


# Rel.RMSE <- sqrt(mean((test.pred.forest-test$BIOMASS)^2)) / sd(test$BIOMASS)
Rel.RMSErf <- sqrt(mean((test.pred.forest-test$BIOMASS)^2)) / diff(range(test$BIOMASS))
Rel.RMSErf
### ref Remote Sens. 2016, 8, 583; doi:10.3390/rs8070583 
R2RF <- 1 - sum((test$BIOMASS-test.pred.forest)^2)/sum((test$BIOMASS-mean(test$BIOMASS))^2)


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
predictedY <- predict(SVM_H, test2)
predictedY <- predictedY[1:nrow(test)]
plot(test$BIOMASS, predictedY, col = "red", pch=4)

error <- train$BIOMASS - predictedY
svrPredictionRMSE <- rmse(error)
svrPredictionRMSE

tuneResult <- tune(svm, BIOMASS ~., data = train, kernel="radial", ranges=list(cost=10^(-1:2), gamma=c(.5,1,2)))

print(tuneResult)
# Draw the tuning graph
plot(tuneResult)

##tune the SVM model
tunedModel <- tuneResult$best.model
test.pred.svm <- predict(tunedModel, test2) 
test.pred.svm <- test.pred.svm[1:nrow(test)]
error <- test$BIOMASS - test.pred.svm  


##compute the performance metrics
tunedModelRMSE <- rmse(error)
tunedModelRMSE


RMSE.svm_H_with <- sqrt(mean((test.pred.svm-test$BIOMASS)^2))
RMSE.svm_H_with


MAE.svm_H_with <- mean(abs(test.pred.svm-test$BIOMASS))
MAE.svm_H_with

Rel.RMSEsvm <- sqrt(mean((test.pred.svm-test$BIOMASS)^2)) / diff(range(test$BIOMASS))
Rel.RMSEsvm
# SVM
R2svm <- 1 - sum((test$BIOMASS-test.pred.svm)^2)/sum((test$BIOMASS-mean(test$BIOMASS))^2)




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
all.predictions <- data.frame(actual = test$BIOMASS,
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

all.predictions$RMSE <- rep(accuracy_H_with$RMSE, each=77)
all.predictions$MAE <- rep(accuracy_H_with$MAE, each=77)
all.predictions$relRMSE <- rep(round(accuracy_H_with$relRMSE, digit=3), each=77)



####first plot(base theme)

Map_1 <- ggplot(data = all.predictions,aes(x = actual, y = predictions)) + ###,color= V4
  geom_point(colour = "blue") + #colour = "blue"
  geom_abline(intercept = 0, slope = 1, colour = "black") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
  facet_wrap(~ model,ncol = 2)  +xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
  ggtitle("Biomass, Predicted vs. Actual, by model")+ coord_fixed(ratio = 1)+theme_Publication()+ 
  xlab("Actual Biomass [t/ha]") +
  ylab("Predicted Biomass [t/ha]")+
  geom_text(aes(x=12, y=200, label= paste("RMSE= ", RMSE ), hjust=0))+
  geom_text(aes(x=12, y=185, label= paste("MAE= ", MAE ), hjust=0))+
  geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0)) # data= accuracy_H_with,


Map_1


##save as png
namepng <- paste('Plot1.png', sep = '')
png(namepng, width=2800, height=2800, units='px', res=300)
print(Map_1)
dev.off()



######alternative plot with viridis library (i.e. fancy colors)
Map_2 <- ggplot(data = all.predictions,aes(x = actual, y = predictions)) + ###,color= V4
  geom_abline(intercept = 0, slope = 1, colour = "black",linetype="dashed") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
  geom_point(aes(size=actual, fill=actual), colour="black",pch=21,show.legend = FALSE) + #colour = "blue"
  facet_wrap(~ model,ncol = 2)  +xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
  ggtitle("Biomass, Predicted vs. Actual, by model")+ coord_fixed(ratio = 1)+theme_Publication()+ 
  xlab("Actual Biomass [t/ha]") +
  ylab("Predicted Biomass [t/ha]")+
  geom_text(aes(x=12, y=200, label= paste("RMSE= ", RMSE ), hjust=0))+
  geom_text(aes(x=12, y=185, label= paste("MAE= ", MAE ), hjust=0))+
  geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0))+scale_fill_viridis() # data= accuracy_H_with,


Map_2
##save as png

namepng <- paste('Plot2.png', sep = '')
png(namepng, width=2800, height=2800, units='px', res=300)
print(Map_2)
dev.off()




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
Map_4 <- ggplot(data = all.predictions,aes(x = actual, y = predictions)) +
  geom_smooth(method=lm)+# data= accuracy_H_with, ###,color= V4
   #colour = "blue"
  geom_abline(intercept = 0, slope = 1, colour = "black",linetype="dashed") + #geom_vline(xintercept = 23, colour = "green", linetype = "dashed") +
  geom_point(aes(size=actual, fill=actual), colour="black",pch=21,show.legend = T) +
  facet_wrap(~ model,ncol = 2)  +xlim (c(0,max(all.predictions$actual))) +ylim (c(0,max(all.predictions$actual))) +
  ggtitle("AboveGround Biomass, Predicted vs. Actual, by model")+ coord_fixed(ratio = 1)+theme_Publication()+ 
  xlab("Actual AboveGround Biomass [t/ha]") +
  ylab("Predicted AboveGround Biomass [t/ha]")+
  geom_text(aes(x=12, y=200, label= paste("RMSE= ", RMSE ), hjust=0))+
  geom_text(aes(x=12, y=185, label= paste("MAE= ", MAE ), hjust=0))+
  # geom_text(aes(x=12, y=170, label= paste("relRMSE= ", relRMSE ), hjust=0))+
  scale_fill_viridis()


Map_4
namepng <- paste('Plot4.png', sep = '')
png(namepng, width=2800, height=2800, units='px', res=300)
print(Map_4)
dev.off()




##################################################################################################################
#  RANDOM FOREST permultation to find the statistical distribution of the variable importance
#################################################################################################################


##### use the entire dataset: DATA_BIOMASS
rf_H_with <- randomForest(BIOMASS ~., data = DATA_BIOMASS,importance = TRUE, proximity = T,ntree=500)


# Using the importance()  function to calculate the importance of each variable
imp_H_with <- as.data.frame(sort(importance(rf_H_with)[,1],decreasing = TRUE),optional = T,col.names=T)
names(imp_H_with) <- "% Inc MSE"
imp_H_with


actual <- DATA_BIOMASS$BIOMASS
predicted <- unname(predict(rf_H_with, DATA_BIOMASS))
R2 <- 1 - (sum((actual-predicted)^2)/sum((actual-mean(actual))^2))


##### start permutations

var_imp_range <- as.data.frame(importance(rf_H_with)[,1])
var_impSD_range <- as.data.frame(rf_H_with$importanceSD)
for (i in 1:100){
  rf<-randomForest(BIOMASS ~., data = DATA_BIOMASS,importance = TRUE, proximity = T,ntree=500, importanceSD=T, keep.inbag=T, mse=T, rsq=T)
  imp.df<- as.data.frame(importance(rf))
  var_imp_range <- cbind(var_imp_range, imp.df[,1])
  impSD.df<- as.data.frame(rf$importanceSD)
  var_impSD_range<- cbind(var_impSD_range, impSD.df[,1])
}


# 

var_imp_range

## first plot to have a general idea
boxplot(var_imp_range, las = 2, horizontal = T)



#reorder the data and comput the statistical distribution needed for the boxplot
newdata <- var_imp_range[order(var_imp_range$`importance(rf_H_with)[, 1]`),] 
##create a dataframe need for printing the results
newdata<- newdata[3:26,]

df <- data.frame(
  name= rownames(newdata),
  y0 = apply(newdata,1,min),
  y25 = apply(newdata, 1, quantile, probs = c(0.25),  na.rm = TRUE),
  y50 = apply(newdata,1,median),  
  y75 = apply(newdata, 1, quantile, probs = c(0.75),  na.rm = TRUE),
  y100 = apply(newdata,1,max)   
)

### take the most important variables
df <- df[5:23,]

### change names: reflectances is indicated with both LX and RX.
levels(df$name)[levels(df$name)=="L5SD5"] <- "R5SD5"
levels(df$name)[levels(df$name)=="L3SD3"] <- "R3SD3"
levels(df$name)[levels(df$name)=="L4SD4"] <- "R4SD4"
levels(df$name)[levels(df$name)=="L1SD1"] <- "R1SD1"
levels(df$name)[levels(df$name)=="L2SD2"] <- "R2SD2"
levels(df$name)[levels(df$name)=="Rel_shad"] <- "Shadow Index"




##plot boxplot using the ggplot2 library
Map_3 <- ggplot(df, aes(x=reorder(name,y50), y=y50, fill=y50 )) +
  geom_boxplot(
    aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
    stat = "identity"
  )+  coord_flip()+
  scale_fill_viridis(direction = 1, na.value = "gray",  name="% Inc MSE")+  theme( ### trans="reverse
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(linetype="dashed",colour="grey"),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("Variable") +
  ylab("% Inc MSE") +
  ggtitle("RF Variable Importance AboveGround Biomass")

Map_3

##save the results
Map_3
ggsave("Plot3.png", width = 15, height = 20, units = "cm",dpi=300)







##################################################################################################################
#  RANDOM FOREST Split and compare quantiles
#################################################################################################################
## from https://github.com/laresbernardo/lares/blob/ff7378b9f0ba74fc7ee66067f13945d101bdff34/R/model_plots.R 
mplot_splits <- function(tag, score, splits = 5, subtitle = NA, model_name = NA, facet = NA, 
                         save = FALSE, subdir = NA, file_name = "viz_splits.png") {
  
  require(ggplot2)
  require(dplyr)
  require(RColorBrewer)
  
  if (length(tag) != length(score)) {
    message("The tag and score vectors should be the same length.")
    stop(message(paste("Currently, tag has",length(tag),"rows and score has",length(score))))
  }
  
  if (splits > 10) {
    stop("You should try with less splits!")
  }
  
  df <- data.frame(tag, score, facet)
  npersplit <- round(nrow(df)/splits)
  
  # For continuous tag values
  if (length(unique(tag))) {
    names <- df %>% 
      mutate(tag = as.numeric(tag), 
             quantile = ntile(tag, splits)) %>% group_by(quantile) %>%
      summarise(n = n(), 
                max_score = round(max(tag), 1), 
                min_score = round(min(tag), 1)) %>%
      mutate(quantile_tag = paste0(quantile," (",min_score,"-",max_score,")"))
    df <- df %>% 
      #mutate(score = score/100, tag = tag/100) %>%
      mutate(quantile = ntile(tag, splits)) %>%
      left_join(names, by = c("quantile")) %>% mutate(tag = quantile_tag) %>% 
      select(-quantile, -n, -max_score, -min_score)
    
  } else {
    # For categorical tag values
    names <- df %>% 
      mutate(quantile = ntile(score, splits)) %>% group_by(quantile) %>%
      summarise(n = n(), 
                max_score = round(100 * max(score), 1), 
                min_score = round(100 * min(score), 1)) %>%
      mutate(quantile_tag = paste0(quantile," (",min_score,"-",max_score,")")) 
  }
  
  p <- df %>% 
    mutate(quantile = ntile(score, splits)) %>% 
    group_by(quantile, facet, tag) %>% tally() %>%
    ungroup() %>% group_by(facet, tag) %>% 
    arrange(desc(quantile)) %>%
    mutate(p = round(100*n/sum(n),2),
           cum = cumsum(100*n/sum(n))) %>%
    left_join(names, by = c("quantile")) %>%
    ggplot(aes(x = as.character(tag), y = p, label = as.character(p),
               fill = as.character(quantile_tag))) + theme_minimal() +
    geom_col(position = "stack") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), check_overlap = TRUE) +
    xlab("Quantiles AGB") + ylab("Total Percentage by AGB") +
    guides(fill = guide_legend(title=paste0("Predicted Biomass"))) +
    labs(title = "AboveGround Biomass Predicted vs Actual, by Quantiles") +
    scale_fill_brewer(palette = "Spectral")
  
  if(!is.na(subtitle)) {
    p <- p + labs(subtitle = subtitle)
  }  
  
  if(!is.na(model_name)) {
    p <- p + labs(caption = model_name)
  }
  
  if(!is.na(facet)) {
    p <- p + facet_grid(. ~ facet, scales = "free")
  }  
  
  if (!is.na(subdir)) {
    dir.create(file.path(getwd(), subdir))
    file_name <- paste(subdir, file_name, sep="/")
  }
  
  if (save == TRUE) {
    p <- p + ggsave(file_name, width = 6, height = 6)
  }
  
  return(p)
  
}
mplot_splits(tag = test$BIOMASS, 
                    score = test.pred.forest,
                    split = 4)
ggsave("AGB_quantiles.png", width = 25, height = 20, units = "cm",dpi=300)
