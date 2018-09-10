# Forest
Repository for technical data dealing with Above Ground Biomass (AGB) Modelling in Tanzania

A R script (i.e. Code_AGB_RS.R) allows predicting Above Ground Biomass.

We tested four predictive models to relate the remotely sensed parameters to AGB and evaluated their accuracies. Two models using inferential statistics (i.e. a generalised linear model and a generalised exponential model) and two models machine learning (i.e. a Random Forest model and a Support Vector Machine (SVM) model). The predictors of the models are image bands, their textures and spectral indices, while the response variable is the AGB.

Code_AGB_Training.R allows the user to display results using the training dataset.
