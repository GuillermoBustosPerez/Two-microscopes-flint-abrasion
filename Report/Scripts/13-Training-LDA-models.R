load("Report/Data/Data-Both-micro-v3.RData")
library(tidyverse); library(caret)

#### Filter to get data from each microscope ####
Dino.Data <- Data %>% filter(Microscope == "Dinolite.Edge")

# Modal has 0 variance among Sensofar data. Remove
Senso.Data <- Data %>% filter(Microscope == "Sensofar.S.neox.090") %>% 
  select(-c(Modal))

#### Set validation set ####
trControl <- trainControl(method  = "repeatedcv",
                          verboseIter = TRUE,
                          number  = 10,
                          repeats = 50,
                          savePredictions = "final",
                          classProbs = TRUE)


#### Models from Sensofar data ####
# Set formula for complete set of variables
frmla <- as.formula(
  paste("Flake.Time", paste(colnames(Senso.Data[,2:15]), collapse = " + "), sep = " ~ "))

# LDA model on complete variables
set.seed(123)
Senso.LDA <- train(frmla, 
                   Senso.Data,
                   method = "lda",
                   preProc = c("center", "scale"),
                   trControl = trControl)

# LDA model on reduced number of variables
# Reduce number of variables
df.Test = cor(Senso.Data[,c(2:15)])^2
hc = findCorrelation(df.Test, 
                     cutoff = 0.9, 
                     names = TRUE, 
                     exact = TRUE) 
reduced_Data = Senso.Data %>% select(-all_of(hc))

frmla2 <- as.formula(
  paste("Flake.Time", paste(colnames(reduced_Data[,2:10]), collapse = " + "), sep = " ~ "))

set.seed(123)
Senso.Reduced.LDA <- train(frmla2, 
                           Senso.Data,
                           method = "lda",
                           preProc = c("center", "scale"),
                           trControl = trControl)

# LDA on Sensofar reduced dimensionality (PCA)
# PCA to reduce dimensionality
Senso.PCA <- prcomp(Senso.Data[,2:15], scale. = TRUE)
summary(Senso.PCA)

Senso.PC <- as.data.frame(Senso.PCA$x[,1:3])
Senso.PC$Flake.Time <- Senso.Data$Flake.Time

# LDA on PC scores
frmla3 <- as.formula(
  paste("Flake.Time", paste(colnames(Senso.PC[,1:3]), collapse = " + "), sep = " ~ "))

set.seed(123)
Senso.PCA.LDA <- train(frmla3, 
                       Senso.PC,
                       method = "lda",
                       preProc = c("center", "scale"),
                       trControl = trControl)

save(Senso.LDA, Senso.Reduced.LDA, Senso.PCA.LDA, 
     file = "Report/Models/Sensofar-LDA-models.RData")

#### Models from Dinolite data ####
frmla2.1 <- as.formula(
  paste("Flake.Time", paste(colnames(Dino.Data[,2:16]), collapse = " + "), sep = " ~ "))

# LDA model
set.seed(123)
Dino.LDA <- train(frmla2.1, 
                  Dino.Data,
                  method = "lda",
                  preProc = c("center", "scale"),
                  trControl = trControl)

# LDA model on reduced number of variables
# Reduce number of variables
df.Test = cor(Dino.Data[,c(2:16)])^2
hc = findCorrelation(df.Test, 
                     cutoff = 0.9, 
                     names = TRUE, 
                     exact = TRUE) 

reduced_Data = Dino.Data %>% select(-all_of(hc))
head(reduced_Data)

frmla2.2 <- as.formula(
  paste("Flake.Time", paste(colnames(reduced_Data[,2:7]), collapse = " + "), sep = " ~ "))

set.seed(123)
Dino.Reduced.LDA <- train(frmla2.2, 
                          Dino.Data,
                          method = "lda",
                          preProc = c("center", "scale"),
                          trControl = trControl)

# PCA to reduce dimensionality
Dino.PCA <- prcomp(Dino.Data[,2:15], scale. = TRUE)
summary(Dino.PCA)

Dino.PC <- as.data.frame(Dino.PCA$x[,1:3])
Dino.PC$Flake.Time <- Dino.Data$Flake.Time

# LDA on PC scores
frmla2.3 <- as.formula(
  paste("Flake.Time", paste(colnames(Dino.PC[,1:3]), collapse = " + "), sep = " ~ "))

set.seed(123)
Dino.PCA.LDA <- train(frmla2.3, 
                      Dino.PC ,
                      method = "lda",
                      preProc = c("center", "scale"),
                      trControl = trControl)

save(Dino.LDA, Dino.Reduced.LDA, Dino.PCA.LDA, 
     file = "Report/Models/Dinolite-LDA-models.RData")
rm(list = ls())