setwd("...")

#install.packages('ggplot2')
library(ggplot2)
#install.packages(c("ggalt","ggfortify","ggpubr","ggthemes","ztable"))
library(ggalt)
#install.packages("stringi")
library(stringi)
library(ggfortify)
library(dplyr)
library(tidyr)
library(scales)
#install.packages("ggpubr")
library(ggpubr)
library('scales')     
library(survival)
library(cowplot)
#install.packages("survminer")
library(survminer)
library(prodlim)
#install.packages("pec")
library(pec)

library(caret)
library(MASS)
#install.packages('gtools')
#install.packages('caTools')
#install.packages('pROC')
library(pROC)
#install.packages('glmnet')
library(glmnet)

### load and prepare data - merge drugbank pharmacological data and pubmed. 

FDA_clean <- read.csv2("data/FDA2018_by_drugbank2606.latest_clean.txt", 
                       sep = "\t", header = TRUE, dec = ",")
DB_cat <- read.csv2("C:/Users/Pranavi Shekhar/Desktop/RESOLVED2-master/tools/miner/drugbank_mining_category.latest.txt", 
                    sep = "\t", header = TRUE)
DB_target <- read.csv2("C:/Users/Pranavi Shekhar/Desktop/RESOLVED2-master/tools/miner/drugbank_mining_targets.latest.txt", 
                       sep = "\t", header = TRUE)

FDA_clean_ok <- FDA_clean %>% filter(DRUG_TO_KEEP_FOR_ANALYSIS==1)
FDA_clean_ok <- left_join(FDA_clean_ok, DB_cat, by=c("COMMON_DRUGBANK_ALIAS" = "X"))
FDA_clean_ok <- left_join(FDA_clean_ok, DB_target, by=c("COMMON_DRUGBANK_ALIAS" = "X"))

# select columns of interest
colnames(FDA_clean_ok[,1:30])
FDA_clean_ok_drugbank <- FDA_clean_ok[,c(1,17,18,27,31,33,37,39,42,43,47,49,58:ncol(FDA_clean_ok))]
colnames(FDA_clean_ok_drugbank[,1:30])
# keep only rows without NA
FDA_clean_ok_drugbank_noNA <- FDA_clean_ok_drugbank %>% filter(Clinical_activity_detected_as.following_recurrent_responder_OR_clinical_activity_OR_prolonged_stability != "NA")
FDA_clean_ok_drugbank_noNA <- FDA_clean_ok_drugbank_noNA %>% filter(DLT_identified_or_MTD_reached != "NA")
# check
sum(is.na(FDA_clean_ok_drugbank_noNA))
FDA_clean_ok_drugbank <-FDA_clean_ok_drugbank_noNA
nrow(FDA_clean_ok_drugbank)
dim(FDA_clean_ok_drugbank)

#### For review 1 use till here ##############################################################


# remove unrelevant categories
drops <- c("category.Antineoplastic.Agents", 
           "category.Myelosuppressive.Agents", 
           "category.Antineoplastic.and.Immunomodulating.Agents", 
           "category.Pharmacologic.Actions",             
           "category.Therapeutic.Uses", 
           "category.Chemical.Actions.and.Uses",
           "category.Immunosuppressive.Agents")
FDA_clean_ok_drugbank <- FDA_clean_ok_drugbank[ , !(names(FDA_clean_ok_drugbank) %in% drops)]
names(FDA_clean_ok_drugbank[,1:15])
dim(FDA_clean_ok_drugbank)


# Split the data in training and test sets 70:30
partition = createDataPartition(FDA_clean_ok_drugbank[,'FDA_APPROVAL'], 
                                times = 1, p = 0.7, list = FALSE)
training = FDA_clean_ok_drugbank[partition,] # Create the training sample
dim(training)
# save original training set that will be used for model fitting
write.table(training, "training_set.txt", sep = "\t")

test = FDA_clean_ok_drugbank[-partition,] # Create the test sample
dim(test)
# save original test set that will be used for model evaluation
write.table(test, "test_set.txt", sep = "\t")

# for further verifications: load directly the training and test sets used for model fit
# test <- read.csv2("data/test_set.txt" , sep="\t", header=TRUE, dec = ".")
# training <- read.csv2("data/training_set.txt", sep="\t", header=TRUE, dec = ".")
# test <- test[,2:1416]
# training <- training[,2:1416]

dim(test)
dim(training)

# View(training[1:10,])
# sum(is.na(test))
# sum(is.na(training))
# FDA_clean_ok$OLDEST_DATE_OF_PUBLICATION_manually_cured

# check date distribution in test and training

FDA_clean_ok_date_of_pub <- FDA_clean_ok[,c("COMMON_DRUGBANK_ALIAS","OLDEST_DATE_OF_PUBLICATION_manually_cured")]

test_date <- left_join(test, FDA_clean_ok_date_of_pub, by="COMMON_DRUGBANK_ALIAS")
training_date <- left_join(training, FDA_clean_ok_date_of_pub, by="COMMON_DRUGBANK_ALIAS")

test_date$date <- as.Date(test_date$OLDEST_DATE_OF_PUBLICATION_manually_cured, "%d/%m/%Y")
test_date$YEAR <- as.Date(cut(test_date$date,
                              breaks = "year"))

training_date$date <- as.Date(training_date$OLDEST_DATE_OF_PUBLICATION_manually_cured, "%d/%m/%Y")
training_date$YEAR <- as.Date(cut(training_date$date,
                                  breaks = "year"))




# 1. Fetaure selection

variables_matrix = as.matrix(training[,5:ncol(training)])
y <- Surv(training$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, training$FDA_APPROVAL)

# Can be used to get the same folds each time while comparing the model
fold=rep(1:100,length.out=nrow(training))

# Fit Cox regression using lasso with 100-fold cross-validation, in order to better features
cox.lasso=cv.glmnet(variables_matrix
                    ,y
                    ,family="cox",
                    foldid=fold, 
                    alpha = 1)#alpha = 1 means lasso



# Value of lambda that gives the minimum of the cross-validation curve
cox.lasso$lambda.min

# fit model with lambda min identified by CV. The above was done to identify minimum lambda.
model_final = glmnet(variables_matrix
                     ,y
                     ,family="cox"
                     ,lambda = cox.lasso$lambda.min
                     ,alpha = 1)

# Get estimated positive coefficients:
coefficients=coef(model_final, s=model_final$lambda.min)
coefficients_pos= coefficients[coefficients[,1]>0,]

# model_final_coef <- as.data.frame(coefficients_pos)
# model_final_coef$category <- rownames(model_final_coef)
# model_final_coef$category <- gsub("category.", "", model_final_coef$category)
# model_final_coef$category <- gsub("\\.", " ", model_final_coef$category)



# get positive coefficients and put it in a dataframe

#coefficients <- coef(model_final, s=model_final$lambda.min)
coefficients_nonzero =  coefficients[coefficients[,1]>0,]
fit_coef <- as.data.frame(coefficients_nonzero)

# save the model for further use
write.table(coefficients_nonzero, "Lasso-penalized Beta coefficients.txt", sep = "\t")

# for further validation steps: load the model
coefficients_nonzero <- read.csv2("Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE, dec = ".")
#coefficients_nonzero = read.csv2("model_coefficients.csv", sep=",",header=TRUE, dec = ".")
dim(coefficients_nonzero)
fit_coef <- as.data.frame(coefficients_nonzero)

# get features names
variables <- rownames(fit_coef)

# prepare test set with only selected features and prepare a matrix
test_select <- test %>% dplyr::select(variables)
test_mat <- as.matrix(test_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)

# compute approval score using HR
test$probs <- exp(test_mat%*%coefficients_nonzero)


# compute basic cindex
#install.packages("survcomp")
#BiocManager::install("survcomp")
library(survcomp)
ConcIndex<-concordance.index(x=test$probs, 
                             surv.time=Surv(test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                            test$FDA_APPROVAL), 
                             surv.event=test$FDA_APPROVAL, 
                             method="noether")

output <- data.frame(C.index=ConcIndex$c.index,
                     Lower=ConcIndex$lower,
                     Higher=ConcIndex$upper,
                     pvalue=ConcIndex$p.value)
output$a <- "["; output$b <- "-"; output$c <- "]"
output <- unite(output, 'C-index [95% Confidence interval]', c(C.index, a , Lower, b, Higher , c), sep = "", remove = T)
print(output)
write.table(output, "cindex.txt", sep = "\t")

# compute IPCW
library(pec)
#install.packages("htmlTable")
library(htmlTable)
#install.packages("sandwich")
library(sandwich)
library(survival)

# compute cindex
variables <- read.csv2("Lasso-penalized Beta coefficients.txt", sep="\t", header=TRUE)
features <- rownames(variables)
test_features = as.data.frame(test[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])

# inverse of the probability of censoring weigthed estimate of the concordance probability to adjust for right censoring
cindex <- pec::cindex(as.matrix(1/test$probs),
                      formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                     FDA_APPROVAL)~.,
                      data=test,
                      cens.model="marginal")
cindex


#2. Compare wiht EffTox
# clinical activity, complete response and toxicity
# select features & put variables into a matrix
# this model will be named "EffTOX" model

variables_matrix_alter = as.matrix(training[,c("Clinical_activity_detected_as.following_recurrent_responder_OR_clinical_activity_OR_prolonged_stability","Complete_response_reported","DLT_identified_or_MTD_reached")])

y <- Surv(training$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, training$FDA_APPROVAL)
fold=rep(1:100,length.out=nrow(training))

cox.lasso=cv.glmnet(variables_matrix_alter
                    ,y
                    ,family="cox",
                    foldid=fold, 
                    alpha = 1)

cox.lasso$lambda.min

model_final_alter = glmnet(variables_matrix_alter
                           ,y
                           ,family="cox"
                           ,lambda = cox.lasso$lambda.min
                           ,alpha = 1
)

coefficients=coef.glmnet(model_final_alter)
coefficients_pos=coefficients[coefficients[,1]!=0,]

write.table(as.data.frame(as.matrix(coefficients_pos)), "efficacity x tox model coefficients.txt", sep = "\t")

model_final_alter_coef <- as.data.frame(as.matrix(coefficients_pos))

model_final_alter_coef$category <- rownames(model_final_alter_coef)
model_final_alter_coef$category <- gsub("category.", "", model_final_alter_coef$category)
model_final_alter_coef$category <- gsub("\\.", " ", model_final_alter_coef$category)
model_final_alter_coef$category <- gsub("Clinical_activity_detected_as following_recurrent_responder_OR_clinical_activity_OR_prolonged_stability"
                                        , "Clinical activity", model_final_alter_coef$category)
model_final_alter_coef$category <- gsub("_", " ", model_final_alter_coef$category)
model_final_alter_coef$Coefficients <- model_final_alter_coef$s0



# compute approvsl score

coefficients_nonzero <- read.csv2("efficacity x tox model coefficients.txt", sep="\t", header=TRUE, dec = ".")
fit_coef <- as.data.frame(coefficients_nonzero)
variables <- rownames(fit_coef)
test_select <- test %>% dplyr::select(variables)
test_mat <- as.matrix(test_select)
coefficients_nonzero <- as.matrix(coefficients_nonzero)
test$probs_efftox <- exp(test_mat%*%coefficients_nonzero)


# compute basic cindex
library(survcomp)
ConcIndex_effTox<-concordance.index(x=test$probs_efftox, 
                                    surv.time=Surv(test$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                                   test$FDA_APPROVAL), 
                                    surv.event=test$FDA_APPROVAL, 
                                    method="noether")
ConcIndex_effTox$c.index

output_efftox <- data.frame(C.index=ConcIndex_effTox$c.index,
                            Lower=ConcIndex_effTox$lower,
                            Higher=ConcIndex_effTox$upper,
                            pvalue=ConcIndex_effTox$p.value)
output_efftox$a <- "["; output_efftox$b <- "-"; output_efftox$c <- "]"
output_efftox <- unite(output_efftox, 'C-index [95% Confidence interval]', c(C.index, a , Lower, b, Higher , c), sep = "", remove = T)
print(output_efftox)
write.table(output_efftox, "cindex_from_effXtox_model.txt", sep = "\t")

# compute IPCW
library(pec)
library(htmlTable)
library(sandwich)

variables <- read.csv2("efficacity x tox model coefficients.txt", sep="\t", header=TRUE)
features <- rownames(variables)
test_features = as.data.frame(test[,c(features,"DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN", "FDA_APPROVAL")])

# inverse of the probability of censoring weigthed estimate of the concordance probability to adjust for right censoring
cindex_IPCW_effTox <- pec::cindex(as.matrix(1/test$probs_efftox),
                                  formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, 
                                                 FDA_APPROVAL)~.,
                                  data=test,
                                  cens.model="marginal"
                                  #           splitMethod = "BootCv"
)
cindex_IPCW_effTox = as.data.frame(cindex_IPCW_effTox$AppCindex)



# 3. Compare with RF
test <- read.csv2("test_set.txt" , sep="\t", header=TRUE, dec = ".")
training <- read.csv2("data/tset.txt", sep="\t", header=TRUE, dec = ".")
test <- test[,2:1415]
training <- training[,2:1415]

dim(FDA_clean_ok_drugbank)

colnames(FDA_clean_ok_drugbank[,1:10])
FDA_clean_ok_drugbank_ok <- FDA_clean_ok_drugbank[,3:1415]
colnames(FDA_clean_ok_drugbank_ok[,1:10])
#install.packages("randomForestSRC")

library(randomForestSRC)


train <- sample(1:nrow(FDA_clean_ok_drugbank), 
                round(nrow(FDA_clean_ok_drugbank_ok) * 0.70))

RF <- rfsrc(Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN,
                 FDA_APPROVAL) ~ ., data = FDA_clean_ok_drugbank_ok[train, ], 
            ntree = 500, importance = T)

RF_predictions <- predict(RF_RESOLVED2, FDA_clean_ok_drugbank_ok[-train , ])

print(RF_predictions)

cindex_IPCW_RF <- pec::cindex(RF_RESOLVED2,
                              formula = Surv(DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN, FDA_APPROVAL)~.,data=FDA_clean_ok_drugbank_ok[-train , ],cens.model="marginal")

cindex_IPCW_RF

###########################  REVIEW 3   ########s######################
#THRESHOLD CALIBRATION
score_model = function(df, threshold){
  df$score = ifelse(df$probs < threshold, 'predicted non-approved', 'predicted approved')
  df
}
# first compute predicted probabilities on training set
fit_coef <- as.data.frame(coefficients_nonzero)
rownames(fit_coef)=coefficients_nonzero$Feature
 
# get features names
variables <- rownames(fit_coef)
# prepare test set with only selected features and prepare a matrix
training_select <- training %>% dplyr::select(variables)
training_mat <- as.matrix(training_select)
coefficients_nonzero <- as.matrix(fit_coef$`coefficients_nonzero$Coeff`)
# compute predicted probabilities
training$probs <- exp(training_mat%*%coefficients_nonzero)

lseq <- function(from=1, to=100000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

log_sequence <- lseq(from = min(training$probs)+0.1, to = 20, length.out = 50)

data <- training

output_threshold <- data.frame(matrix(ncol = 3, nrow = length(log_sequence)))
colnames(output_threshold) <- c("threshold","p","chisq")
output_threshold$threshold <- log_sequence

for (i in 1:length(log_sequence)) {
  
  data = score_model(data, output_threshold$threshold[i])
  
  FDA_FS_threshold_i <- survdiff(Surv(data$DELAY_FROM_OLDEST_PMID_TO_FDA_APPROVAL_or_DDN/12, 
                                      data$FDA_APPROVAL) ~ data$score)
  
  output_threshold[i,2] <- pchisq(FDA_FS_threshold_i$chisq, length(FDA_FS_threshold_i$n)-1, lower.tail = FALSE)
  output_threshold[i,3] <- FDA_FS_threshold_i$chisq
}
sorted_by_p <- output_threshold[order(output_threshold$p),]
optimal_threshold = sorted_by_p$threshold[1]


# Code to test
test_select <- test %>% dplyr::select(variables)
test_mat <- matrix(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1),nrow=1)
coefficients_nonzero <- as.matrix(coefficients_nonzero)

# compute HR
myprobs <- exp(test_mat%*%coefficients_nonzero)

