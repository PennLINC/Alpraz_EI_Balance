## Load Libraries
require(ggplot2)
require(dplyr)
library(e1071)
library(randomForest)
library(matrixTests)
library(doParallel)
theme_set(theme_classic())
setwd("/cbica/projects/alpraz_EI/scripts/")

extract_matrix2 <- function(subid,sesid,atlasname,gsr,subdivide=F, subdivision=NULL){
  # This function loads the connectivity matrix for a subj, filters features (if needed; subdivide=T and specify subdivision), and converts to a vector
  if (gsr == "GSR") {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK_GSR/%s/%s/%s_CC_GSR_000.netcc",subid,sesid,atlasname)
  } else if (gsr == "rest_GSR") {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/NO_TASK_REGRESS_GSR/%s/%s/%s_CC_rest_GSR_000.netcc",subid,sesid,atlasname)
  } else if (gsr == "rest_noGSR") {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/NO_TASK_REGRESS/%s/%s/%s_CC_rest_noGSR_000.netcc",subid,sesid,atlasname)
  } else if (gsr == "PNC") {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/PNC/IDEMO_ACOMPCOR_GSR/%s/%s/%s_CC_000.netcc",subid,sesid,atlasname)
  }else {
    fname <- sprintf("/cbica/projects/alpraz_EI/data/TASK/%s/%s/%s_CC_noGSR_000.netcc",subid,sesid,atlasname)
  }
  print(fname)
  thisCC <- as.matrix(read.table(fname,header = T,skip = 4))
  # thisCC_vector = as.vector(thisCC)
  
  # Select specific features if desired.
  if (subdivide==T) {
    templateCC <- thisCC
    
    # Get labels
    labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_labels.csv',atlasname),header=F)
    if (atlasname %in% c("schaefer400x7_aal","schaefer200x7_aal")) {
      community_labels <- read.csv(sprintf('/cbica/projects/alpraz_EI/input/atlases/%s_node_labels.csv',
                                           gsub(x = atlasname,pattern = "_aal",replacement = ""))
                                   ,header=T)
      if (atlasname == "schaefer200x7_aal") {
        trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef200_transmodality7.csv',header = F)
      } else{
        trans_labels <- read.csv('/cbica/projects/alpraz_EI/input/atlases/Schaef400_transmodality7.csv',header = F)
    }

      colnames(trans_labels) <- c('name','transmodality')
      trans_labels$name <- gsub(trans_labels$name,pattern = "7Networks_",replacement = "")
      trans_labels$nums = as.numeric(rownames(trans_labels))

      
      mat_labels <- data.frame(nums = colnames(templateCC))
      mat_labels$nums <- as.numeric(gsub(mat_labels$nums,pattern="X",replace=""))
      mat_labels <- mat_labels %>% left_join(labels,by = c("nums"="V1"))
      mat_labels <- mat_labels %>% left_join(community_labels,by = c("nums"="index"))
      mat_labels <- mat_labels %>%
        mutate(community_name = as.character(community_name),V2 = as.character(V2),node_names=as.character(node_names))
      mat_labels$community_name[is.na(mat_labels$community_name)]="AAL_subcortex"
      mat_labels$node_names[is.na(mat_labels$node_names)]= mat_labels$V2[is.na(mat_labels$node_names)]
      
      # Check if we are using transmodality scores
      if (subdivision %in% c("transmodal","unimodal","transmodal25","unimodal25")) {
        #add transmodality scores
        # mat_labels <- mat_labels %>%
        #   left_join(trans_labels,by = c("node_names"="name"))
        mat_labels <- mat_labels %>%
          left_join(trans_labels,by = "nums")
        
        #cut transmodality into quantiles
        if (grepl(x=subdivision,pattern = "25")) {
          cuts <- quantile(mat_labels$transmodality,c(0,1/4,3/4,1),na.rm = T)
        } else{
          cuts <- quantile(mat_labels$transmodality,c(0,1/3,2/3,1),na.rm = T)
        }
        mat_labels <- mat_labels %>% 
          mutate(trans_category=cut(transmodality, breaks=cuts, labels=c("unimodal","middle","transmodal"),include.lowest = T))
        # Set names to transmodality names

        colnames(thisCC)=mat_labels$trans_category
        rownames(thisCC)=mat_labels$trans_category
        
        
        if (grepl(x=subdivision,pattern = "transmodal")) {
          communities="transmodal"
          commCC <- thisCC[rownames(thisCC) %in% communities,]
          vals <- as.vector(commCC)
        } else if (grepl(x=subdivision,pattern = "unimodal")) {
          communities="unimodal"
          commCC <- thisCC[rownames(thisCC) %in% communities,]
          vals <- as.vector(commCC)
          length(vals)
        } 
      } else if (subdivision=="regional"){
        longCC <- as.data.frame(thisCC)%>%
          mutate(row_label = colnames(thisCC))%>% 
          pivot_longer(cols = contains("X"))
        longCC$subid = subid
        longCC$sesid = sesid
        vals=longCC
        
      } else {
        # Set names to community names
        colnames(thisCC)=mat_labels$community_name
        rownames(thisCC)=mat_labels$community_name
        if (subdivision=="association") {
          communities=c("default","frontoparietalControl")  #"salienceVentralAttention",
          commCC <- thisCC[rownames(thisCC) %in% communities,]
          vals <- as.vector(commCC)
          } else if (subdivision=="sensorimotor") {
            communities=c("visual","somatomotor")
            commCC <- thisCC[rownames(thisCC) %in% communities,]
            vals <- as.vector(commCC)
            length(vals)
          } else if (subdivision=="association_plus") {
            communities=c("default","frontoparietalControl","dorsalAttention")
            commCC <- thisCC[rownames(thisCC) %in% communities,]
            vals <- as.vector(commCC)
            length(vals)
          }
      } 
    }
  } 
  else{
      thisCC_vector = thisCC[upper.tri(thisCC,diag = F)]
      vals <-thisCC_vector
  }
  return(vals)
}


SVM_2class <- function(df,folds,feature_selection = F,feature_proportion = .1){
  #SVM 2-class classifier
  
  # set up folds for CV
  if (folds == "LOO") {
    num_folds = length(unique(df$subid))
  } else {
    num_folds = folds
  }
  foldIdxs <- data.frame(subid=unique(df$subid))
  foldIdxs$foldID <- row_number(foldIdxs$subid)
  foldIdxs$foldID <- ntile(foldIdxs$foldID,num_folds)
  # foldIdxs$subid <- sample(foldIdxs$subid)
  pred_data<-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("subid", "sesid", "drug","model.pred"))
  
  W = matrix(nrow = num_folds,ncol = length(4:dim(df)[2]))
  colnames(W)<-colnames(df)[4:dim(df)[2]]
  
  # Loop over folds
  for (fold in 1:num_folds) {
    #training data
    trainingIDs <- as.matrix(foldIdxs %>% filter(foldID != fold) %>% select(subid))
    trainingIndex <- df$subid %in% trainingIDs # indices for training subs
    trainingData <- df[trainingIndex, 3:dim(df)[2] ] # Training data. Take columns 3:end (Removes subid and sesid).
    
    #testing data
    testData <- df[!trainingIndex, 4:dim(df)[2]] # test data. Take columns 4:end (Removes subid sesid drug).
    testLabels <- data.frame(df[!trainingIndex,c(1:3) ]) # Labels for test data
    
    #feature extraction if needed
    if (feature_selection == TRUE) {
      feature_extracted_data <- featureExtraction(trainingData,feature_proportion = feature_proportion)
      trainingData <- feature_extracted_data[[1]]
      feature_names <- feature_extracted_data[[2]]
      
      testData <- testData %>% select(feature_names)
      cat(sprintf("Retained %d features\n",length(feature_names)))
    }
    
    # svm
    x <- as.matrix(trainingData[, 2:dim(trainingData)[2]])
    y <- as.factor(as.matrix(trainingData[,1]))
    svm.model <- svm(x =x, y = y, 
                     cost = 100, kernel = "linear",type = "C-classification",scale = F)
    svm.pred <- predict(svm.model, as.matrix(testData))
    
    W[fold,colnames(x)] <- t(svm.model$coefs) %*% svm.model$SV
    w <- t(svm.model$coefs) %*% svm.model$SV
    num_features <- dim(x)[2]
    decisionValues <- w %*% t(testData)-svm.model$rho
    distance <- decisionValues/norm(w)
    testLabels$decisionValues <- t(decisionValues)
    testLabels$distance <- t(distance)
    testLabels$model.pred = svm.pred
    pred_data <- rbind(pred_data,testLabels)
  }
  final_data<-df[, 3:dim(df)[2] ]
  x <- as.matrix(final_data[, 2:dim(final_data)[2]])
  y <- as.factor(as.matrix(final_data[,1]))
  svm.model <- svm(x =x, y = y, 
                   cost = 100, kernel = "linear",type = "C-classification",scale = F)
  svm_results <- list(pred_data,W,svm.model)
  return(svm_results)
}

RF_2class <- function(df,folds,feature_selection = F,feature_proportion = .1){
  # Random forest classifier
  
  # set up folds for CV
  if (folds == "LOO") {
    num_folds = length(unique(df$subid))
  } else {
    num_folds = folds
  }
  foldIdxs <- data.frame(subid=unique(df$subid))
  foldIdxs$foldID <- row_number(foldIdxs$subid)
  foldIdxs$foldID <- ntile(foldIdxs$foldID,num_folds)
  # foldIdxs$subid <- sample(foldIdxs$subid)
  pred_data<-setNames(data.frame(matrix(ncol = 4, nrow = 0)), c("subid", "sesid", "drug","model.pred"))
  
  W = matrix(nrow = num_folds,ncol = length(4:dim(df)[2]))
  colnames(W)<-colnames(df)[4:dim(df)[2]]
  
  # Loop over folds
  for (fold in 1:num_folds) {
    #training data
    trainingIDs <- as.matrix(foldIdxs %>% filter(foldID != fold) %>% select(subid))
    trainingIndex <- df$subid %in% trainingIDs # indices for training subs
    trainingData <- trainingData <- df[trainingIndex, 3:dim(df)[2] ] # Training data. Take columns 3:end (Removes subid and sesid).
    
    #testing data
    testData <- df[!trainingIndex, 4:dim(df)[2]] # test data. Take columns 4:end (Removes subid sesid drug).
    testLabels <- data.frame(df[!trainingIndex,c(1:3) ]) # Labels for test data
    
    #feature extraction if needed
    if (feature_selection == TRUE) {
      feature_extracted_data <- featureExtraction(trainingData,feature_proportion = feature_proportion)
      trainingData <- feature_extracted_data[[1]]
      feature_names <- feature_extracted_data[[2]]
      
      testData <- testData %>% select(feature_names)
      cat(sprintf("Retained %d features\n",length(feature_names)))
    }
    
    # randomForest
    x <- as.matrix(trainingData[, 2:dim(trainingData)[2]])
    y <- as.factor(as.matrix(trainingData[,1]))
    rf.model <- randomForest(x = x, y = y,importance = T,mtry = 2*sqrt(dim(x)[2]))
    rf.pred <- predict(rf.model, as.matrix(testData))
    
    imp <- importance(rf.model)
    W[fold,colnames(x)] <- imp[,"MeanDecreaseAccuracy"]
    
    testLabels$model.pred = rf.pred
    # print(sprintf("Fold %d, acc = %1.3f",fold,sum(testLabels$drug==testLabels$model.pred)/nrow(testLabels)))
    pred_data <- rbind(pred_data,testLabels)
    
  }
  rf_results <- list(pred_data,W)
  return(rf_results)
}

run_model <- function(df,folds,feature_selection = F,feature_proportion = .1,permutation_test = F, num_permutations = 10000,type = "svm"){
  ## This sends the data to the desired classifier and other functions
  if (feature_selection == T) {
    # cat(sprintf("\nPerforming feature extraction: extracting the top %1.3f features\n",feature_proportion))
  }

  if (type == "svm") {
    print("Performing classification using SVM")
    model_results <- SVM_2class(df,folds,feature_selection = feature_selection,feature_proportion = feature_proportion)
  }else if (type == "rf") {
    print("Performing classification using Random Forest")
    model_results <- RF_2class(df,folds,feature_selection = feature_selection,feature_proportion = feature_proportion)
  }
  pred_data <- model_results[[1]]
  W_folds <- model_results[[2]]
  svm.model <- model_results[[3]]
  num_features <- sum(!is.na(W_folds[1,]))
  meanW <- colMeans(W_folds,na.rm=T)
  if (feature_selection==T) {
    feature_counts <- colSums(!is.na(W_folds))
    W <- rbind(meanW,feature_counts)
  } else{
    W <- meanW
  }
  
  accuracy <- sum(pred_data$model.pred==pred_data$drug)/dim(pred_data)[1]
  b <- binom.test(sum(pred_data$model.pred==pred_data$drug),dim(pred_data)[1],.5)
  b$pred_data <- pred_data
  
  if (permutation_test == T) {
    print(sprintf("Permuting %d times...",num_permutations))
    ptm = proc.time()
    cl<-makeCluster(12)
    registerDoParallel(cl)
    nw <- getDoParWorkers()
    # # perms <- shuffleSet(n = dim(df)[1],nset = num_permutations)
    perm_acc <- matrix(nrow = num_permutations)
    perm_W <- matrix(nrow = num_permutations,ncol = num_features)
    perm_list <- list()
    perm_list = foreach(perm_chunk = idiv(num_permutations,chunks = nw),
                        .combine=c,
                        .export = c("featureExtraction","SVM_2class"),
                        .packages = c("dplyr","e1071","matrixTests")) %dopar% {
                          # pacc = numeric(perm_chunk)
                          # pW = matrix(nrow = perm_chunk,ncol = num_features)
                          perm_result_list=list()
                          for (p in 1:perm_chunk){
                            # thisPerm <- perms[p,]
                            thisDF <- df
                            # thisDF$drug <- df$drug[thisPerm]
                            permuted <- df %>% select(subid,drug) %>% group_by(subid) %>% mutate(perm_drug=sample(drug))
                            thisDF$drug <- permuted$perm_drug
                            perm_pred_result <- SVM_2class(thisDF,folds,feature_selection = feature_selection,feature_proportion = feature_proportion)
                            perm_pred_data <- perm_pred_result[[1]]
                            perm_W_folds <- perm_pred_result[[2]]
                            perm_num_features <- sum(!is.na(perm_W_folds[1,]))
                            perm_meanW <- colMeans(perm_W_folds,na.rm=T)
                            pW <- perm_meanW
                            pacc<-sum(perm_pred_data$model.pred==perm_pred_data$drug)/dim(perm_pred_data)[1]
                            perm_result_list[p]=list(list(acc=pacc,W=pW))
                          }
                          
                          perm_result_list
                        }
    print("done")
    stopCluster(cl)
    
    print(proc.time()-ptm)
    perm_df <- as.data.frame(perm_list)
    perm_acc <- perm_df %>% select(contains("acc")) %>% summarize_all(mean)
    perm_plot <- ggplot(data = data.frame(perm_acc=t(perm_acc)),aes(x = perm_acc))+geom_histogram()+geom_vline(xintercept = accuracy)
    ggsave(perm_plot,filename = "permutaton_plot.png",device = "png")
    perm_p <- sum(perm_acc>accuracy)/length(perm_acc)
    cat(sprintf("\nPermutation p-value =  %1.3f\n",perm_p))
    perm_W <- perm_df %>% select(contains("W"))
    W_test <- sapply(perm_W, FUN = function(x) {
      abs(x)>abs(W)
    })
    W_sig <- rowMeans(W_test)
    b$perm_p <- perm_p
    b$perm_W_sig <- W_sig
    b$perm_W <- perm_W
    b$perm_results = perm_acc
    print(b$perm_p)
  }
  
  print(sprintf("Overall Accuracy: %1.3f; p = %.5f\n\n",sum(pred_data$model.pred==pred_data$drug)/dim(pred_data)[1],b$p.value))
  
  return(list(b,W,num_features,svm.model))
}

featureExtraction <- function(trainingData, feature_proportion = .1){
  # This is no longer used, but it will perform data driven feature selection.
  features <- as.matrix(trainingData %>% select(-drug))
  drug <- features[trainingData$drug==0,]
  placebo <- features[trainingData$drug==1,]
  
  tstats<-col_t_paired(placebo,drug)
  newdata=t(tstats[order(-abs(tstats$statistic)),])
  keep_features = colnames(newdata)[1:round(feature_proportion*dim(newdata)[2])]
  
  newTrainingData <- trainingData %>% select(drug,keep_features)
  return(list(newTrainingData,keep_features))
}

FD_thresh = .3
subData <- read.csv('/cbica/projects/alpraz_EI/input/alpraz_sublist_FD.csv')
subInfo <- subData %>% filter(exists==1 & motion_pass==1)
# Make sure whole subj is removed
# goodIDs <- subData %>%
#   na.omit()%>%
#   group_by(subid)%>%
#   summarise(motion_pass = max(FD,na.rm=T)<FD_thresh)%>%
#   filter(motion_pass==TRUE)%>%
#   select(subid)
# subInfo <- subData %>% filter(exists==1 & subid %in% goodIDs$subid)

PNCData <- read.csv('/cbica/projects/alpraz_EI/input/PNC_sublist_FD.csv')
PNCInfo <- PNCData %>% filter(exists==1 & FD < FD_thresh)

#### SET OPTIONS ####
classifier = "svm"
GSR="GSR"

# Permutation test?
# Do we want to use a permutation test for significance testing?
perm_test="permute_off"

# Do we want to pull out only certain brain areas?
# Set subdivide = TRUE if this is desired
# Set subdivision desired:
## "transmodal25" = top 25% most transmodal regions
## "unimodal25" = top 25% most unimodal regions
## "regional" = perform classification separately for each region in the atlas.
subdivide = TRUE
subdivision = "transmodal25"
cat(sprintf("\nsubdivision = %s\n",subdivision))

# Atlas and FE  
## Create a list of atlases that we want to use for classification. Now we are using only schaefer 400.
## "schaefer200x7_aal" can be used to do schaefer200 as well.
atlas_list = c("schaefer400x7_aal") 
fe_list=c(1)
###############

atlas_acc <- matrix(ncol = 6)
colnames(atlas_acc)=c("accuracy", "p.value", "fe","perm.p","atlas","num_features")

## Loop over atlases and run classification
for (atlasname in atlas_list){
  cat(atlasname)
  ## Load the data
  if (file.exists(sprintf("/cbica/projects/alpraz_EI/input/CorMats/%s_%s_%s.rdss",atlasname,GSR,subdivision))) {
    # The data has already been compiled into a dataframe, just load it.
    df <- readRDS(sprintf("/cbica/projects/alpraz_EI/input/CorMats/%s_%s_%s.rds",atlasname,GSR,subdivision))
  } else {
    # Compile the data if we haven't already.
    gm<- subInfo %>%
      group_by(subid, sesid) %>%
      mutate(mat = list(extract_matrix2(subid,sesid,atlasname,GSR,subdivide=subdivide,subdivision = subdivision)))
    
    if (subdivision=="regional") {
      mat <- do.call(rbind, lapply(gm$mat, function(x) data.frame(x)))
      df <- gm %>% full_join(mat,by=c("subid","sesid"))%>%select(-mat)
    } else {
      mat <- data.frame(matrix(unlist(gm$mat), nrow=dim(gm[1]), byrow=T))
      df <- cbind(gm %>%
                    ungroup() %>%
                    select(subid,sesid,drug)
                  ,mat)
    }
    
    print("saving...")
    saveRDS(df,file = sprintf("/cbica/projects/alpraz_EI/input/CorMats/%s_%s_%s.rds",atlasname,GSR,subdivision))
    print("saved")
    rm(gm)
  }
  
  ## Now run the SVM for a set of feature extraction levels or regions.
  acc = matrix(nrow = length(fe_list),ncol = 5)
  for (n in 1:length(fe_list)) {
    num <- fe_list[n]
    if (file.exists(
      sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rdss",
              atlasname,GSR,subdivision,classifier,as.character(fe_list[n]),perm_test
      ))
    ) {
      results <- readRDS(
        file = sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rds",
                       atlasname,GSR,subdivision,classifier,as.character(fe_list[n]),perm_test))
    } else{
      if (num == 1) {
        feature_selection = FALSE
      } else {
        feature_selection =TRUE
      }

      if (perm_test == "permute_on") {
        permutation_test = TRUE
      } else {
        permutation_test = FALSE
      }
      
      if (subdivision=="regional") {
        send_to_svm <- function(a_df){
          wide_df <- a_df %>%pivot_wider(id_cols = c(subid,sesid,drug),names_from = name,values_from = value)
          results <- run_model(df = wide_df,folds = 10,
                               feature_selection = feature_selection,feature_proportion = num,
                               permutation_test = permutation_test,num_permutations = 1000,
                               type = classifier)
        }
        nest_df <- df %>% group_by(row_label) %>% group_nest()
        r <- nest_df %>% mutate(results = purrr::map(data,send_to_svm))
        rdf <- r %>% 
          select(-data) %>% 
          mutate(acc = map(.x = results,.f = ~ unlist(.x,recursive = F)[["estimate"]]))%>%
          mutate(acc = unlist(acc))
        results <-rdf

      } else {
        results <- run_model(df = df,folds = 10,
                             feature_selection = feature_selection,feature_proportion = num,
                             permutation_test = permutation_test,num_permutations = 1000,
                             type = classifier)
        results[[5]]=atlasname
      }
      

      saveRDS(results,
              file = sprintf("/cbica/projects/alpraz_EI/output/drug_classification/%s_%s_%s_%s_%s_%s_results.rds",
                             atlasname,GSR,subdivision,classifier,as.character(fe_list[n]),perm_test))
    }

    # b <- results[[1]]
  }

# get critical values
# n=b$parameter #number of observations
# p01 <- qbinom(.01/2,size = n,prob = .5,lower.tail = F) +1 # get num correct that corresponds to two-tailed p-value of .01. Add one to get make the value represent prob of x >= p01 instead of x>p01.
# p01 = p01/n #convert to proprotion
# p001 <- (qbinom(.001/2,size = n,prob = .5,lower.tail = F) +1)/n


## Now we use the trained model on the PNC data.
## Load the PNC data
if (file.exists(sprintf("/cbica/projects/alpraz_EI/input/CorMats/PNC/%s_%s.rds",atlasname,subdivision))) {
  df <- readRDS(sprintf("/cbica/projects/alpraz_EI/input/CorMats/PNC/%s_%s.rds",atlasname,subdivision))
} else {
  gm<- PNCInfo %>% 
    group_by(subid,sesid) %>%
    mutate(mat = list(extract_matrix2(subid,sesid,atlasname,"PNC",subdivide=subdivide,subdivision = subdivision)))
  if (atlasname=="gordon333_aal") {
    feat_count=gm%>%summarise(n=length(unlist(mat)))
    good_bbl=feat_count%>%filter(n==19701)%>%select(subid)
    gm <- gm %>% filter(subid%in%good_bbl$subid)
  }
  mat <- data.frame(matrix(unlist(gm$mat), nrow=dim(gm[1]), byrow=T))
  df <- cbind(gm %>% 
                ungroup() %>% 
                select(subid,sesid)
              ,mat)
  print("saving...")
  saveRDS(df,file = sprintf("/cbica/projects/alpraz_EI/input/CorMats/PNC/%s_%s.rds",atlasname,subdivision))
  print("saved")
  rm(gm)
}
svm.model <- results[[4]]
testData <- as.matrix(df[, 3:dim(df)[2]])
svm.pred <- predict(svm.model, as.matrix(testData))

w <- t(svm.model$coefs) %*% svm.model$SV
decisionValues <- w %*% t(as.matrix(testData))-svm.model$rho
distance <- decisionValues/norm(w)
df$decisionValues <- t(decisionValues)
df$distance <- t(distance)
df$pred <-svm.pred
#save the output
saveRDS(df%>%select(subid,sesid,pred,distance,decisionValues),file = sprintf("/cbica/projects/alpraz_EI/output/PNC_predictions/%s_%s.rds",atlasname,subdivision))
}
