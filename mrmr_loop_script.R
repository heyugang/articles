#KNN/SVM/NB
setwd("/public/home/06025/WORK/lianghui/data_classification")
#Data storage path and name
library(dplyr)
library(caret)
library(class)
library(e1071)
library(mRMRe)
#library(CORElearn)

#准备一些可以用到的向量
trash<-c("FID","PAT","MAT","SEX","PHENOTYPE")
non_numeric_columns<-c("FID","IID","PAT","MAT","SEX","PHENOTYPE")
snp_num<-c(100,250,500,800,1000,2000,3000,5000,10000)             #这个重要，要改成自己需要的组别，现在已经是全部的
excludevars <- c("ID","IID","num","group")
excludevars_all<-c("ID","IID","num","group","grouplevels")
ctrl <- trainControl(method="cv", number=10)
grid <- expand.grid(.laplace = seq(0, 1, 0.1), .usekernel = c(TRUE, FALSE), .adjust = 1)

#读取文件read files
data<-read.table(file='Buffalo_187_HC_pass_simplify_num_qc_indep_fix_4k.raw',header=T,sep='')
group<-read.table(file="buffalo_3col.list",header=F,sep="")
names(group)<-c("IID","ID","group")
group$num<-1:nrow(group)
group$group<-as.factor(group$group)               #转换列格式
load("trainset_list.RData")        #读取已经抽取出来的list列表文件，里面是随机数抽取的训练集数字，后续会对应上的
trainingsets<-list()
testingsets<-list()
for (i in 1:5) {
    # 提取 training set
    train <- number_list[[i]]
    trainingset <- group[group$num %in% train, ]
    trainingset <- trainingset[order(trainingset$num, decreasing = FALSE), ]
    # 提取 testing set
    testingset <- subset(group, !(num %in% trainingset$num))
    # 存储训练集和测试集到列表中
    trainingsets[[i]] <- trainingset
    testingsets[[i]] <- testingset
}

#文件处理
data1<-data
numeric_columns_tmp <- setdiff(names(data), non_numeric_columns)
## 将这些数字列的类型从int转换为numeric
data1[numeric_columns_tmp] <- lapply(data1[numeric_columns_tmp], function(x) as.numeric(as.character(x)))
data1<-data1[,!(names(data1) %in% trash)]      #到这就只留IID和其他列了
names<-colnames(data1)
names1<-gsub(pattern = "X(\\d+)\\.(\\d+)_.*", replacement = "\\1:\\2", names)
colnames(data1)<-names1

#======================================
###第一回合
#准备储存分类准确率的向量（*要根据前面的snp_num改成对应的数字）（可以根据自己需要改名）
acc_mrmr1_svm<-numeric(9)
acc_mrmr1_knn<-numeric(9)
acc_mrmr1_nb<-numeric(9)
acc_mrmr1_rf<-numeric(9)
mrmr_all_1list<-list()            #跑出来之后可以通过mrmr_all_1list$这个命令提取对应的SNP

#循环进行
trainings<-trainingsets[[1]]
testings<-testingsets[[1]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)
trainings$grouplevels<-as.factor(trainings$group)
levels(trainings$grouplevels)<-1:nlevels(trainings$grouplevels)
trainings$grouplevels<-as.ordered(trainings$grouplevels)
try<-merge(data1,trainings,by="IID")
try<-try[,!(names(try) %in% excludevars)]
try <- try %>% mutate_at(vars(-ncol(try)),~as.factor(as.numeric(.) + 1))
try <- try %>% mutate_at(vars(-ncol(try)),~ if (is.factor(.)) ordered(.) else .)
mrmr_data <- mRMR.data(data = try)
for (n in 1:9) {
  mrmr.solution <- mRMR.classic(data = mrmr_data, target_indices = 39997, feature_count = snp_num[[n]])
  selected_features <- solutions(mrmr.solution)[[1]]
  selected_try <- try[, selected_features]
  snp<-colnames(selected_try)
  snp<-append("IID",snp)
  tryall<-data1[,colnames(data1) %in% snp]
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  #SVM
  set.seed(77)
  tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
  set.seed(77)
  svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
  pred_svm <- predict(svm,testdata_features)
  acc_mrmr1_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
  #KNN
  folds <- createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
  best_k <- 1
  best_accuracy <- 0
  for (k in 1:15) {
    fold_accuracy <- 0
    for (i in 1:15) {
      train_indices <- folds[[i]]
      train_data <- traindata_features[train_indices, ]
      train_labels <- traindata_labels[train_indices]
      test_data <- traindata_features[-train_indices, ]
      test_labels <- traindata_labels[-train_indices]
      predictions <- class::knn(train = train_data, test = test_data, cl = train_labels, k = k)
      fold_accuracy <- fold_accuracy + sum(predictions == test_labels) / length(test_labels)
    }
    fold_accuracy <- fold_accuracy / 10
    if (fold_accuracy > best_accuracy) {
      best_k <- k
      best_accuracy <- fold_accuracy
    }
  }
  k_value <- best_k
  set.seed(77)
  pred_knn <- class::knn(traindata_features, testdata_features, traindata_labels, k=k_value)
  acc_mrmr1_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
  #NB
  tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  set.seed(77)
  model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
  set.seed(77)
  nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
  pred_nb<- predict(nb, testdata_features)
  acc_mrmr1_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
  #保存到list文件
  mrmr_all_1list[[paste0("snp_mrmr_", n)]] <- snp
  #跑完了流程就删掉，节省空间
  rm(mrmr.solution)
}

#save保存SNP，就不用自己再跑了
acc_mrmr1_9<-data.frame(snp_num,acc_mrmr1_knn,acc_mrmr1_rf,acc_mrmr1_svm,acc_mrmr1_nb)
save(acc_mrmr1_9, file = "acc_mrmr1_9.RData")
save(mrmr_all_1list,file="snp_mrmr_all_1list.RData")
#最后可以在新的任务中使用“mrmr_1list$snp_mrmr_1”之类的代码提取我在之前跑出来的SNP


#======================================
###第二回合
#准备储存分类准确率的向量（*要根据前面的snp_num改成对应的数字）（可以根据自己需要改名）
acc_mrmr2_svm<-numeric(9)
acc_mrmr2_knn<-numeric(9)
acc_mrmr2_nb<-numeric(9)
acc_mrmr2_rf<-numeric(9)
mrmr_all_2list<-list()            #跑出来之后可以通过mrmr_all_1list$这个命令提取对应的SNP

#循环进行
trainings<-trainingsets[[2]]
testings<-testingsets[[2]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)
trainings$grouplevels<-as.factor(trainings$group)
levels(trainings$grouplevels)<-1:nlevels(trainings$grouplevels)
trainings$grouplevels<-as.ordered(trainings$grouplevels)
try<-merge(data1,trainings,by="IID")
try<-try[,!(names(try) %in% excludevars)]
try <- try %>% mutate_at(vars(-ncol(try)),~as.factor(as.numeric(.) + 1))
try <- try %>% mutate_at(vars(-ncol(try)),~ if (is.factor(.)) ordered(.) else .)
mrmr_data <- mRMR.data(data = try)
for (n in 1:9) {
  mrmr.solution <- mRMR.classic(data = mrmr_data, target_indices = 39997, feature_count = snp_num[[n]])
  selected_features <- solutions(mrmr.solution)[[1]]
  selected_try <- try[, selected_features]
  snp<-colnames(selected_try)
  snp<-append("IID",snp)
  tryall<-data1[,colnames(data1) %in% snp]
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  #SVM
  set.seed(77)
  tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
  set.seed(77)
  svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
  pred_svm <- predict(svm,testdata_features)
  acc_mrmr2_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
  #KNN
  folds <- createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
  best_k <- 1
  best_accuracy <- 0
  for (k in 1:15) {
    fold_accuracy <- 0
    for (i in 1:15) {
      train_indices <- folds[[i]]
      train_data <- traindata_features[train_indices, ]
      train_labels <- traindata_labels[train_indices]
      test_data <- traindata_features[-train_indices, ]
      test_labels <- traindata_labels[-train_indices]
      predictions <- class::knn(train = train_data, test = test_data, cl = train_labels, k = k)
      fold_accuracy <- fold_accuracy + sum(predictions == test_labels) / length(test_labels)
    }
    fold_accuracy <- fold_accuracy / 10
    if (fold_accuracy > best_accuracy) {
      best_k <- k
      best_accuracy <- fold_accuracy
    }
  }
  k_value <- best_k
  set.seed(77)
  pred_knn <- class::knn(traindata_features, testdata_features, traindata_labels, k=k_value)
  acc_mrmr2_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
  #NB
  tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  set.seed(77)
  model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
  set.seed(77)
  nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
  pred_nb<- predict(nb, testdata_features)
  acc_mrmr2_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
  #保存到list文件
  mrmr_all_2list[[paste0("snp_mrmr_", n)]] <- snp
  #跑完了流程就删掉，节省空间
  rm(mrmr.solution)
}

#save保存SNP，就不用自己再跑了
acc_mrmr2_9<-data.frame(snp_num,acc_mrmr2_knn,acc_mrmr2_rf,acc_mrmr2_svm,acc_mrmr2_nb)
save(acc_mrmr2_9, file = "acc_mrmr2_9.RData")
save(mrmr_all_2list,file="snp_mrmr_all_2list.RData")
#最后可以在新的任务中使用“mrmr_1list$snp_mrmr_1”之类的代码提取我在之前跑出来的SNP


#======================================
###第三回合
#准备储存分类准确率的向量（*要根据前面的snp_num改成对应的数字）（可以根据自己需要改名）
acc_mrmr3_svm<-numeric(9)
acc_mrmr3_knn<-numeric(9)
acc_mrmr3_nb<-numeric(9)
acc_mrmr3_rf<-numeric(9)
mrmr_all_3list<-list()            #跑出来之后可以通过mrmr_all_1list$这个命令提取对应的SNP

#循环进行
trainings<-trainingsets[[3]]
testings<-testingsets[[3]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)
trainings$grouplevels<-as.factor(trainings$group)
levels(trainings$grouplevels)<-1:nlevels(trainings$grouplevels)
trainings$grouplevels<-as.ordered(trainings$grouplevels)
try<-merge(data1,trainings,by="IID")
try<-try[,!(names(try) %in% excludevars)]
try <- try %>% mutate_at(vars(-ncol(try)),~as.factor(as.numeric(.) + 1))
try <- try %>% mutate_at(vars(-ncol(try)),~ if (is.factor(.)) ordered(.) else .)
mrmr_data <- mRMR.data(data = try)
for (n in 1:9) {
  mrmr.solution <- mRMR.classic(data = mrmr_data, target_indices = 39997, feature_count = snp_num[[n]])
  selected_features <- solutions(mrmr.solution)[[1]]
  selected_try <- try[, selected_features]
  snp<-colnames(selected_try)
  snp<-append("IID",snp)
  tryall<-data1[,colnames(data1) %in% snp]
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  #SVM
  set.seed(77)
  tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
  set.seed(77)
  svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
  pred_svm <- predict(svm,testdata_features)
  acc_mrmr3_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
  #KNN
  folds <- createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
  best_k <- 1
  best_accuracy <- 0
  for (k in 1:15) {
    fold_accuracy <- 0
    for (i in 1:15) {
      train_indices <- folds[[i]]
      train_data <- traindata_features[train_indices, ]
      train_labels <- traindata_labels[train_indices]
      test_data <- traindata_features[-train_indices, ]
      test_labels <- traindata_labels[-train_indices]
      predictions <- class::knn(train = train_data, test = test_data, cl = train_labels, k = k)
      fold_accuracy <- fold_accuracy + sum(predictions == test_labels) / length(test_labels)
    }
    fold_accuracy <- fold_accuracy / 10
    if (fold_accuracy > best_accuracy) {
      best_k <- k
      best_accuracy <- fold_accuracy
    }
  }
  k_value <- best_k
  set.seed(77)
  pred_knn <- class::knn(traindata_features, testdata_features, traindata_labels, k=k_value)
  acc_mrmr3_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
  #NB
  tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  set.seed(77)
  model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
  set.seed(77)
  nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
  pred_nb<- predict(nb, testdata_features)
  acc_mrmr3_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
  #保存到list文件
  mrmr_all_3list[[paste0("snp_mrmr_", n)]] <- snp
  #跑完了流程就删掉，节省空间
  rm(mrmr.solution)
}

#save保存SNP，就不用自己再跑了
acc_mrmr3_9<-data.frame(snp_num,acc_mrmr3_knn,acc_mrmr3_rf,acc_mrmr3_svm,acc_mrmr3_nb)
save(acc_mrmr3_9, file = "acc_mrmr3_9.RData")
save(mrmr_all_3list,file="snp_mrmr_all_3list.RData")
#最后可以在新的任务中使用“mrmr_1list$snp_mrmr_1”之类的代码提取我在之前跑出来的SNP


#======================================
###第四回合
#准备储存分类准确率的向量（*要根据前面的snp_num改成对应的数字）（可以根据自己需要改名）
acc_mrmr4_svm<-numeric(9)
acc_mrmr4_knn<-numeric(9)
acc_mrmr4_nb<-numeric(9)
acc_mrmr4_rf<-numeric(9)
mrmr_all_4list<-list()            #跑出来之后可以通过mrmr_all_1list$这个命令提取对应的SNP

#循环进行
trainings<-trainingsets[[4]]
testings<-testingsets[[4]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)
trainings$grouplevels<-as.factor(trainings$group)
levels(trainings$grouplevels)<-1:nlevels(trainings$grouplevels)
trainings$grouplevels<-as.ordered(trainings$grouplevels)
try<-merge(data1,trainings,by="IID")
try<-try[,!(names(try) %in% excludevars)]
try <- try %>% mutate_at(vars(-ncol(try)),~as.factor(as.numeric(.) + 1))
try <- try %>% mutate_at(vars(-ncol(try)),~ if (is.factor(.)) ordered(.) else .)
mrmr_data <- mRMR.data(data = try)
for (n in 1:9) {
  mrmr.solution <- mRMR.classic(data = mrmr_data, target_indices = 39997, feature_count = snp_num[[n]])
  selected_features <- solutions(mrmr.solution)[[1]]
  selected_try <- try[, selected_features]
  snp<-colnames(selected_try)
  snp<-append("IID",snp)
  tryall<-data1[,colnames(data1) %in% snp]
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  #SVM
  set.seed(77)
  tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
  set.seed(77)
  svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
  pred_svm <- predict(svm,testdata_features)
  acc_mrmr4_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
  #KNN
  folds <- createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
  best_k <- 1
  best_accuracy <- 0
  for (k in 1:15) {
    fold_accuracy <- 0
    for (i in 1:15) {
      train_indices <- folds[[i]]
      train_data <- traindata_features[train_indices, ]
      train_labels <- traindata_labels[train_indices]
      test_data <- traindata_features[-train_indices, ]
      test_labels <- traindata_labels[-train_indices]
      predictions <- class::knn(train = train_data, test = test_data, cl = train_labels, k = k)
      fold_accuracy <- fold_accuracy + sum(predictions == test_labels) / length(test_labels)
    }
    fold_accuracy <- fold_accuracy / 10
    if (fold_accuracy > best_accuracy) {
      best_k <- k
      best_accuracy <- fold_accuracy
    }
  }
  k_value <- best_k
  set.seed(77)
  pred_knn <- class::knn(traindata_features, testdata_features, traindata_labels, k=k_value)
  acc_mrmr4_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
  #NB
  tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  set.seed(77)
  model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
  set.seed(77)
  nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
  pred_nb<- predict(nb, testdata_features)
  acc_mrmr4_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
  #保存到list文件
  mrmr_all_4list[[paste0("snp_mrmr_", n)]] <- snp
  #跑完了流程就删掉，节省空间
  rm(mrmr.solution)
}

#save保存SNP，就不用自己再跑了
acc_mrmr4_9<-data.frame(snp_num,acc_mrmr4_knn,acc_mrmr4_rf,acc_mrmr4_svm,acc_mrmr4_nb)
save(acc_mrmr4_9, file = "acc_mrmr4_9.RData")
save(mrmr_all_4list,file="snp_mrmr_all_4list.RData")
#最后可以在新的任务中使用“mrmr_1list$snp_mrmr_1”之类的代码提取我在之前跑出来的SNP


#======================================
###第五回合
#准备储存分类准确率的向量（*要根据前面的snp_num改成对应的数字）（可以根据自己需要改名）
acc_mrmr5_svm<-numeric(9)
acc_mrmr5_knn<-numeric(9)
acc_mrmr5_nb<-numeric(9)
acc_mrmr5_rf<-numeric(9)
mrmr_all_5list<-list()            #跑出来之后可以通过mrmr_all_1list$这个命令提取对应的SNP

#循环进行
trainings<-trainingsets[[5]]
testings<-testingsets[[5]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)
trainings$grouplevels<-as.factor(trainings$group)
levels(trainings$grouplevels)<-1:nlevels(trainings$grouplevels)
trainings$grouplevels<-as.ordered(trainings$grouplevels)
try<-merge(data1,trainings,by="IID")
try<-try[,!(names(try) %in% excludevars)]
try <- try %>% mutate_at(vars(-ncol(try)),~as.factor(as.numeric(.) + 1))
try <- try %>% mutate_at(vars(-ncol(try)),~ if (is.factor(.)) ordered(.) else .)
mrmr_data <- mRMR.data(data = try)
for (n in 1:9) {
  mrmr.solution <- mRMR.classic(data = mrmr_data, target_indices = 39997, feature_count = snp_num[[n]])
  selected_features <- solutions(mrmr.solution)[[1]]
  selected_try <- try[, selected_features]
  snp<-colnames(selected_try)
  snp<-append("IID",snp)
  tryall<-data1[,colnames(data1) %in% snp]
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  #SVM
  set.seed(77)
  tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
  set.seed(77)
  svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
  pred_svm <- predict(svm,testdata_features)
  acc_mrmr5_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
  #KNN
  folds <- createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
  best_k <- 1
  best_accuracy <- 0
  for (k in 1:15) {
    fold_accuracy <- 0
    for (i in 1:15) {
      train_indices <- folds[[i]]
      train_data <- traindata_features[train_indices, ]
      train_labels <- traindata_labels[train_indices]
      test_data <- traindata_features[-train_indices, ]
      test_labels <- traindata_labels[-train_indices]
      predictions <- class::knn(train = train_data, test = test_data, cl = train_labels, k = k)
      fold_accuracy <- fold_accuracy + sum(predictions == test_labels) / length(test_labels)
    }
    fold_accuracy <- fold_accuracy / 10
    if (fold_accuracy > best_accuracy) {
      best_k <- k
      best_accuracy <- fold_accuracy
    }
  }
  k_value <- best_k
  set.seed(77)
  pred_knn <- class::knn(traindata_features, testdata_features, traindata_labels, k=k_value)
  acc_mrmr5_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
  #NB
  tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
  traindata<-merge(tryall,trainings,by="IID")
  testdata<-merge(tryall,testings,by="IID")
  traindata_features <- traindata[, !(names(traindata) %in% excludevars_all)]
  traindata_labels <- traindata[["group"]]
  testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
  testdata_labels <- testdata[["group"]]
  set.seed(77)
  model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
  set.seed(77)
  nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
  pred_nb<- predict(nb, testdata_features)
  acc_mrmr5_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
  #保存到list文件
  mrmr_all_5list[[paste0("snp_mrmr_", n)]] <- snp
  #跑完了流程就删掉，节省空间
  rm(mrmr.solution)
}

#save保存SNP，就不用自己再跑了
acc_mrmr5_9<-data.frame(snp_num,acc_mrmr5_knn,acc_mrmr5_rf,acc_mrmr5_svm,acc_mrmr5_nb)
save(acc_mrmr5_9, file = "acc_mrmr5_9.RData")
save(mrmr_all_5list,file="snp_mrmr_all_5list.RData")
#最后可以在新的任务中使用“mrmr_1list$snp_mrmr_1”之类的代码提取我在之前跑出来的SNP

