#Rscript.R for Linux, including KNN/SVM/NB

setwd("/public/home/06025/WORK/lianghui/data_classification")

#Data storage path and name
library(dplyr)
#library(dplyr, lib.loc="/public/home/06025/WORK/Software/mambaforge-pypy3/envs/lianghui/lib/R/library")
library(caret)
#library(caret, lib.loc="/public/home/06025/WORK/Software/mambaforge-pypy3/envs/lianghui/lib/R/library")
library(class)
#library(class, lib.loc="/public/home/06025/WORK/Software/mambaforge-pypy3/envs/lianghui/lib/R/library")
library(e1071)
#library(e1071, lib.loc="/public/home/06025/WORK/Software/mambaforge-pypy3/envs/lianghui/lib/R/library")
#library(mRMRe)
#library(CORElearn)

#prepared values需要准备的东西
ctrl <- trainControl(method="cv", number=10)
grid <- expand.grid(.laplace = seq(0, 1, 0.1), .usekernel = c(TRUE, FALSE), .adjust = 1)
trash<-c("FID","PAT","MAT","SEX","PHENOTYPE")
excludevars <- c("IID","ID","num","group")
#col_del<-c("IID","FID","num")
snp_num<-c(100,250,500,800,1000,2000,3000,5000,10000,15000)

#data file and processing读取文件
data<-read.table(file='Buffalo_187_HC_pass_simplify_num_qc_indep_fix_4k.raw',header=T,sep='')
group<-read.table(file="buffalo_3col.list",header=F,sep="")
names(group)<-c("IID","ID","group")
group$num<-1:nrow(group)
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

#数据文件处理
group$group<-as.factor(group$group)
data<-data[,!(names(data) %in% trash)]      #到这就只留IID和其他列了
names<-colnames(data)
names1<-gsub(pattern = "X(\\d+)\\.(\\d+)_.*", replacement = "\\1:\\2", names)
colnames(data)<-names1
#data<-data %>% select_if(~!any(is.na(.)))           #delete nas
#sum(is.na(data1))
fst<-read.table(file="Buffalo_187_HC_pass_simplify_num_qc_indep_fix_out.fst",header=T,sep="")
fst1 <- fst[apply(fst, 1, function(row) any(row %in% names1)), ]

#对于单个训练集，fst方法做特征选择和KNN、SVM、NB三种分类方法的结果
fst2<-fst1[order(fst1$FST,decreasing = T),]
extracted_data_list <- list()
for (len in snp_num){
  # 提取前len条数据
  extracted_data <- fst2[1:len, ]
  # 将提取的数据保存到列表中
  extracted_data_list[[length(extracted_data_list) + 1]] <- extracted_data
}

#======================================================
#第一回合兼流程测试
trainings<-trainingsets[[1]]
testings<-testingsets[[1]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)

acc_fst1_svm<-numeric(10)
acc_fst1_knn<-numeric(10)
acc_fst1_nb<-numeric(10)
acc_fst1_rf<-numeric(10)
fst_1list<-list()            #这个是生成保存的snplist的RData，方便后续我拉出来跑RF。这里也区分了12345次训练集。

for (n in 1:10){
    snp<-as.data.frame(extracted_data_list[[n]])         #从此处开始循环
    list_snp<-snp[,2]
    list_snp<-append("IID",list_snp)
    tryall<- data[,colnames(data) %in% list_snp]
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    #SVM
    set.seed(77)
    tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
    set.seed(77)
    svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
    pred_svm <- predict(svm,testdata_features)
    acc_fst1_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
    #KNN
    folds<-createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
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
    acc_fst1_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
    #NB
    tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    set.seed(77)
    model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
    set.seed(77)
    nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
    pred_nb<- predict(nb, testdata_features)
    acc_fst1_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
    #保存输出的SNP
    fst_1list[[paste0("snp_fst_", n)]] <- list_snp
}

#save保存SNP，懒得自己再跑了
acc_fst1_10<-data.frame(snp_num,acc_fst1_knn,acc_fst1_rf,acc_fst1_svm,acc_fst1_nb)
save(acc_fst1_10, file = "acc_fst1_10.RData")
save(fst_1list,file="snp_fst1_list.RData")
#最后可以在新的任务中使用“fst_1list$snp_fst_1”之类的代码提取我在之前跑出来的SNP


#======================================================
#第二回合
trainings<-trainingsets[[2]]
testings<-testingsets[[2]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)

acc_fst2_svm<-numeric(10)
acc_fst2_knn<-numeric(10)
acc_fst2_nb<-numeric(10)
acc_fst2_rf<-numeric(10)
fst_2list<-list()            #这个是生成保存的snplist的RData，方便后续我拉出来跑RF。这里也区分了12345次训练集。

for (n in 1:10){
    snp<-as.data.frame(extracted_data_list[[n]])         #从此处开始循环
    list_snp<-snp[,2]
    list_snp<-append("IID",list_snp)
    tryall<- data[,colnames(data) %in% list_snp]
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    #SVM
    set.seed(77)
    tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
    set.seed(77)
    svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
    pred_svm <- predict(svm,testdata_features)
    acc_fst2_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
    #KNN
    folds<-createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
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
    acc_fst2_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
    #NB
    tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    set.seed(77)
    model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
    set.seed(77)
    nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
    pred_nb<- predict(nb, testdata_features)
    acc_fst2_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
    #保存输出的SNP
    fst_2list[[paste0("snp_fst_", n)]] <- list_snp
}

#save保存SNP，懒得自己再跑了
acc_fst2_10<-data.frame(snp_num,acc_fst2_knn,acc_fst2_rf,acc_fst2_svm,acc_fst2_nb)
save(acc_fst2_10, file = "acc_fst2_10.RData")
save(fst_2list,file="snp_fst2_list.RData")
#最后可以在新的任务中使用“fst_2list$snp_fst_1”之类的代码提取我在之前跑出来的SNP


#======================================================
#第三回合
trainings<-trainingsets[[3]]
testings<-testingsets[[3]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)

acc_fst3_svm<-numeric(10)
acc_fst3_knn<-numeric(10)
acc_fst3_nb<-numeric(10)
acc_fst3_rf<-numeric(10)
fst_3list<-list()            #这个是生成保存的snplist的RData，方便后续我拉出来跑RF。这里也区分了12345次训练集。

for (n in 1:10){
    snp<-as.data.frame(extracted_data_list[[n]])         #从此处开始循环
    list_snp<-snp[,2]
    list_snp<-append("IID",list_snp)
    tryall<- data[,colnames(data) %in% list_snp]
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    #SVM
    set.seed(77)
    tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
    set.seed(77)
    svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
    pred_svm <- predict(svm,testdata_features)
    acc_fst3_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
    #KNN
    folds<-createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
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
    acc_fst3_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
    #NB
    tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    set.seed(77)
    model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
    set.seed(77)
    nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
    pred_nb<- predict(nb, testdata_features)
    acc_fst3_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
    #保存输出的SNP
    fst_3list[[paste0("snp_fst_", n)]] <- list_snp
}

#save保存SNP，懒得自己再跑了
acc_fst3_10<-data.frame(snp_num,acc_fst3_knn,acc_fst3_rf,acc_fst3_svm,acc_fst3_nb)
save(acc_fst3_10, file = "acc_fst3_10.RData")
save(fst_3list,file="snp_fst3_list.RData")
#最后可以在新的任务中使用“fst_3list$snp_fst_1”之类的代码提取我在之前跑出来的SNP


#======================================================
#第四回合
trainings<-trainingsets[[4]]
testings<-testingsets[[4]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)

acc_fst4_svm<-numeric(10)
acc_fst4_knn<-numeric(10)
acc_fst4_nb<-numeric(10)
acc_fst4_rf<-numeric(10)
fst_4list<-list()            #这个是生成保存的snplist的RData，方便后续我拉出来跑RF。这里也区分了12345次训练集。

for (n in 1:10){
    snp<-as.data.frame(extracted_data_list[[n]])         #从此处开始循环
    list_snp<-snp[,2]
    list_snp<-append("IID",list_snp)
    tryall<- data[,colnames(data) %in% list_snp]
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    #SVM
    set.seed(77)
    tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
    set.seed(77)
    svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
    pred_svm <- predict(svm,testdata_features)
    acc_fst4_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
    #KNN
    folds<-createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
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
    acc_fst4_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
    #NB
    tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    set.seed(77)
    model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
    set.seed(77)
    nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
    pred_nb<- predict(nb, testdata_features)
    acc_fst4_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
    #保存输出的SNP
    fst_4list[[paste0("snp_fst_", n)]] <- list_snp
}

#save保存SNP，懒得自己再跑了
acc_fst4_10<-data.frame(snp_num,acc_fst4_knn,acc_fst4_rf,acc_fst4_svm,acc_fst4_nb)
save(acc_fst4_10, file = "acc_fst4_10.RData")
save(fst_4list,file="snp_fst4_list.RData")
#最后可以在新的任务中使用“fst_4list$snp_fst_1”之类的代码提取我在之前跑出来的SNP


#======================================================
#第五回合
trainings<-trainingsets[[5]]
testings<-testingsets[[5]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)

acc_fst5_svm<-numeric(10)
acc_fst5_knn<-numeric(10)
acc_fst5_nb<-numeric(10)
acc_fst5_rf<-numeric(10)
fst_5list<-list()            #这个是生成保存的snplist的RData，方便后续我拉出来跑RF。这里也区分了12345次训练集。

for (n in 1:10){
    snp<-as.data.frame(extracted_data_list[[n]])         #从此处开始循环
    list_snp<-snp[,2]
    list_snp<-append("IID",list_snp)
    tryall<- data[,colnames(data) %in% list_snp]
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    #SVM
    set.seed(77)
    tunemodel<-tune.svm(x=traindata_features,y=traindata_labels,data=traindata,kernel="linear",gamma = 10^(-5:5),cost=2^(1:9))
    set.seed(77)
    svm<-svm(x=traindata_features,y=traindata_labels,kernel = "linear",gamma = tunemodel$best.parameters$gamma,cost=tunemodel$best.parameters$cost)
    pred_svm <- predict(svm,testdata_features)
    acc_fst5_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
    #KNN
    folds<-createFolds(traindata_labels, k = 15, list = TRUE, returnTrain = TRUE)
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
    acc_fst5_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
    #NB
    tryall <- tryall %>% mutate_all(~ if (is.numeric(.)) as.factor(.) else .)
    traindata<-merge(tryall,trainings,by="IID")
    testdata<-merge(tryall,testings,by="IID")
    traindata_features <- traindata[, !(names(traindata) %in% excludevars)]
    traindata_labels <- traindata[["group"]]
    testdata_features <- testdata[, !(names(testdata) %in% excludevars)]
    testdata_labels <- testdata[["group"]]
    set.seed(77)
    model <- train(x=traindata_features,y=traindata_labels, method="naive_bayes", trControl=ctrl, tuneGrid=grid)
    set.seed(77)
    nb <- e1071::naiveBayes(x=traindata_features,y=traindata_labels, laplace = model$bestTune$laplace)
    pred_nb<- predict(nb, testdata_features)
    acc_fst5_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
    #保存输出的SNP
    fst_5list[[paste0("snp_fst_", n)]] <- list_snp
}

#save保存SNP，懒得自己再跑了
acc_fst5_10<-data.frame(snp_num,acc_fst5_knn,acc_fst5_rf,acc_fst5_svm,acc_fst5_nb)
save(acc_fst5_10, file = "acc_fst5_10.RData")
save(fst_5list,file="snp_fst5_list.RData")
#最后可以在新的任务中使用“fst_5list$snp_fst_1”之类的代码提取我在之前跑出来的SNP

