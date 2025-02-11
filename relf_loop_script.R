#Rscript.R for Linux, including KNN/SVM/NB

setwd("/public/home/06025/WORK/lianghui/data_classification")

#Data storage path and name
library(dplyr)
library(caret)
library(class)
library(e1071)
#library(mRMRe)
library(CORElearn)


#prepared values需要准备的东西
ctrl <- trainControl(method="cv", number=10)
grid <- expand.grid(.laplace = seq(0, 1, 0.1), .usekernel = c(TRUE, FALSE), .adjust = 1)
trash<-c("FID","PAT","MAT","SEX","PHENOTYPE")
excludevars <- c("ID","IID","num","group")
col_del<-c("IID","ID","num")
snp_num<-c(100,250,500,800,1000,2000,3000,5000,10000,15000)

#data file and processing读取文件
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


#数据处理
data<-data[,!(names(data) %in% trash)]      #到这就只留IID和其他列了
names<-colnames(data)
names1<-gsub(pattern = "X(\\d+)\\.(\\d+)_.*", replacement = "\\1:\\2", names)
colnames(data)<-names1
#data<-data %>% select_if(~!any(is.na(.)))           #delete nas
#sum(is.na(data1))


#======================================
###第一回合
#设置对应的list向量
acc_relf1_svm<-numeric(10)
acc_relf1_knn<-numeric(10)
acc_relf1_nb<-numeric(10)
acc_relf1_rf<-numeric(10)
relf_1list<-list()

#提训练集
trainings<-trainingsets[[1]]
testings<-testingsets[[1]]
trainings$group<-as.factor(trainings$group)
testings$group<-as.factor(testings$group)

#对于单个训练集，relief-F方法做特征选择和KNN、SVM、NB三种分类方法，循环进行
try<-merge(data,trainings,by="IID")
try<-try[,!(names(try) %in% col_del)]
result <- attrEval(group ~ ., data = try, estimator="ReliefFbestK")
for (n in 1:10){
    N<-snp_num[[n]]
    top_features <- names(sort(result, decreasing = TRUE)[1:N])
    top_features<-append("IID",top_features)
    tryall <- data[, top_features]
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
    acc_relf1_svm[n]<-sum(pred_svm== testdata_labels) / length(testdata_labels)
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
    acc_relf1_knn[n]<-sum(pred_knn== testdata_labels) / length(testdata_labels)
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
    acc_relf1_nb[n]<-sum(pred_nb== testdata_labels) / length(testdata_labels)
    #保存输出的SNP
    relf_1list[[paste0("snp_relf_", n)]] <- top_features
}

#save保存SNP，懒得自己再跑了
acc_relf1_10<-data.frame(snp_num,acc_relf1_knn,acc_relf1_rf,acc_relf1_svm,acc_relf1_nb)
save(acc_relf1_10, file = "acc_relf1_10.RData")
save(relf_1list,file="snp_relf_1list.RData")
#最后可以在新的任务中使用“relf_1list$snp_relf_1”之类的代码提取我在之前跑出来的SNP


