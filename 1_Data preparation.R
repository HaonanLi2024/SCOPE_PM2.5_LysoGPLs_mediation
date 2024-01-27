setwd("C:\\Users\\admin\\Desktop\\李不怕呀\\Data & Writing\\LPA et al\\SCOPE\\20231124 update_文章1修改稿")

rawdata_a<-read.csv("2. scope  imp and health - 李浩楠.csv",header=TRUE,as.is=c("date"))#访视数据
str(rawdata_a)#检查变量数据类型（especially日期），核对observation、variable的数量，缺失值情况!
rawdata_a$date<-as.Date(rawdata_a$date,"%m/%d/%Y")#"input format"

rawdata_b<-read.csv("3. LPA et al.csv",header=TRUE)#代谢物数据
str(rawdata_b)
rawdata_b$date<-as.Date(rawdata_b$date,"%Y-%m-%d")

rawdata_c<-read.csv("1. air - 李浩楠.csv",header=TRUE)#污染物数据
str(rawdata_c)



#根据受试者id和访视date合并a_visit raw data和b_LPA et al raw data,未匹配的也保留
merge_a<-merge(rawdata_a, rawdata_b,by=c("id","date"),all=TRUE)



#匹配Lag的污染物数据
rownames(rawdata_c)[which(rownames(rawdata_c)=="date")]<-"before_xd"
rawdata_c$before_xd<-as.Date(rawdata_c$before_xd,"%m/%d/%Y")

library(lubridate)
merge_b<-merge_a
merge_b$before_xd<-merge_b$date
a<-1
rownumber<-ncol(merge_b)+1
for(i in 1:30){#计算访视前1-30天
  merge_b$before_xd<-merge_b$before_xd-rep(days(1),times=length(merge_b$before_xd));
  merge_b<-merge(merge_b, rawdata_c,by="before_xd",all.x=TRUE);
  names(merge_b)[rownumber:(rownumber+ncol(rawdata_c)-2)]<-paste(names(merge_b)[rownumber:(rownumber+ncol(rawdata_c)-2)],"_Lag",a,"d",sep="");
  a<-a+1;
  rownumber<-rownumber+ncol(rawdata_c)-1;
}



#匹配Moving Average的污染物数据
merge_c<-merge_b#merge_c初始从ncol(merge_a)+2到ncol(merge_b)是污染物Lag数据

proportionofnoneNA<-function(x){#自定义函数计算NA比例
  z<-sum(!is.na(x))/length(x)
  return(z)
}

b<-1#b代表第几种污染物
for(j in 1:(ncol(rawdata_c)-1)){
  pollutantsX<-data.frame(matrix(NA,nrow(merge_c),1))
  c<-1#c代表average前几天
  for(k in 1:30){
    noneNApropotion<-apply(merge_c[seq(ncol(merge_a)+1+b,ncol(merge_a)+1+b+(ncol(rawdata_c)-1)*(c-1),ncol(rawdata_c)-1)],1,proportionofnoneNA)
    meanpollutantsX<-apply(merge_c[seq(ncol(merge_a)+1+b,ncol(merge_a)+1+b+(ncol(rawdata_c)-1)*(c-1),ncol(rawdata_c)-1)],1,mean,na.rm=TRUE)
    #对于有效数据＜2/3的，赋值为NA
    d<-1
    for(m in 1:nrow(merge_c)){
      if(noneNApropotion[d]<2/3)meanpollutantsX[d]<-NA
      d<-d+1
    }
    pollutantsX<-cbind(pollutantsX,meanpollutantsX)
    c<-c+1
  }
  pollutantsX<-pollutantsX[,-1]
  names(pollutantsX)[1:30]<-paste(names(rawdata_c)[b+1],"_MA",1:30,"d",sep="")
  merge_c<-cbind(merge_c,pollutantsX)
  b<-b+1
}



write.csv(merge_c,"rawdata_20231124.csv")
