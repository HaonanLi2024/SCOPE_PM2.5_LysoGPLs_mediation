setwd("D:\\Desktop\\李不怕呀\\Data & Writing\\LPA et al\\SCOPE\\20231124 update_文章1修改稿")
rawdata<-read.csv("rawdata_20231124.csv",header=TRUE,encoding="UTF-8")
str(rawdata)

rawdata<-rawdata[-which(rawdata$nvisit==1),]#只完成一次访视的删除



###### Figure 1_线性混合效应
library(lme4)
library(lmerTest)
library(splines)

pollutants<-"PM2.5_F"
biomarkers<-c("cPA.16.0","cPA.18.1",
              "LPA.16.0","LPA.18.0","LPA.18.1","LPA.18.2","LPA.20.3","LPA.20.4","LPA.22.6",
              "LysoPAF.C16","LysoPAF.C18","LysoPAF.C18.1",
              "LPG.16.0","LPG.16.1","LPG.18.0","LPG.18.1",
              "LysoPS.18.0","LysoPS.18.1",
              "IL_8_5p","MCP_1_5p","sCD40L_5p",
              "ifn_s",
              "plt.x","mpv.x","pdw.x","p_lcr","TxB2",
              "nonhdl")
resultofLMM<-data.frame(matrix(NA,length(biomarkers)*length(pollutants)*30,8))
names(resultofLMM)<-c("Biomarker","Pollutant","Time window (MA day)","per pollutant(ug/m3)","effect:mean","effect:lower confidence limit","effect:upper confidence limit","p-value")
a<-1#第几个生物指标
for(i in 1:length(biomarkers)){
  b<-1#第几种污染物
  for(j in 1:length(pollutants)){
    c<-1#第几个时间窗口
    for(k in 1:30){
      LMM<-lmer(log(get(biomarkers[a]),exp(1))~get(paste(pollutants[b],"_MA",c,"d",sep=""))+(1|id)+dow+ns(Temp_F_MA7d,1)+ns(RH_F_MA7d,1)+season+gender+age+bmi_mean+marry+edu_v2+income_v2+smoke+secdsmoke+fredrink,data=rawdata,na.action=na.omit)
      numofcolumn<-which(names(data.frame(summary(LMM)$coefficients))=="Estimate")
      numofcolumn2<-which(names(data.frame(summary(LMM)$coefficients))=="Std..Error")
      numofcolumn3<-which(names(data.frame(summary(LMM)$coefficients))=="Pr...t..")
      numofraw<-which(row.names(data.frame(summary(LMM)$coefficients))=="get(paste(pollutants[b], \"_MA\", c, \"d\", sep = \"\"))")
      beta<-data.frame(summary(LMM)$coefficients)[numofraw,numofcolumn]
      seofbeta<-data.frame(summary(LMM)$coefficients)[numofraw,numofcolumn2]
      p<-data.frame(summary(LMM)$coefficients)[numofraw,numofcolumn3]
      effect<-exp(beta*10)-1
      effect_l<-exp((beta-1.96*seofbeta)*10)-1
      effect_h<-exp((beta+1.96*seofbeta)*10)-1
      resultofLMM[(a-1)*length(pollutants)*30+(b-1)*30+c,]<-c(biomarkers[a],pollutants[b],c,"10",effect,effect_l,effect_h,p)
      c<-c+1
    }
    b<-b+1
  }
  a<-a+1
}

#为结果增加是否显著列
resultofLMM<-data.frame(resultofLMM)
resultofLMM$effect.mean<-as.numeric(resultofLMM$effect.mean)
resultofLMM$effect.lower.confidence.limit<-as.numeric(resultofLMM$effect.lower.confidence.limit)
resultofLMM$effect.upper.confidence.limit<-as.numeric(resultofLMM$effect.upper.confidence.limit)
resultofLMM$significant<-ifelse((resultofLMM$effect.lower.confidence.limit)*(resultofLMM$effect.upper.confidence.limit)<0,"Nonsignificant",ifelse(resultofLMM$effect.mean>0,"Significant Positive","Significant Negative"))
write.csv(resultofLMM,"result_PM2.5_F(MA) and all Biomarkers(one-pollutant)_1124.csv")



######Figure 2_中介效应分析
library(lme4)
detach("package:lmerTest")#卸载lmerTest包，否则mediate函数无法实现
library(splines)
library(mediation)

#部分时间窗口(outcome显著效应最高的时间窗口/效应最高的时间窗口），仅LPA
pollutants<-c("PM2.5_F")
biomarkers_mediators<-c("LPA.16.0","LPA.18.0","LPA.18.1","LPA.18.2","LPA.20.4","LPA.22.6",
                        "LysoPAF.C16","LysoPAF.C18","LysoPAF.C18.1",
                        "LPG.16.0","LPG.16.1","LPG.18.0",
                        "LPA.20.3","LysoPS.18.0")
biomarkers_outcomes<-c("IL_8_5p","MCP_1_5p","sCD40L_5p",
                       "ifn_s",
                       "plt.x","mpv.x","pdw.x","p_lcr","TxB2",
                       "nonhdl")
timewindow<-c(14,14,14,
              14,
              14,14,14,14,14,
              14)
resultofMediation<-data.frame(matrix(NA,length(biomarkers_outcomes)*length(biomarkers_mediators)*length(pollutants),20))
names(resultofMediation)<-c("Pollutant","Outcome","Mediator","Time","Mediation Proportion","Mediation Proportion 25%","Mediation Proportion 75%","Mediation Proportion p-Value"
                            ,"ACME (causal mediation effect)","ACME 25%","ACME 75%","ACME p-Value"
                            ,"ADE","ADE 25%","ADE 75%","ADE p-Value"
                            ,"Total Effect","Total Effect 25%","Total Effect 75%","Total Effect p-Value")
a<-1#第几个outcome
for(i in 1:length(biomarkers_outcomes)){
  b<-1#第几个mediator
  for(j in 1:length(biomarkers_mediators)){
    c<-1#第几个污染物
    for(k in 1:length(pollutants)){
      b_rawdata<-rawdata[-which(is.na(rawdata[,which(names(rawdata)==biomarkers_outcomes[a])])),]#剔除含有缺失值的数据使得med.fit和out.fit使用同样的观测（仅针对两模型不共同的变量即可）
      b_rawdata[,ncol(b_rawdata)+1]<-log(b_rawdata[which(names(b_rawdata)==biomarkers_outcomes[a])],exp(1))
      b_rawdata[,ncol(b_rawdata)+1]<-log(b_rawdata[which(names(b_rawdata)==biomarkers_mediators[b])],exp(1))#ln转化
      outcome<-paste("ln",biomarkers_outcomes[a],sep="")
      mediator<-paste("ln",biomarkers_mediators[b],sep="")
      pollutant<-paste(pollutants[c],"_MA",timewindow[a],"d",sep="")#给变量赋值
      names(b_rawdata)[(ncol(b_rawdata)-1):ncol(b_rawdata)]<-c(outcome,mediator)#给变量赋名
      med.fit<-lmer(get(mediator)~get(pollutant)+(1|id)+dow+ns(Temp_F_MA7d,1)+ns(RH_F_MA7d,1)+season+gender+age+bmi_mean+marry+edu_v2+income_v2+smoke+secdsmoke+fredrink,data=b_rawdata,na.action=na.omit)
      out.fit<-lmer(get(outcome)~get(mediator)+get(pollutant)+(1|id)+dow+ns(Temp_F_MA7d,1)+ns(RH_F_MA7d,1)+season+gender+age+bmi_mean+marry+edu_v2+income_v2+smoke+secdsmoke+fredrink,data=b_rawdata,na.action=na.omit)
      med.out<-mediate(med.fit,out.fit,treat='get(pollutant)',mediator='get(mediator)',treat.value=IQR(b_rawdata[,which(names(b_rawdata)==paste(pollutants[c],"_MA1d",sep=""))],na.rm=TRUE),control.value=0,sims = 1000)
      a1<-summary(med.out)$n.avg
      a2<-summary(med.out)$n.avg.ci[1]#置信区间25%四分位数
      a3<-summary(med.out)$n.avg.ci[2]#置信区间75%四分位数
      a4<-summary(med.out)$n.avg.p#Prop. Mediated
      b1<-summary(med.out)$d.avg
      b2<-summary(med.out)$d.avg.ci[1]
      b3<-summary(med.out)$d.avg.ci[2]
      b4<-summary(med.out)$d.avg.p#ACME (causal mediation effect)
      c1<-summary(med.out)$z.avg
      c2<-summary(med.out)$z.avg.ci[1]
      c3<-summary(med.out)$z.avg.ci[2]
      c4<-summary(med.out)$z.avg.p#ADE
      d1<-summary(med.out)$tau.coef
      d2<-summary(med.out)$tau.ci[1]
      d3<-summary(med.out)$tau.ci[2]
      d4<-summary(med.out)$tau.p#Total effect
      resultofMediation[(a-1)*length(biomarkers_mediators)*length(pollutants)+(b-1)*length(pollutants)+c,]<-c(pollutants[c],biomarkers_outcomes[a],biomarkers_mediators[b],paste("MA",timewindow[a],"d",sep="")
                                                                                                              ,a1,a2,a3,a4
                                                                                                              ,b1,b2,b3,b4
                                                                                                              ,c1,c2,c3,c4
                                                                                                              ,d1,d2,d3,d4)
      c<-c+1
    }
    b<-b+1
  }
  a<-a+1
}

#FDR校正
resultofMediation<-resultofMediation[order(resultofMediation$"Mediation Proportion p-Value"),]
resultofMediation$BH_proportion<-p.adjust(resultofMediation$"Mediation Proportion p-Value",method="BH")
write.csv(resultofMediation,"result_mediation(PM2.5_F to LPLs to pro- biomarkers_14d)_20231124.csv")
