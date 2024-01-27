#设置工作空间，正斜杠please
setwd("D:\\Desktop\\李不怕呀\\Data & Writing\\LPA et al\\SCOPE\\20231124 update_文章1修改稿")
#读入CSV数据文件），首行是标题
rawdata<-read.csv("rawdata_20231124.csv",header=TRUE,encoding="UTF-8")
#检查变量数据类型（especially日期，字符型向量character和因子factor的区别），核对observation、variable的数量，缺失值情况
rawdata$date<-as.Date(rawdata$date,"%Y-%m-%d")
str(rawdata)
#只完成一次访视的删除
rawdata<-rawdata[-which(rawdata$nvisit==1),]


###correlation
biomarkers<-c("cPA.16.0","cPA.18.1",
              "LPA.16.0","LPA.18.0","LPA.18.1","LPA.18.2","LPA.20.3","LPA.20.4","LPA.22.6",
              "LysoPAF.C16","LysoPAF.C18","LysoPAF.C18.1",
              "LPG.16.0","LPG.16.1","LPG.18.0","LPG.18.1",
              "LysoPS.18.0","LysoPS.18.1",
              "IL_8_5p","MCP_1_5p","sCD40L_5p","ifn_s","nonhdl")
MET<-data.frame(matrix(NA,nrow(rawdata),1))
MET<-MET[,-1]

a<-1#第几个溶血甘油磷脂
for(i in 1:length(biomarkers)){
  METx<-rawdata[which(names(rawdata)==biomarkers[a])]
  MET<-cbind(MET,METx)
  a<-a+1
}
b<-cor(MET,method='spearman',use='pairwise.complete.obs')
library(corrplot)
c<-cor.mtest(MET,method='spearman',use='pairwise.complete.obs')$p
COLx<-colorRampPalette(c("#02546F","white","red"))
library(Cairo)
CairoPNG(file="FigureS2_20231130.png",
         width=10.5,
         height=10.5,
         units="in",
         dpi=500)
corrplot(b,type="upper",p.mat=c,sig.level=0.05,insig="blank",col=COLx(100),tl.col="black",tl.pos="tl")
corrplot(b,add=TRUE,type="lower",method="number",p.mat=c,sig.level=0.05,insig="blank",col=COLx(100),tl.pos="n",cl.pos="n")
dev.off()



##########敏感性分析
##只在距离监测站3km内受试者中

#Figure S3_LME(3km)
#设置工作空间，正斜杠please
setwd("D:/Desktop/李不怕呀/Data & Writing/LPA et al/SCOPE/20231124 update_文章1修改稿")
#读入CSV数据文件），首行是标题
rawdata<-read.csv("rawdata_20231124.csv",header=TRUE,encoding="UTF-8")
#检查变量数据类型（especially日期，字符型向量character和因子factor的区别），核对observation、variable的数量，缺失值情况
rawdata$date<-as.Date(rawdata$date,"%Y-%m-%d")
str(rawdata)
#只完成一次访视的删除
rawdata<-rawdata[-which(rawdata$nvisit==1),]
rawdata<-rawdata[which(rawdata$distance<=3),]

pollutants<-"PM2.5_F"
biomarkers<-c("cPA.16.0","cPA.18.1",
              "LPA.16.0","LPA.18.0","LPA.18.1","LPA.18.2","LPA.20.3","LPA.20.4","LPA.22.6",
              "LysoPAF.C16","LysoPAF.C18","LysoPAF.C18.1",
              "LPG.16.0","LPG.16.1","LPG.18.0","LPG.18.1",
              "LysoPS.18.0","LysoPS.18.1",
              "IL_8_5p","MCP_1_5p","sCD40L_5p",
              "ifn_s",
              "nonhdl")
timewindow<-c(14,30)
resultofLMM<-data.frame(matrix(NA,length(biomarkers)*length(pollutants),8))
names(resultofLMM)<-c("Biomarker","Pollutant","Time window (MA day)","per pollutant(ug/m3)","effect:mean","effect:lower confidence limit","effect:upper confidence limit","p-value")
a<-1#第几个生物指标
for(i in 1:length(biomarkers)){
  b<-1#第几种污染物
  for(j in 1:length(pollutants)){
    c<-1#第几个时间窗口
    for(k in 1:length(timewindow)){
      LMM<-lmer(log(get(biomarkers[a]),exp(1))~get(paste(pollutants[b],"_MA",timewindow[c],"d",sep=""))+(1|id)+dow+ns(Temp_F_MA7d,1)+ns(RH_F_MA7d,1)+season+gender+age+bmi_mean+marry+edu_v2+income_v2+smoke+secdsmoke+fredrink,data=rawdata,na.action=na.omit)
      numofcolumn<-which(names(data.frame(summary(LMM)$coefficients))=="Estimate")
      numofcolumn2<-which(names(data.frame(summary(LMM)$coefficients))=="Std..Error")
      numofcolumn3<-which(names(data.frame(summary(LMM)$coefficients))=="Pr...t..")
      numofraw<-which(row.names(data.frame(summary(LMM)$coefficients))=="get(paste(pollutants[b], \"_MA\", timewindow[c], \"d\", sep = \"\"))")
      beta<-data.frame(summary(LMM)$coefficients)[numofraw,numofcolumn]
      seofbeta<-data.frame(summary(LMM)$coefficients)[numofraw,numofcolumn2]
      p<-data.frame(summary(LMM)$coefficients)[numofraw,numofcolumn3]
      effect<-exp(beta*10)-1
      effect_l<-exp((beta-1.96*seofbeta)*10)-1
      effect_h<-exp((beta+1.96*seofbeta)*10)-1
      resultofLMM[(a-1)*length(pollutants)*length(timewindow)+(b-1)*length(timewindow)+c,]<-c(biomarkers[a],pollutants[b],timewindow[c],"10",effect,effect_l,effect_h,p)
      c<-c+1
    }
    b<-b+1
  }
  a<-a+1
}

#更改结果表格中biomarkers的名字
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.22.6"]<-"LPA 22:6"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.20.3"]<-"LPA 20:3"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.20.4"]<-"LPA 20:4"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.18.0"]<-"LPA 18:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.18.1"]<-"LPA 18:1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.18.2"]<-"LPA 18:2"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPA.16.0"]<-"LPA 16:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="cPA.18.1"]<-"cPA 18:1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="cPA.16.0"]<-"cPA 16:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LysoPAF.C18"]<-"LPC(O) 18:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LysoPAF.C18.1"]<-"LPC(O) 18:1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LysoPAF.C16"]<-"LPC(O) 16:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPG.18.0"]<-"LPG 18:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPG.18.1"]<-"LPG 18:1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPG.16.0"]<-"LPG 16:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LPG.16.1"]<-"LPG 16:1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LysoPS.18.0"]<-"LysoPS 18:0"
resultofLMM$Biomarker[resultofLMM$Biomarker=="LysoPS.18.1"]<-"LysoPS 18:1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="IL_8_5p"]<-"IL-8"
resultofLMM$Biomarker[resultofLMM$Biomarker=="MCP_1_5p"]<-"MCP-1"
resultofLMM$Biomarker[resultofLMM$Biomarker=="ifn_s"]<-"IFN-γ"
resultofLMM$Biomarker[resultofLMM$Biomarker=="sCD40L_5p"]<-"sCD40L"
resultofLMM$Biomarker[resultofLMM$Biomarker=="nonhdl"]<-"Non-HDL-C"

#为结果增加是否显著列
resultofLMM<-data.frame(resultofLMM)
resultofLMM$effect.mean<-as.numeric(resultofLMM$effect.mean)
resultofLMM$effect.lower.confidence.limit<-as.numeric(resultofLMM$effect.lower.confidence.limit)
resultofLMM$effect.upper.confidence.limit<-as.numeric(resultofLMM$effect.upper.confidence.limit)
resultofLMM$significant<-ifelse((resultofLMM$effect.lower.confidence.limit)*(resultofLMM$effect.upper.confidence.limit)<0,"Nonsignificant",ifelse(resultofLMM$effect.mean>0,"Significant Positive","Significant Negative"))
resultofLMM$group<-rep("Only participants residing within 3 km of the station",length(biomarkers)*2)

#合并
yy<-read.csv("result_PM2.5_F(MA) and all Biomarkers(one-pollutant)_1124.csv",header=TRUE)
yy<-yy[c(which(yy$Time.window..MA.day.==14),which(yy$Time.window..MA.day.==30)),]
yy$group<-rep("All participants",length(biomarkers)*2)
yy<-yy[,-1]
resultofLMM<-rbind(resultofLMM,yy)
write.csv(resultofLMM,"result_PM2.5_F(MA) and all Biomarkers(one-pollutant)_Only participants residing within 3 km.csv")

#ggplot画图:x=MA day，y=effect，facet=biomarker
a_rawdata<-read.csv("result_PM2.5_F(MA) and all Biomarkers(one-pollutant)_Only participants residing within 3 km.csv",header=TRUE,as.is=FALSE)
biomarkers<-rev(c("IL-8","MCP-1","sCD40L",
                  "IFN-γ","Non-HDL-C",
                  "cPA 16:0","cPA 18:1",
                  "LPA 16:0","LPA 18:0","LPA 18:1","LPA 18:2","LPA 20:3","LPA 20:4","LPA 22:6",
                  "LPC(O) 16:0","LPC(O) 18:0","LPC(O) 18:1",
                  "LPG 16:0","LPG 16:1","LPG 18:0","LPG 18:1",
                  "LysoPS 18:0","LysoPS 18:1"))
a_rawdata$Biomarker<-factor(a_rawdata$Biomarker,levels=biomarkers)
a_rawdata$Time.window..MA.day.[a_rawdata$Time.window..MA.day.==14]<-"14-d time window"
a_rawdata$Time.window..MA.day.[a_rawdata$Time.window..MA.day.==30]<-"30-d time window"
a_rawdata$group<-factor(a_rawdata$group,levels=rev(c("All participants","Only participants residing within 3 km of the station")))

library(ggplot2)
library(ggthemes)

ggplot(a_rawdata,aes(effect.mean*100,Biomarker,shape=group,colour=significant))+
  geom_vline(aes(xintercept=0),linetype="dashed")+
  geom_point(size=2,position=position_dodge(width=0.6))+
  geom_errorbar(aes(xmin=effect.lower.confidence.limit*100,xmax=effect.upper.confidence.limit*100),
                width=0.27,position=position_dodge(width=0.6))+
  facet_wrap(~Time.window..MA.day.,scales="free_x",ncol=2)+
  scale_colour_manual(values=c("#B3B8BC","black","black"),guide="none")+
  scale_shape_manual(values=c("All participants"=16,
                              "Only participants residing within 3 km of the station"=1))+
  labs(x=expression(paste("Percent change (%) in biomarkers per 10 μg/m"^"3"," increment of PM"["2.5"])),
       y=NULL)+
  theme_few()+
  theme(legend.position="top",
        legend.title=element_blank(),
        strip.background=element_rect(fill="grey95",colour="black"),
        axis.title.x=element_text(vjust=-2),
        panel.spacing.y=unit(0.16,"in"),
        plot.margin=margin(20,20,20,20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text=element_text(size=12),
        legend.text=element_text(size=12))#字号设置为12 points

ggsave("FigureS3.tiff", width=10.5,height=11)



##Figure S4(mediation_residing within 3 km)
library(lme4)
library(splines)
library(mediation)
detach("package:lmerTest")#卸载lmerTest包，否则mediate函数无法实现

#部分时间窗口(outcome显著效应最高的时间窗口/效应最高的时间窗口），仅LPA
pollutants<-c("PM2.5_F")
biomarkers_mediators<-c("LPA.16.0","LPA.18.0","LPA.18.1","LPA.18.2","LPA.20.4","LPA.22.6",
                        "LysoPAF.C16","LysoPAF.C18","LysoPAF.C18.1",
                        "LPG.16.0","LPG.16.1","LPG.18.0",
                        "LPA.20.3","LysoPS.18.0")
biomarkers_outcomes<-c("IL_8_5p","MCP_1_5p","sCD40L_5p",
                       "ifn_s",
                       "nonhdl")
timewindow<-c(14,14,14,
              14,
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

resultofMediation$group<-rep("Only participants residing within 3 km of the station",nrow(resultofMediation))
#合并
yy<-read.csv("result_mediation(PM2.5_F to LPLs to inflammatory biomarkers_30d)_20231124.csv",header=TRUE)
yy$group<-rep("All participants",nrow(yy))
yy<-yy[,-1]
resultofMediation<-data.frame(resultofMediation)
resultofMediation<-rbind(resultofMediation,yy)
write.csv(resultofMediation,"result_mediation(PM2.5_F to LPLs to inflammatory biomarkers_30d)_Only participants residing within 3 km.csv")

bb<-read.csv("result_mediation(PM2.5_F to LPLs to inflammatory biomarkers_14d)_Only participants residing within 3 km.csv",header=TRUE)
bb<-bb[,-1]
resultofMediation<-rbind(resultofMediation,bb)

write.csv(resultofMediation,"result_mediation(PM2.5_F to LPLs to inflammatory biomarkers_14&30d)_Only participants residing within 3 km.csv")


#ggplot画图:x=MA day，y=effect，facet=biomarker
a_rawdata<-read.csv("result_mediation(PM2.5_F to LPLs to inflammatory biomarkers_14&30d)_Only participants residing within 3 km.csv",header=TRUE,as.is=FALSE)
Outcomes<-c("IL-8","MCP-1","sCD40L","IFN-γ","Non-HDL-C")
Mediators<-rev(c("LPA 16:0","LPA 18:0","LPA 18:1","LPA 18:2","LPA 20:4","LPA 22:6",
                 "LPC(O) 16:0","LPC(O) 18:0","LPC(O) 18:1",
                 "LPG 16:0","LPG 16:1","LPG 18:0",
                 "LPA 20:3","LysoPS 18:0"))
a_rawdata$Outcome<-factor(a_rawdata$Outcome,levels=Outcomes)
a_rawdata$Mediator<-factor(a_rawdata$Mediator,levels=Mediators)
a_rawdata$group<-factor(a_rawdata$group,levels=rev(c("All participants","Only participants residing within 3 km of the station")))
a_rawdata$significant<-ifelse(a_rawdata$Mediation.Proportion.p.Value>0.05,"Nonsignificant",
                              ifelse(a_rawdata$Mediation.Proportion>0,"Significant Positive","Significant Negative"))


library(ggplot2)
library(ggthemes)

ggplot(a_rawdata,aes(Mediation.Proportion*100,Mediator,shape=group,colour=significant))+
  geom_vline(aes(xintercept=0),linetype="dashed")+
  geom_point(size=2,position=position_dodge(width=0.6))+
  geom_errorbar(aes(xmin=Mediation.Proportion.25.*100,xmax=Mediation.Proportion.75.*100),
                width=0.27,position=position_dodge(width=0.6))+
  facet_grid(Time~Outcome,scales="free_x")+
  scale_colour_manual(values=c("#B3B8BC","black","black"),guide="none")+
  scale_shape_manual(values=c("All participants"=16,
                              "Only participants residing within 3 km of the station"=1))+
  labs(x="Mediation Proportion (%)",
       y=NULL)+
  theme_few()+
  theme(legend.position="top",
        legend.title=element_blank(),
        strip.background=element_rect(fill="grey95",colour="black"),
        axis.title.x=element_text(vjust=-2),
        panel.spacing.y=unit(0.16,"in"),
        plot.margin=margin(20,20,20,20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text=element_text(size=12),
        legend.text=element_text(size=12))#字号设置为12 points
ggsave("FigureS4.tiff", width=10.5,height=11)
