#设置工作空间，正斜杠please
setwd("D:\\Desktop\\李不怕呀\\Data & Writing\\LPA et al\\SCOPE\\20231124 update_文章1修改稿")
#ggplot画图:x=MA day，y=effect，facet=biomarker
a_rawdata<-read.csv("result_PM2.5_F(MA) and all Biomarkers(one-pollutant)_1124.csv",header=TRUE,as.is=FALSE)

biomarkers<-c("IL-8","MCP-1","sCD40L",
              "IFN-γ","Non-HDL-C","LysoPS 18:1",
              "cPA 16:0","cPA 18:1",
              "LPA 16:0","LPA 18:0","LPA 18:1","LPA 18:2","LPA 20:3","LPA 20:4","LPA 22:6",
              "LPC(O) 16:0","LPC(O) 18:0","LPC(O) 18:1",
              "LPG 16:0","LPG 16:1","LPG 18:0","LPG 18:1",
              "LysoPS 18:0")
a_rawdata$Biomarker<-factor(a_rawdata$Biomarker,levels=biomarkers)
a_rawdata<-a_rawdata[order(a_rawdata$p.value),]
a_rawdata$BH<-p.adjust(a_rawdata$p.value,method="BH")
a_rawdata$text<-ifelse(a_rawdata$BH<0.05,"*","")

library(ggplot2)
library(ggthemes)

ggplot(a_rawdata,aes(Time.window..MA.day.,effect.mean*100,colour=significant))+
  geom_hline(aes(yintercept=0),linetype="dashed")+
  geom_rect(aes(xmin=13.5,xmax=30.5,ymin=-Inf,ymax=Inf),colour="white",fill="#CD534C",alpha=0.002)+
  geom_point(shape=15,size=2)+
  geom_errorbar(aes(ymin=effect.lower.confidence.limit*100,ymax=effect.upper.confidence.limit*100),
                width=0.27)+
  geom_text(aes(y=ifelse(effect.mean<0,effect.upper.confidence.limit*100-effect.lower.confidence.limit*15,effect.lower.confidence.limit*100-effect.upper.confidence.limit*15),label=text),col="black",size=5)+
  facet_wrap(~Biomarker,scales="free_y",ncol=6)+
  scale_x_continuous(breaks=c(1,7,14,30))+
  scale_colour_manual(values=c("#B3B8BC","black","black"))+
  labs(x="Avg days",
       y=expression(paste("Percent change (%) in biomarkers per 10 μg/m"^"3"," increment of PM"["2.5"])))+
  theme_few()+
  theme(legend.position="none",
        axis.title.x=element_text(vjust=-6),
        axis.title.y=element_text(vjust=8),
        axis.text.x=element_text(angle=90,vjust=0.5),
        panel.spacing.y=unit(0.16,"in"),
        plot.margin=margin(20,20,20,20),
        axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        strip.text=element_text(size=12),
        legend.text=element_text(size=12))#字号设置为12 points
ggsave("Figure1_draft_20231130.pdf", width=10.5,height=10.5)