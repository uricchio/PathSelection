library(ggplot2)
library(cowplot)

ret_sig<-function(l0,l1,l2) {
  return(l0[sign(l0$V3) == sign(l0$V4) & sign(l1$V3) == sign(l1$V4) & sign(l2$V3) == sign(l2$V4)])
}

myPlot<-function(tval) {
  pl<-ggplot(tval,aes(V1,V2,color=as.factor(V5)))+geom_point()+geom_errorbar(data=tval,aes(ymin=V3,ymax=V4,color=as.factor(V5)))+xlab("position")+ylab(expression(italic("2Ns")))+ylim(c(-200,200))+scale_colour_manual(values = c("black", "red", "blue","gray"),labels=c("NS","+","-"),name="")
  return(pl)
}


read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp2.1.traj.param.data")->t2.1
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp2.1.traj.param.data")->t2.1I
t2.1<-rbind(t2.1,t2.1I)
plA<-myPlot(t2.1)

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp2.2.traj.param.data")->t2.2
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp2.2.traj.param.data")->t2.2I
t2.2<-rbind(t2.2,t2.2I)
plB<-myPlot(t2.2)

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp2.3.traj.param.data")->t2.3
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp2.3.traj.param.data")->t2.3I
t2.3<-rbind(t2.3,t2.3I)
plC<-myPlot(t2.3)

pl2<-plot_grid(plA,plB,plC,ncol=1)
pl2
ggsave("~/projects/elisa/plots/condBOTH2.pdf")

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp9.1.traj.param.data")->t9.1
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp9.1.traj.param.data")->t9.1I
t9.1<-rbind(t9.1,t9.1I)
plA<-myPlot(t9.1)

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp9.2.traj.param.data")->t9.2
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp9.2.traj.param.data")->t9.2I
t9.2<-rbind(t9.2,t9.2I)
plB<-myPlot(t9.2)

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp9.3.traj.param.data")->t9.3
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp9.3.traj.param.data")->t9.3I
t9.3<-rbind(t9.3,t9.3I)
plC<-myPlot(t9.3)

pl9<-plot_grid(plA,plB,plC,ncol=1)
pl9
ggsave("~/projects/elisa/plots/condBOTH9.pdf")



read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp17.1.traj.param.data")->t17.1
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp17.1.traj.param.data")->t17.1I
t17.1<-rbind(t17.1,t17.1I)
plA<-myPlot(t17.1)

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp17.2.traj.param.data")->t17.2
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp17.2.traj.param.data")->t17.2I
t17.2<-rbind(t17.2,t17.2I)
plB<-myPlot(t17.2)

read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeSNP.samp17.3.traj.param.data")->t17.3
read.table("~/projects/elisa/data/analysis/MCMCinit35.rate5.0.mutTypeINDEL.samp17.3.traj.param.data")->t17.3I
t17.3<-rbind(t17.3,t17.3I)
plC<-myPlot(t17.3)

plC<-myPlot(t17.3)

pl17<-plot_grid(plA,plB,plC,ncol=1)
pl17
ggsave("~/projects/elisa/plots/condBOTH17.pdf")


