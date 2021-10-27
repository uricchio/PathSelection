library(ggplot2)
library(cowplot)

ret_sig<-function(l0,l1,l2) {
  return(l0[sign(l0$V3) == sign(l0$V4) & sign(l1$V3) == sign(l1$V4) & sign(l2$V3) == sign(l2$V4)])
}

myPlot<-function(tval) {
  pl<-ggplot(tval,aes(V1,V2,color=as.factor(V5)))+geom_point()+geom_errorbar(data=tval,aes(ymin=V3,ymax=V4,color=as.factor(V5)))+xlab("position")+ylab(expression(italic("2Ns")))+ylim(c(-200,200))+scale_colour_manual(values = c("black", "red", "blue","gray"),labels=c("NS","+","-"),name="")
  return(pl)
}



read.table("~/projects/elisa/data/analysis/MCMC.samp2.1.traj.param.data")->t2.1
plA<-myPlot(t2.1)

read.table("~/projects/elisa/data/analysis/MCMC.samp2.2.traj.param.data")->t2.2
plB<-myPlot(t2.2)

read.table("~/projects/elisa/data/analysis/MCMC.samp2.3.traj.param.data")->t2.3
plC<-myPlot(t2.3)

pl2<-plot_grid(plA,plB,plC,ncol=1)
ggsave("~/projects/elisa/plots/cond2.pdf")

read.table("~/projects/elisa/data/analysis/MCMC.samp9.1.traj.param.data")->t9.1
plA<-myPlot(t9.1)

read.table("~/projects/elisa/data/analysis/MCMC.samp9.2.traj.param.data")->t9.2
plB<-myPlot(t9.2)

read.table("~/projects/elisa/data/analysis/MCMC.samp9.3.traj.param.data")->t9.3
plC<-myPlot(t9.3)

pl9<-plot_grid(plA,plB,plC,ncol=1)
ggsave("~/projects/elisa/plots/cond9.pdf")

read.table("~/projects/elisa/data/analysis/MCMC.samp17.1.traj.param.data")->t17.1
plA<-myPlot(t17.1)

read.table("~/projects/elisa/data/analysis/MCMC.samp17.2.traj.param.data")->t17.2
plB<-myPlot(t17.2)

read.table("~/projects/elisa/data/analysis/MCMC.samp17.3.traj.param.data")->t17.3
plC<-myPlot(t17.3)

pl17<-plot_grid(plA,plB,plC,ncol=1)
ggsave("~/projects/elisa/plots/cond17.pdf")

# plotting the one SNP that is same sign and sig acrross all three experiments
# actually just looking at the trajectories (in ../data/traj/pos27870.samp9.1.traj etc) can see that this is nonsense and not useful

