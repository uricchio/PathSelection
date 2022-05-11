library(cowplot)
library(ggplot2)
library(wesanderson)
library(dplyr)

read.table("~/projects/elisa/scriptsWSensTest/comp_across_reps.txt")->compRep
read.table("~/projects/elisa/scriptsWSensTest/comp_num_hits.txt")->compHits
read.table("~/projects/elisa/scriptsWSensTest/total_hits.txt")->totalHits

compRep %>% rowwise() %>% mutate(sameSignSig = ifelse(sign(V2)==sign(V3) & V6 > 0,1,0)) -> newCompRep


plA<-ggplot(data=newCompRep,aes(V2,V3,col=as.factor(sameSignSig)))+geom_point()+xlim(-55,55)+ylim(-55,55)+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                    xlab(expression(italic("2Ns (N = 92)")))+ylab(expression(italic("2Ns (N > 92)")))+
                    scale_color_manual(values=wes_palette("Royal1"))+theme(legend.position = "None")
                
plB<-ggplot(data=totalHits,aes(V1,V2))+geom_point()+scale_x_log10()+scale_y_log10()+xlab(expression(italic(N)))+ylab(expression(italic("Number of selected loci")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_smooth(method="lm")

plot_grid(plA,plB,labels=c("A","B"))

ggsave("~/projects/elisa/scriptsWSensTest/Figures/FigS1.png",width=9,height=3)


# plot traj

plot_traj<-function(pos,samp,pass) {
	
  read.table(paste("~/projects/elisa/MCMCout/MCMCinit35rate5.0mutTypeINDEL/pos",pos,".samp", samp,".",pass,".traj.traj.gz",sep=""))->traj
  read.table(paste("~/projects/elisa/MCMCout/MCMCinit35rate5.0mutTypeINDEL/pos",pos,".samp", samp,".",pass,".traj.time.gz",sep=""))->time

  df<-data.frame(time=c(t(time[2,]))[2:length(c(t(time[2,])))],traj=c(t(traj[2,]))[2:length(c(t(traj[2,])))],r=1)
  i <-3
  while ( i< length(traj$V1)){
      df<-rbind(df,data.frame(time=c(t(time[i,]))[2:length(c(t(time[i,])))],traj=c(t(traj[i,]))[2:length(c(t(traj[i,])))],r=i-1))
      i <- i + 1
  }
  
  read.table(paste("/Users/uricchio/projects/elisa/data/trajexp5.0.init35.mutTypeINDEL/pos",pos,".samp", samp, ".",pass,".traj",sep=""))->trajObs
  #print(trajObs)

  trajObs<-cbind(trajObs,r=rep(1000,length(trajObs$V1)))
 
  pl<-ggplot(df,aes(time,traj,group=r))+geom_line(alpha=0.3,col='gray')+xlab("time")+ylab("allele frequency")
  pl<-pl+geom_point(data=trajObs,aes(V3,V1/V2,group=r),col='red')

  return(pl)

}

plA<-plot_traj(89033,9,3)
plB<-plot_traj(89033,9,2)
plC<-plot_traj(89033,9,1)


plot_grid(plA,plB,plC,labels=c("A","B","C"),ncol=3)
ggsave("~/projects/elisa/scriptsWSensTest/Figures/FigS5.png",width=9,height=2.5)
