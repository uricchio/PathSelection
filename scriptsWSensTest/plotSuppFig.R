library(cowplot)
library(ggplot2)
library(wesanderson)

read.table("~/projects/elisa/scriptsWSensTest/comp_across_reps.txt")->compRep
read.table("~/projects/elisa/scriptsWSensTest/comp_num_hits.txt")->compHits
read.table("~/projects/elisa/scriptsWSensTest/total_hits.txt")->totalHits

compRep %>% rowwise() %>% mutate(sameSignSig = ifelse(sign(V2)==sign(V3) & V6 > 0,1,0)) -> newCompRep


plA<-ggplot(data=newCompRep,aes(V2,V3,col=as.factor(sameSignSig)))+geom_point()+xlim(-100,100)+ylim(-100,100)+
              theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
                    xlab(expression(italic("2Ns (N = 425)")))+ylab(expression(italic("2Ns (N > 425)")))+
                    scale_color_manual(values=wes_palette("Royal1"))+theme(legend.position = "None")
                


plB<-ggplot(data=totalHits,aes(V1,V2))+geom_point()+scale_x_log10()+scale_y_log10()+xlab(expression(italic(N)))+ylab(expression(italic("Number of selected loci")))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+geom_smooth(method="lm")

plot_grid(plA,plB,labels=c("A","B"))

ggsave("~/projects/elisa/scriptsWSensTest/Figures/FigS1.png",width=9,height=3)


read.table("~/projects")
