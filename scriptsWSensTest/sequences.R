###Load Libraries
library(ggplot2)

###set working directory and load data
setwd("/Users/uricchio/projects/elisa/data")
data<-read.csv(file='Specialization SNPs.csv', header=T, stringsAsFactors = FALSE)
View(data)

###adding passage numbers, replicate IDs, and evolution line IDs
data$passage<-NA
data$passage[which(data$SAMPLE == 'Sample37')] <- 0
data$passage[which((data$SAMPLE == 'Sample1') |(data$SAMPLE == 'Sample2') | (data$SAMPLE =='Sample3') | (data$SAMPLE =='Sample4') | (data$SAMPLE =='Sample5') | (data$SAMPLE =='Sample6') | (data$SAMPLE =='Sample7') | (data$SAMPLE =='Sample8') | (data$SAMPLE =='Sample9'))] <- 1
data$passage[which((data$SAMPLE == 'Sample10') |(data$SAMPLE == 'Sample11') | (data$SAMPLE =='Sample12') | (data$SAMPLE =='Sample13') | (data$SAMPLE =='Sample14') | (data$SAMPLE =='Sample15') | (data$SAMPLE =='Sample16') | (data$SAMPLE =='Sample17') | (data$SAMPLE =='Sample18'))] <- 2
data$passage[which((data$SAMPLE == 'Sample19') |(data$SAMPLE == 'Sample20') | (data$SAMPLE =='Sample21') | (data$SAMPLE =='Sample22') | (data$SAMPLE =='Sample23') | (data$SAMPLE =='Sample24') | (data$SAMPLE =='Sample25') | (data$SAMPLE =='Sample26') | (data$SAMPLE =='Sample27'))] <- 3
data$passage[which((data$SAMPLE == 'Sample28') |(data$SAMPLE == 'Sample29') | (data$SAMPLE =='Sample30') | (data$SAMPLE =='Sample31') | (data$SAMPLE =='Sample32') | (data$SAMPLE =='Sample33') | (data$SAMPLE =='Sample34') | (data$SAMPLE =='Sample35') | (data$SAMPLE =='Sample36'))] <- 4
data$rep<-NA
data$rep[which((data$SAMPLE == 'Sample1') |(data$SAMPLE == 'Sample10') | (data$SAMPLE =='Sample19') | (data$SAMPLE =='Sample28'))] <- 2.1
data$rep[which((data$SAMPLE == 'Sample2') |(data$SAMPLE == 'Sample11') | (data$SAMPLE =='Sample20') | (data$SAMPLE =='Sample29'))] <- 2.2
data$rep[which((data$SAMPLE == 'Sample3') |(data$SAMPLE == 'Sample12') | (data$SAMPLE =='Sample21') | (data$SAMPLE =='Sample30'))] <- 2.3
data$rep[which((data$SAMPLE == 'Sample4') |(data$SAMPLE == 'Sample13') | (data$SAMPLE =='Sample22') | (data$SAMPLE =='Sample31'))] <- 9.1
data$rep[which((data$SAMPLE == 'Sample5') |(data$SAMPLE == 'Sample14') | (data$SAMPLE =='Sample23') | (data$SAMPLE =='Sample32'))] <- 9.2
data$rep[which((data$SAMPLE == 'Sample6') |(data$SAMPLE == 'Sample15') | (data$SAMPLE =='Sample24') | (data$SAMPLE =='Sample33'))] <- 9.3
data$rep[which((data$SAMPLE == 'Sample7') |(data$SAMPLE == 'Sample16') | (data$SAMPLE =='Sample25') | (data$SAMPLE =='Sample34'))] <- 17.1
data$rep[which((data$SAMPLE == 'Sample8') |(data$SAMPLE == 'Sample17') | (data$SAMPLE =='Sample26') | (data$SAMPLE =='Sample35'))] <- 17.2
data$rep[which((data$SAMPLE == 'Sample9') |(data$SAMPLE == 'Sample18') | (data$SAMPLE =='Sample27') | (data$SAMPLE =='Sample36'))] <- 17.3

data$Line<-NA
data$Line[which((data$rep == '2.1') |(data$rep == '2.2') | (data$rep == '2.3'))] <- 2
data$Line[which((data$rep == '9.1') |(data$rep == '9.2') | (data$rep == '9.3'))] <- 9
data$Line[which((data$rep == '17.1') |(data$rep == '17.2') | (data$rep == '17.3'))] <- 17

#making frequency, and passage numeric
data$FREQ<-gsub("%", "", paste(data$FREQ))
data$FREQ<-as.numeric(data$FREQ)
data$passage<-as.numeric(data$passage)

#make unique SNP by replicate IDs so that they can be looked at over time
data$ID<-paste(data$POS,data$rep)

#make a graph that shows where the variable regions are over time
figure3<-ggplot(data=data, aes(x=POS, y=FREQ))+facet_wrap(~passage)+geom_point(aes(color=rep))+theme_light()

#just look at that region of high variability (and get rid of fixed alleles)
datanofixed<-data[which(data$FREQ <98),]
datazoom1<-datanofixed[which(datanofixed$POS>32677),]
datazoom<-datazoom1[which(datazoom1$POS<35697),]
figurezoom<-ggplot(data=datazoom, aes(x=passage, y=FREQ, group=ID))+geom_line(aes(color=factor(rep)))+theme_light()+scale_x_discrete(breaks=c(0,1,4,6,9))
datazoomred<-datazoom[which(datazoom$FREQ>10),]
figurezoomr<-ggplot(data=datazoomred, aes(x=passage, y=FREQ, group=ID))+geom_line(aes(color=factor(POS)))+theme_light()+scale_x_discrete(breaks=c(0,1,4,6,9))
 
#i should also check the indels since all the papers mostly talk about function through indels
#add all the same labels as above
datain<-read.csv(file='Specialization indel.csv', header=T, stringsAsFactors = FALSE)
datain$passage<-NA
datain$passage[which(datain$SAMPLE == 'Sample37')] <- 0
datain$passage[which((datain$SAMPLE == 'Sample1') |(datain$SAMPLE == 'Sample2') | (datain$SAMPLE =='Sample3') | (datain$SAMPLE =='Sample4') | (datain$SAMPLE =='Sample5') | (datain$SAMPLE =='Sample6') | (datain$SAMPLE =='Sample7') | (datain$SAMPLE =='Sample8') | (datain$SAMPLE =='Sample9'))] <- 1
datain$passage[which((datain$SAMPLE == 'Sample10') |(datain$SAMPLE == 'Sample11') | (datain$SAMPLE =='Sample12') | (datain$SAMPLE =='Sample13') | (datain$SAMPLE =='Sample14') | (datain$SAMPLE =='Sample15') | (datain$SAMPLE =='Sample16') | (datain$SAMPLE =='Sample17') | (datain$SAMPLE =='Sample18'))] <- 4
datain$passage[which((datain$SAMPLE == 'Sample19') |(datain$SAMPLE == 'Sample20') | (datain$SAMPLE =='Sample21') | (datain$SAMPLE =='Sample22') | (datain$SAMPLE =='Sample23') | (datain$SAMPLE =='Sample24') | (datain$SAMPLE =='Sample25') | (datain$SAMPLE =='Sample26') | (datain$SAMPLE =='Sample27'))] <- 6
datain$passage[which((datain$SAMPLE == 'Sample28') |(datain$SAMPLE == 'Sample29') | (datain$SAMPLE =='Sample30') | (datain$SAMPLE =='Sample31') | (datain$SAMPLE =='Sample32') | (datain$SAMPLE =='Sample33') | (datain$SAMPLE =='Sample34') | (datain$SAMPLE =='Sample35') | (datain$SAMPLE =='Sample36'))] <- 9
datain$rep<-NA
datain$rep[which((datain$SAMPLE == 'Sample1') |(datain$SAMPLE == 'Sample10') | (datain$SAMPLE =='Sample19') | (datain$SAMPLE =='Sample28'))] <- 2.1
datain$rep[which((datain$SAMPLE == 'Sample2') |(datain$SAMPLE == 'Sample11') | (datain$SAMPLE =='Sample20') | (datain$SAMPLE =='Sample29'))] <- 2.2
datain$rep[which((datain$SAMPLE == 'Sample3') |(datain$SAMPLE == 'Sample12') | (datain$SAMPLE =='Sample21') | (datain$SAMPLE =='Sample30'))] <- 2.3
datain$rep[which((datain$SAMPLE == 'Sample4') |(datain$SAMPLE == 'Sample13') | (datain$SAMPLE =='Sample22') | (datain$SAMPLE =='Sample31'))] <- 9.1
datain$rep[which((datain$SAMPLE == 'Sample5') |(datain$SAMPLE == 'Sample14') | (datain$SAMPLE =='Sample23') | (datain$SAMPLE =='Sample32'))] <- 9.2
datain$rep[which((datain$SAMPLE == 'Sample6') |(datain$SAMPLE == 'Sample15') | (datain$SAMPLE =='Sample24') | (datain$SAMPLE =='Sample33'))] <- 9.3
datain$rep[which((datain$SAMPLE == 'Sample7') |(datain$SAMPLE == 'Sample16') | (datain$SAMPLE =='Sample25') | (datain$SAMPLE =='Sample34'))] <- 17.1
datain$rep[which((datain$SAMPLE == 'Sample8') |(datain$SAMPLE == 'Sample17') | (datain$SAMPLE =='Sample26') | (datain$SAMPLE =='Sample35'))] <- 17.2
datain$rep[which((datain$SAMPLE == 'Sample9') |(datain$SAMPLE == 'Sample18') | (datain$SAMPLE =='Sample27') | (datain$SAMPLE =='Sample36'))] <- 17.3

datain$Line<-NA
datain$Line[which((datain$rep == '2.1') |(datain$rep == '2.2') | (datain$rep == '2.3'))] <- 2
datain$Line[which((datain$rep == '9.1') |(datain$rep == '9.2') | (datain$rep == '9.3'))] <- 9
datain$Line[which((datain$rep == '17.1') |(datain$rep == '17.2') | (datain$rep == '17.3'))] <- 17
datain$FREQ<-gsub("%", "", paste(datain$FREQ))
datain$FREQ<-as.numeric(datain$FREQ)
datain$ID<-paste(datain$POS,datain$rep)
figure2in<-ggplot(data=datain, aes(x=POS, y=FREQ))+facet_wrap(~passage)+geom_point()+theme_light()
figure2in

####to see the frequencies better, I wanted to make a column that calculates the change in frequency
####since the previous time point. This loop takes each position in each replicate and takes the difference
####in that position since the previous replicate. (If previous was NA bc poor sequence quality, it takes the 
####difference since the last good value). This also adds in the 0 passage data as the starting point.It also
####adds a sumdelta column which is the sum of the absolute values of deltafreq (change in frequency). This is 
####so that we can filter

####WARNING: for loop to add deltas starts here. It takes hours. 
datapass<-subset(data, passage>0)
u<-unique(datapass$ID)
data$deltafreq<- NA
data$sumdelta<-NA
datawdelta2<-data.frame()

for(i in 1:length(u)){
  z1<-subset(data, ID==u[i])
  z2<-subset(data, passage==0 & POS==z1$POS[1])
  z3<-rbind(z1, z2)
  z4 <- z3[order(z3$passage),]
  z <-z4[!is.na(z4$FREQ),]
  if (sum(z$FREQ, na.rm=TRUE)>0){
    z$deltafreq[1]<-0
    if (length(z$ID)>1){
  for(y in 2:length(z$ID)){
    z$deltafreq[y] <- z$FREQ[y] - z$FREQ[y-1]}}
    z$sumdelta<-sum(abs(z$deltafreq))
    datawdelta2 <-rbind(datawdelta2,z)}
    
    print(i)
}

#making a figure that plots deltafreqs for the ones with sumdeltas >10 (i.e. changes at least a little over time)
datachange1<-subset(datawdelta2, sumdelta>10)
datachange<-subset(datachange1, !is.na(rep))
moonrise9<-c("#e0a295", "#b08e23", "#bcb46d", "#e9ba4e", "#d0e4e2", "#15B2D3", "#9C2706", "#deb95e", "#DE5086")
figurechange<-ggplot(data=datachange, aes(x=passage, y=deltafreq, group=ID))+geom_line(aes(color=factor(rep), alpha=.8), show.legend = FALSE)+theme_minimal()+facet_wrap(~rep)+ylab("Change in Frequency")+xlab("Passage Number")+scale_color_manual(values=moonrise9)
figurechange
moonrise10<-c("#e0a295", "#b08e23", "#bcb46d", "#e9ba4e", "#d0e4e2", "#15B2D3", "#9C2706", "#deb95e", "#DE5086", "#e0a295")
data2<-data
data2$rep[is.na(data2$rep)] <- 0
figure2<-ggplot(data=data2, aes(x=POS, y=FREQ))+facet_wrap(~passage)+geom_point(aes(color=factor(rep)))+theme_minimal()+scale_color_manual(values=moonrise10)+xlab("Genome Position")+ylab("Frequency")

