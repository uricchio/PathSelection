
# Load packages

library(ggplot2)
library(lme4)
library(afex)
library(emmeans)

# Load full data


IFD <- read.csv(file = 'InfectivityDataFull2.csv',
                header=T, stringsAsFactors = FALSE)


# Check read-in

#head(IFD)

# Create a binary preditor of whether evolved line = assay line 
# (infection is familiar or foriegn genotype for pathogen population)

IFD$Self <- IFD$PlodiaAssayLine == IFD$EvolvedLine

#Create an average count metric that calculates the number of virions per infected cadaver on average
IFD$AveCount <-IFD$PoolCount/IFD$NumberInfected
# note this creates 2 typed of NAs:
# NA where data is truly missing (was never measured or lost)
# NA where none were infected, so poolcount couldn't be counted

#make a proportion infected metric
IFD$PropInf <-IFD$NumberInfected/(IFD$NumberUninfected+IFD$NumberInfected)

#make an general fitness metric by multiplying average count and proportion infected
IFD$VirFit <-IFD$AveCount*IFD$PropInf
# for the 'second type NAs (see above) change fitness = 0
# this is still broadly biologically correct and necessary for downstream modelling 
IFD$VirFit[which(IFD$NumberInfected == 0)] <- 0


# Now convert appropriately to factors where needed
  
IFD$PlodiaAssayLine <- factor(IFD$PlodiaAssayLine, ordered = FALSE)

IFD$EvolvedLine <- factor(IFD$EvolvedLine, ordered = FALSE)

IFD$VirusLine <- factor(IFD$VirusLine, ordered = FALSE)

# Convert the count data to be discrete counts as needed

IFD$PoolCount <- round(IFD$PoolCount, digits = 0)
IFD$AveCount <- round(IFD$AveCount, digits = 0)
IFD$VirFit <- round(IFD$VirFit, digits = 0)

################
## Begin with 'end of evolution' analysis (Figs. 1-3 main manuscript)
################

# Create plotting colours
sunset5<-c("goldenrod1", "darkorange", "orangered3", "deeppink4", "gray4" )
sunset3<-c("goldenrod1",  "orangered3",  "violetred4" )

# Pull out last-passage data for first batch of analysis ('end of evolution')

IFDE <- IFD[which(IFD$PassageNumber == max(IFD$PassageNumber)),]

# First hypothesis: viruses evolved on a line are more likely to be able to infect that line than those evolved on a different line

# 3 replicate lines (VLine) nested under each 'evolved on' treatment.
# every replicate tested on every assay line

# Assay line as a fixed effect (host factor - we know these populations differ in resistance)
# Dose as a fixed effect (for obvious reason)
# 'Self' is the fixed effect of main interest
# 'evolved line' is also a fixed effect of secondary interest
# Random effects: Virus line nested under evolved line as they're replicate populations under same treatment

mixed(cbind(NumberInfected, NumberUninfected) ~ Dose + PlodiaAssayLine + EvolvedLine + Self + (1|EvolvedLine/VirusLine) , 
      data = IFDE, family = binomial, method = 'LRT')

# Significant effects of Dose, Assay Line (both anticipated / expected)
# No significant effect of evolved line (no treatment leads to universally more or less infective virus)
# Significant effect of self - check direction

MMod1 <- glmer(cbind(NumberInfected, NumberUninfected) ~ Dose + PlodiaAssayLine + EvolvedLine + Self + (1|EvolvedLine/VirusLine) , 
               data = IFDE, family = binomial)

summary(MMod1)

# 'SelfTRUE' has positive coefficient (see Dose also) - hypothesis confirmed
# viruses evolved on a line are more likely to be able to infect that line than those evolved on a different line

#lets plot this using emmeans, make model with Evolved Line and Assay Line interaction to be able to graph it 
MMod1a <- glmer(cbind(NumberInfected, NumberUninfected) ~ Dose + EvolvedLine*PlodiaAssayLine + (1|EvolvedLine/VirusLine) , 
               data = IFDE, family = binomial)
graphmod1 <- emmeans(MMod1a, ~ EvolvedLine*PlodiaAssayLine)
graphmod1 <-as.data.frame(graphmod1)

#change EvolvedOn Labels for facet wrap
graphmod1$LabEvolvedLine<- NA
graphmod1$LabEvolvedLine[which(graphmod1$EvolvedLine == '2')]<-'Evolved on Host Genotype 2'
graphmod1$LabEvolvedLine[which(graphmod1$EvolvedLine == '9')]<-'Evolved on Host Genotype 9'
graphmod1$LabEvolvedLine[which(graphmod1$EvolvedLine == '17')]<-'Evolved on Host Genotype 17'

#Figure for End of Evolution Proportion Infected
figure1<-ggplot(data=graphmod1, aes(x=PlodiaAssayLine, y=emmean))+geom_point(aes(color=factor(PlodiaAssayLine)))+geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL, color=factor(PlodiaAssayLine)))+
  scale_color_manual(values=sunset3)+facet_wrap(~LabEvolvedLine)+
  theme_light()+labs(x="Assay Genotype", y="Effect on Proportion Infected", color="Assay Genotype") + guides(fill = FALSE, alpha=FALSE)

figure1

# Second hypothesis: viruses evolved on a line produce more virions on that line than those evolved on a different line
# Restrict analysis to max dose data as virions counts only for those doses
# Otherwise, same anlysis as above, different response variable (and family)

IFDEMD <- IFDE[which(IFDE$Dose == max(IFDE$Dose)),]

mixed(AveCount ~ PlodiaAssayLine + EvolvedLine + Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = poisson, method = 'LRT')

# Significant effect of Assay Line (also arguably expected)
# Possible effect of evolved line
# Significant effect of self - check direction

MMod2 <- glmer(AveCount ~ PlodiaAssayLine + EvolvedLine + Self + (1|EvolvedLine/VirusLine) , 
               data = IFDEMD, family = poisson)

summary(MMod2)

# 'SelfTRUE' has positive value - hypothesis confirmed, ELine 17 does have higher productivity
# viruses evolved on a line produce more virions on that line than those evolved on a different line
#lets plot this using emmeans
MMod2a <- glmer(AveCount ~ PlodiaAssayLine*EvolvedLine + (1|EvolvedLine/VirusLine) , 
                data = IFDEMD, family = poisson)
graphmod2 <- emmeans(MMod2a, ~ PlodiaAssayLine*EvolvedLine)
graphmod2 <-as.data.frame(graphmod2)

#change evolved on labels for facet wrap 
graphmod2$LabEvolvedLine<- NA
graphmod2$LabEvolvedLine[which(graphmod2$EvolvedLine == '2')]<-'Evolved on Host Genotype 2'
graphmod2$LabEvolvedLine[which(graphmod2$EvolvedLine == '9')]<-'Evolved on Host Genotype 9'
graphmod2$LabEvolvedLine[which(graphmod2$EvolvedLine == '17')]<-'Evolved on Host Genotype 17'

#figure for end of evolution viral productivity
figure2<-ggplot(data=graphmod2, aes(x=PlodiaAssayLine, y=emmean))+geom_point(aes(color=factor(PlodiaAssayLine)))+geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL, color=factor(PlodiaAssayLine)))+
  scale_color_manual(values=sunset3)+facet_wrap(~LabEvolvedLine)+
  theme_light()+labs(x="Assay Genotype", y="Effect on Viral Productivity", color="Assay Genotype") + guides(fill = FALSE, alpha=FALSE)

figure2

# hypothesis 3 - confirmation that by virtue of hypothesis 1 and 2, evolution on a line leads to higher fitness on that line
# same analysis as for hypothesis 2

mixed(VirFit ~ PlodiaAssayLine + EvolvedLine + Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = poisson, method = 'LRT')

# Significant effect of Assay Line (also arguably expected)
# No significant effect of evolved line
# Significant effect of self - check direction

MMod3 <- glmer(VirFit ~ PlodiaAssayLine + EvolvedLine + Self + (1|EvolvedLine/VirusLine) , 
               data = IFDEMD, family = poisson)

summary(MMod3)

# As expected, extremely strong positive effect of 'Self' on fitness
#lets plot this using emmeans
MMod3a <- glmer(VirFit ~ PlodiaAssayLine*EvolvedLine + (1|EvolvedLine/VirusLine) , 
                data = IFDEMD, family = poisson)
graphmod3 <- emmeans(MMod3a, ~ PlodiaAssayLine*EvolvedLine)
graphmod3 <-as.data.frame(graphmod3)
#changing labels... again
graphmod3$LabEvolvedLine<- NA
graphmod3$LabEvolvedLine[which(graphmod3$EvolvedLine == '2')]<-'Evolved on Host Genotype 2'
graphmod3$LabEvolvedLine[which(graphmod3$EvolvedLine == '9')]<-'Evolved on Host Genotype 9'
graphmod3$LabEvolvedLine[which(graphmod3$EvolvedLine == '17')]<-'Evolved on Host Genotype 17'

#figure for viral fitness end of evolution
figure3<-ggplot(data=graphmod3, aes(x=PlodiaAssayLine, y=emmean))+geom_point(aes(color=factor(PlodiaAssayLine)))+geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL, color=factor(PlodiaAssayLine)))+
  scale_color_manual(values=sunset3)+facet_wrap(~LabEvolvedLine)+
  theme_light()+labs(x="Assay Genotype", y="Effect on Viral Fitness", color="Assay Genotype") + guides(fill = FALSE, alpha=FALSE)

figure3

# Further investigation - 
# do we see variation between the replicate lines beyond just the 'self' effects and
# did each replicate population differ in its ability to evolve higher fitness on its own Plodia line?
# allow fixed effects of assay line, self, and Virus Line, and allow interaction between the two

mixed(VirFit ~ PlodiaAssayLine + VirusLine*Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = poisson, method = 'LRT')

# Yes, extremely strong evidence of differences in the effect of 'self' depending on which replicate line
 
# MMod4 <- glmer(VirFit ~ PlodiaAssayLine + VirusLine*Self + (1|EvolvedLine/VirusLine) , 
#         data = IFDEMD, family = poisson)
# 
# summary(MMod4)

# Is this effect apparent in both constituent parts of the fitness metric?

mixed(AveCount ~ PlodiaAssayLine + VirusLine*Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = poisson, method = 'LRT')


mixed(cbind(NumberInfected, NumberUninfected) ~ PlodiaAssayLine + VirusLine*Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = binomial, method = 'LRT')

mixed(cbind(NumberInfected, NumberUninfected) ~ PlodiaAssayLine + VirusLine + Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = binomial, method = 'LRT')

# Driven entirely by the effect on virion production following successful infection
# No interaction effect between VL and Self and no effect of virus line alone on infection

#see if there is also an interaction with EvolvedLine
mixed(VirFit ~ PlodiaAssayLine + EvolvedLine*Self + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = poisson, method = 'LRT')
#Interaction is also significant here, indicating that the evolutinary background drives some of these trends, not just virus line

# Adjusted pairwise comparisons of the interacting effects to see what might be driving the fitness interaction
# use bonferroni correction (harsh)

#emmeans package abuse 

MMod4 <- glmer(VirFit ~ PlodiaAssayLine + VirusLine*Self + (1|EvolvedLine/VirusLine) , 
data = IFDEMD, family = poisson)

PWC1 <- contrast(emmeans(MMod4, ~ VirusLine * Self), interaction = TRUE, method = "pairwise", by = "Self", adjust = 'bonferroni')

PWC1

PWC2 <- contrast(emmeans(MMod4, ~ VirusLine * Self), interaction = TRUE, method = "pairwise", by = "VirusLine", adjust = 'bonferroni')

PWC2


# interpretation #
# look down p values and at 'estimate' to get significance of comparison and direction
# all of this controls for other factors, with what is being controlled for either above the comps or elsewhere in the model

# Within-lines (comparing effect of 'self' for each line) THIS IS THE WAY MORE INTERESTING ONE
# Line 9.1 did significantly worse on a hypothetical average ('h.a.') familiar line than a foreign one (p < 0.001)
# Lines 9.2 and 2.2 did not significantly differ in fitness between h.a. familiar and foreign lines(p = 0.50 and p = 0.28 resp.)
# All other lines did significantly better on h.a. familiar lines than a foreign ones (p = 0.0016 for 9.3, p < 0.001 for 2.1,2.3, 17.1, 17.2, 17.3)
# TREND OF SELF = BETTER NOT TRUE FOR 3 LINES, REVERSED FOR 1

# Within familiar/foreign THIS BASICALLY JUST LOOKS AT HOW MUCH DID THE LINES DIFFER 
# on a h.a. familiar host, only 4 / 36 pairwise comparisons showed no significant difference
# on a h.a. foriegn host, 8/36 pairwise comparisons showed no significant difference
# only one of these pairwise comparisons was the same (2.2 and 9.2) which makes total sense given the above
# ALMOST ALL LINES DIFFERED IN THEIR FITNESS FROM *EVERY OTHER LINE* BASICALLY
# even when lines didn't differ on a h.a. foriegn host genotype, they typically did on a h.a. familiar


# Use these contrasts to investigate a correlation between fitness on foreign genotypes vs fitness on own genotype

# Note this is *extremely brittle* from hard-coding

FitCorr <- data.frame(Line = levels(IFDEMD$VirusLine))

FitCorr$Foreign <- append(x = 0, values = as.data.frame(PWC1)$estimate[1:8], after = 1)

FitCorr$Familiar <- append(x = 0, values = as.data.frame(PWC1)$estimate[37:44], after = 1)

cor.test(FitCorr$Familiar, FitCorr$Foreign, method = 'pearson')

# No correlation apparent here at end of evoution


#### Test if there is a link between how infectious virions are and how productive an infection is
### do the two components of fitness positively correlate?

# Do this by taking one of the models above and adding in the response variable from the other as a predictor

mixed(AveCount~ PlodiaAssayLine + Self + PropInf + VirusLine + (1|EvolvedLine/VirusLine) , 
      data = IFDEMD, family = poisson, method = 'LRT')
MMod6e<-glmer(AveCount~ PlodiaAssayLine + EvolvedLine + Self + PropInf + (1|EvolvedLine/VirusLine) , 
              data = IFDEMD, family = poisson)
summary(MMod6e)

# significant negative correlation / covariation between infectiousness and production at end of evolution and effects of self (positive), ELine 17,and ALines.
# note: >0 larvae infected for all assays in IFDEMD so not driven by 'matching 0s' in infected and fitness or whatever

########
## Let's now attempt to look across passages.
## much of the same analysis, but with a 'time point' variable too.
########

# We want to know if the effect of 'self' changes with passage number
# i.e. did we catch evolution in action?

# Treat passage as a repeated measure within each line and as a factor

# Look at max dose only
IFDMD <- IFD[which(IFD$Dose == max(IFD$Dose)),]

IFDMD$PassageNumber <- factor(as.character(IFDMD$PassageNumber), ordered = F)

mixed(VirFit ~ PlodiaAssayLine +  Self + EvolvedLine + PassageNumber + Self:PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
      data = IFDMD, family = poisson, method = 'LRT')

# Passage number:self interaction significant
# Fitness of lines changes with time
# specialism of lines changes with time.
# check directions & magnitudes

MMod5 <- glmer(VirFit ~ PlodiaAssayLine + EvolvedLine + Self*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
               data = IFDMD, family = poisson)

summary(MMod5)
# Note that fitness is highest at start, lowest at passage 6, and then increases again to approx. stabilise
# BUT effect of 'self' increases continually with passage (see Self:Passage interaction)
graphmod5 <- emmeans(MMod5, ~ Self *PassageNumber)
graphmod5 <-as.data.frame(graphmod5)
figure5<-ggplot(data=graphmod5, aes(x=PassageNumber, y=emmean, group=Self))+geom_smooth(aes(color=factor(Self)))+geom_ribbon(data=graphmod5, aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=.2)+
  theme_light()+labs(x="Passage Number", y="Relative Fitness", color="Familiar Host?")+
  scale_color_manual(values=sunset3)+scale_x_discrete(breaks=c(0,1,4,6,9))
# Look at infectivity over time
MMod5i <- glmer(cbind(NumberInfected, NumberUninfected) ~ PlodiaAssayLine + EvolvedLine + Self*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
                data = IFDMD, family = binomial)
graphmod5i <- emmeans(MMod5i, ~ Self *PassageNumber)
graphmod5i <-as.data.frame(graphmod5i)
figure5i<-ggplot(data=graphmod5i, aes(x=PassageNumber, y=emmean, group=Self))+geom_smooth(aes(color=factor(Self)))+geom_ribbon(data=graphmod5i, aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=.2)+
     theme_light()+labs(x="Passage Number", y="Relative Infectivity", color="Familiar Host?")+
    scale_color_manual(values=sunset3)+scale_x_discrete(breaks=c(0,1,4,6,9))

#look at productivity
MMod5p <- glmer(AveCount ~ PlodiaAssayLine + EvolvedLine + Self*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
                data = IFDMD, family = poisson)
graphmod5p <- emmeans(MMod5p, ~ Self *PassageNumber)
graphmod5p <-as.data.frame(graphmod5p)
figure5p<-ggplot(data=graphmod5p, aes(x=PassageNumber, y=emmean, group=Self))+geom_smooth(aes(color=factor(Self)))+geom_ribbon(data=graphmod5p, aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=.2)+
  theme_light()+labs(x="Passage Number", y="Relative Productivity", color="Familiar Host?")+
  scale_color_manual(values=sunset3)+scale_x_discrete(breaks=c(0,1,4,6,9))

#doing a model with evolved line in the interaction so that I can plot it
#MMod5a <- glmer(VirFit ~ PlodiaAssayLine + EvolvedLine*Self*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
               # data = IFDMD, family = poisson)
#grab out the model effects and convert to dataframe
#graphmod5a <- emmeans(MMod5a, ~ EvolvedLine * Self *PassageNumber)
#graphmod5a <-as.data.frame(graphmod5a)
#graphmod5a$LabEvolvedLine<- NA
#graphmod5a$LabEvolvedLine[which(graphmod5a$EvolvedLine == '2')]<-'Evolved on Host Genotype 2'
#graphmod5a$LabEvolvedLine[which(graphmod5a$EvolvedLine == '9')]<-'Evolved on Host Genotype 9'
#graphmod5a$LabEvolvedLine[which(graphmod5a$EvolvedLine == '17')]<-'Evolved on Host Genotype 17'

#plot the effect of self over time, split by evolved on line to see dynamics better
#figure5a<-ggplot(data=graphmod5a, aes(x=PassageNumber, y=emmean, group=Self))+geom_smooth(aes(color=factor(Self)))+geom_ribbon(data=graphmod5, aes(ymin=asymp.LCL, ymax=asymp.UCL), alpha=.2)+
  #theme_light()+labs(x="Passage Number", y="Relative Fitness", color="Familiar Host?")+
  #scale_color_manual(values=sunset3)+facet_wrap(~LabEvolvedLine)+scale_x_discrete(breaks=c(0,1,4,6,9))

figure5
figure5i
figure5p

#also make a graph split by virus line
MMod5b <- glmer(VirFit ~ PlodiaAssayLine + VirusLine*Self*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
                data = IFDMD, family = poisson)
graphmod5b <- emmeans(MMod5b, ~ VirusLine * Self *PassageNumber)
graphmod5b <-as.data.frame(graphmod5b)
figure5b<-ggplot(data=graphmod5b, aes(x=PassageNumber, y=emmean, group=Self))+geom_smooth(aes(color=factor(Self)))+
  theme_light()+labs(x="Passage Number", y="Relative Fitness", color="Familiar Host?")+
  scale_color_manual(values=sunset3)+facet_wrap(~VirusLine)+scale_x_discrete(breaks=c(0,1,4,6,9))

figure5

#let's check fitness over time w/o P0 since this had different storage/extraction conditions and could be different
IFDMD1<- IFDMD[which(IFDMD$PassageNumber != 0),] 
mixed(VirFit ~ PlodiaAssayLine +  Self + PassageNumber  + (1|EvolvedLine/VirusLine/PassageNumber), 
      data = IFDMD1, family = poisson, method = 'LRT')
MMod7 <- glmer(VirFit ~ PlodiaAssayLine + EvolvedLine + Self + PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber), 
               data = IFDMD1, family = poisson)
summary(MMod7)

# Fitness increases from P1->P4 (effect size + 1.25) and than seems to stay around there (so no change past pass 4). 
# 0.9 and 1.25 aren't so different that I'd be willing to say there's enough of a change

# Investigate if the correlation between infection likelihood and productivity at the end of experiment
# changes in magnitude (or even direction) with time
# e.g. might start positive, become negative

#check for corrolation across the whole data set
mixed(AveCount ~ PlodiaAssayLine + EvolvedLine + Self + PropInf + PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber) , 
      data = IFDMD, family = poisson, method = 'LRT')
# Yes significant correlation / covariation between infectiousness and production across whole dataset, also passage number, self, eline, aline
#check direction
MMod6 <- glmer(AveCount ~ PlodiaAssayLine + EvolvedLine + Self + PropInf + PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber) , 
               data = IFDMD, family = poisson)
summary(MMod6)
#Prop Inf is signifiantly positive across whole dataset
#So corrolation changes from pass 9 only and I see things on plot, now lets look for those interactions that I see on the plot

mixed(AveCount ~ PlodiaAssayLine + EvolvedLine + Self*PropInf*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber) , 
      data = IFDMD, family = poisson, method = 'LRT')
#Self:PropInf:PassageNumber is significant (and everything else)

MMod6b <- glmer(AveCount ~ PlodiaAssayLine + EvolvedLine + Self*PropInf*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber) , 
                data = IFDMD, family = poisson) 
summary(MMod6b)
##The SelfTRUE:AveCount:Pass interactions have significant increasingly negative effect

figure6<-ggplot(data=IFDMD, aes(x=PropInf, y=AveCount, group=PassageNumber)) + geom_point(aes(color=factor(PassageNumber)))+ 
  geom_smooth(method = "lm", aes(color=factor(PassageNumber)), alpha=.25) +theme_light()+labs(x="Viral Infectivity", y="Viral Productivity", color="Passage") +
  facet_wrap(~Self)+scale_color_manual(values=sunset5)
figure6


# Finally, let's look at the correlation between fitness on self=T and fitness on self=F across passages
# at p=9 we found no correlation, but let's look across the whole data as at p=0, there should have been none (by definition)

# use doses 4,6,9 as full data missing for p = 1 and for p = 0 there is no 'self=T' by definition.

head(IFDMD1)

IFDMD2 <- IFDMD1[which(as.numeric(as.character(IFDMD1$PassageNumber)) != 1),]

# IFDMD2 <- IFDMD1[which(as.numeric(as.character(IFDMD1$PassageNumber)) != 1 & as.numeric(as.character(IFDMD1$PassageNumber)) != 4 ),]


# build fitness model with passage number nested in and interacting with the evolution of specificity 
# see above for simpler P=9 equivalent
MMod4b <- glmer(VirFit ~ PlodiaAssayLine + VirusLine*Self*PassageNumber + (1|EvolvedLine/VirusLine/PassageNumber) , 
               data = IFDMD2, family = poisson)

# pairwise comparisons 'emmeans abuse' same as above but with an extra 'PassageLevel' variable
# we don't care about significance here just effect sizes so no correction. 
# simiarly, the warnings mainly relate to error estimates which isn't a huge concern

PWC3 <-contrast(emmeans(MMod4b, ~ VirusLine*Self*PassageNumber),
                interaction = TRUE, method = "pairwise", by = c("Self", "PassageNumber"))

PWC3

PWC4 <- contrast(emmeans(MMod4b, ~ VirusLine*Self*PassageNumber),
                 interaction = TRUE, method = "pairwise", by = c("Self", "VirusLine"))
PWC4


# Note this is *extremely brittle* from hard-coding
# Check structure of PWC3 for this to make sense, 36 pairwise comparisons per passage+self? combo
# PWC4 sets the change in the standard, we fix line 2.1 in passage 4 as the '0' of the relative fitness
# 

FitCorr <- data.frame(Line = rep.int(unique(IFDMD2$VirusLine), times = NROW(unique(IFDMD2$PassageNumber))))

FitCorr$Pass <- sort(rep.int(unique(IFDMD2$PassageNumber), times = NROW(levels(IFDMD2$VirusLine))))

FitCorr$Foreign <- append(x = 0, 
                          values = c(
                            as.data.frame(PWC3)$estimate[1:8],
                            as.data.frame(PWC4)$estimate[1],
                            as.data.frame(PWC3)$estimate[73:80],
                            as.data.frame(PWC4)$estimate[2],
                            as.data.frame(PWC3)$estimate[145:152]
                            ), 
                                            after = 1)
FitCorr$Familiar <- append(x = 0, 
                          values = c(
                            as.data.frame(PWC3)$estimate[37:44],
                            as.data.frame(PWC4)$estimate[28],
                            as.data.frame(PWC3)$estimate[109:116],
                            as.data.frame(PWC4)$estimate[29],
                            as.data.frame(PWC3)$estimate[181:188]
                          ), 
                          after = 1)

## make an evolved on column
FitCorr$EvolvedLine<-NA
FitCorr$EvolvedLine[which(FitCorr$Line=='2.1')]<-'2'
FitCorr$EvolvedLine[which(FitCorr$Line=='2.2')]<-'2'
FitCorr$EvolvedLine[which(FitCorr$Line=='2.3')]<-'2'
FitCorr$EvolvedLine[which(FitCorr$Line=='9.1')]<-'9'
FitCorr$EvolvedLine[which(FitCorr$Line=='9.2')]<-'9'
FitCorr$EvolvedLine[which(FitCorr$Line=='9.3')]<-'9'
FitCorr$EvolvedLine[which(FitCorr$Line=='17.1')]<-'17'
FitCorr$EvolvedLine[which(FitCorr$Line=='17.2')]<-'17'
FitCorr$EvolvedLine[which(FitCorr$Line=='17.3')]<-'17'

figure7<-ggplot(data=FitCorr, aes(x=Foreign, y=Familiar, group=EvolvedLine)) + geom_point(aes(color=factor(EvolvedLine), alpha=Pass))+ 
  geom_smooth(method = "lm", aes(color=factor(EvolvedLine)), se=FALSE) +theme_light()+labs(x="Fitness on Foreign Hosts", y="Fitness on Familiar Hosts", color="Virus Evolved On Host Genotype", alpha="Passage Number") +
  scale_color_manual(values=sunset3)+scale_alpha_manual(values=c(.4,.7,1))

cor.test(FitCorr$Familiar, FitCorr$Foreign, method = 'pearson')
cor.test(FitCorr$Familiar, FitCorr$Foreign, method = 'spearman')

# with 2 outliers removed

FitCorr2 <- FitCorr[c(-7,-9),]

cor.test(FitCorr2$Familiar, FitCorr2$Foreign, method = 'pearson')
cor.test(FitCorr2$Familiar, FitCorr2$Foreign, method = 'spearman')

# nothing in either case
