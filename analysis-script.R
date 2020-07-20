#to view any forest plots or model summary information, uncomment the correspondingly labeled code
#to find them more easily search for "uncomment for", then uncomment only the beginning of the corresponding lines
#set working directory to the folder codesheet2.csv is in

require(ggplot2)
require(Rmisc)
require(meta)
require(metafor)

ds<-read.csv("codesheet2.csv", na.strings = c(999, "NA", ""))


#convert from author and year to character 
#this is to manipulate them in the next section
ds$FirstAuthor<-as.character(ds$FirstAuthor)
ds$Year<-as.character(ds$Year)

#change FirstAuthor to include Year
#this column now serves as a general purpose paper title
#remove Year column as the info is now redundant
ds$FirstAuthor<- paste(ds$FirstAuthor, "et al. (", sep = " ")
ds$FirstAuthor<- paste(ds$FirstAuthor, ds$Year, sep = "")
ds$FirstAuthor<- paste(ds$FirstAuthor, ")", sep = "")
ds$Year<-NULL

#################################################################
##IMPUTE STATS TO BE USED FOR EFFECT SIZE CALCULATIONS###########
#################################################################

#identify cases where SE is missing but we can impute it from CIs
calcSE<-which(is.na(ds$ShamSE) & !is.na(ds$ShamCIUpper))
#calc SE from CI
ds$ShamSE[calcSE]<-
  sqrt(ds$ShamN[calcSE])*( (ds$ShamCIUpper[calcSE] - ds$ShamMean[calcSE]) /1.96)
ds$ActiveStimSE[calcSE]<-
  sqrt(ds$ActiveStimN[calcSE])*( (ds$ActiveStimCIUpper[calcSE] - ds$ActiveStimMean[calcSE]) /1.96)

#identify cases where SD is missing but SE is available
calcSD<-which(!is.na(ds$ShamSE) & is.na(ds$ShamSD))
#calc SD from SE in case we need it for effect size calc
ds$ShamSD[calcSD] <- 
  ds$ShamSE[calcSD]*(sqrt(ds$ShamN[calcSD] - 1))
ds$ActiveStimSD[calcSD] <- 
  ds$ActiveStimSE[calcSD]*(sqrt(ds$ActiveStimN[calcSD]))

#identify cases where we are missing t but can calculate it from F
calct<-which(is.na(ds$t) & ds$Fdf==1)
#calc t from F
ds$t[calct]<-ds$Fval[calct]^2

#identify cases where pooled SD is missing but pooled SE is available
calcpSD_frompSE<-which(!is.na(ds$PooledSE) & is.na(ds$PooledSD))
ds$PooledSD[calcpSD_frompSE]<- ds$PooledSE[calcpSD_frompSE] / 
  (sqrt((1/ds$ShamN[calcpSD_frompSE]) + (1/ds$ActiveStimN[calcpSD_frompSE] ) ) ) #calculate


#identify between subj studies for which we need to calculate pooled SD
#(for within subj we'll be using sham sd to calculate es)
calc_pSD<-which(is.na(ds$PooledSD) & ds$Design=="B" & #missing pooled and between subj
                  !(is.na(ds$ShamSD) & !(is.na(ds$ActiveStimSD)))) #have SDs for both groups
#calculate pooled SD for between subj 
ds$PooledSD[calc_pSD]<- (sqrt((((ds$ActiveStimN[calc_pSD]-1)*(ds$ActiveStimSD[calc_pSD]^2)) + #calculation
                               ((ds$ShamN[calc_pSD]-1)*(ds$ShamSD[calc_pSD]^2)))/
                              (ds$ActiveStimN[calc_pSD] + ds$ShamN[calc_pSD] - 2)))


#identify studies for which we need to calculate mean difference
calcMdiff<-which(is.na(ds$MeanDiff) & !is.na(ds$ShamMean) & !is.na(ds$ActiveStimMean))
ds$MeanDiff[calcMdiff]<-ds$ActiveStimMean[calcMdiff] - ds$ShamMean[calcMdiff] #diff calc


#################################################################
#############EFFECT SIZE CALCULATIONS############################
#################################################################


#calculate ES directly if possible
missingES_B<-which(is.na(ds$ES) & ds$Design=="B")
ds$ES[missingES_B]<- ds$MeanDiff[missingES_B]/ds$PooledSD[missingES_B]

missingES_W<-which(is.na(ds$ES) & ds$Design=="W")
ds$ES[missingES_W]<- ds$MeanDiff[missingES_W]/ds$ShamSD[missingES_W]

#calculate ES from t if we couldn't do it directly
missingES<-which(is.na(ds$ES))

ds$ES[missingES] <- ds$t[missingES]*sqrt((ds$ShamN[missingES] + ds$ActiveStimN[missingES])/
                                           (ds$ShamN[missingES]*ds$ActiveStimN[missingES]))


####################################################################
#########REMOVE STUDIES FOR WHICH WE COULDN'T CALCULATE ES##########
####################################################################

use<-subset(ds, !is.na(ES))
missing<-subset(ds, is.na(ES))

##########################################
#########EFFECT SIZE ADJUSTMENTS##########
##########################################


#Hedges g correction
use$g<-use$ES*(1 - (3/(4*use$TotalN - 9)))


#calculate ES SE
use$g.SE<- sqrt( ( (use$ActiveStimN + use$ShamN) / (use$ActiveStimN*use$ShamN) ) + (use$ES^2/( 2*(use$ActiveStimN + use$ShamN) ) ) )


#account for direction
wrongdirES<-which( (use$ESDir== -1 & use$g > 0) | #should be negative but is positive
                     (use$ESDir== 1 & use$g < 0) ) #should be positive but is negative

#################################################
#####FOREST PLOTS AND MODERATOR ANALYSES#########
#################################################


#sort by N 
use<-use[order(-use$TotalN), ]

#Overall effect size models, split by stimulation type  ####
anodun<-subset(use, StimType == "A" & StimLaterality == "U") #unilateral anodal studies
cathun<-subset(use, StimType == "C" & StimLaterality == "U") #unilateral cathodal studies
bilat<-subset(use, StimLaterality == "B") #bilateral studies
ares<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor) #total model for unilateral anodal studies
cres<-rma(yi = cathun$g, vi = cathun$g.SE, data = cathun, method = "DL", slab = cathun$FirstAuthor) #total model for unilateral cathodal studies
bres<-rma(yi = bilat$g, vi = bilat$g.SE, data = bilat, method = "DL", slab = bilat$FirstAuthor) #total model for bilateral studies
#forest(ares) #uncomment for anodal unilateral plot
#title("Anodal Unilateral Stimulation") #uncomment for anodal unilateral plot
#forest(cres) #uncomment for cathodal unilateral plot
#title("Cathodal Unilateral Stimulation") #uncomment for cathodal unilateral plot
#forest(bres) #uncomment for bilateral plot
#title("Bilateral Stimulation") #uncomment for bilateral plot


#moderator analysis for domain  ####
adomains<-adomains[adomains$Domain != "G",] #removes the few "Gen EF" tasks that did not fit into a single domain
aresdomains<-rma(yi = adomains$g, vi = adomains$g.SE, data = adomains, method = "DL", slab = adomains$FirstAuthor, mods = ~ adomains$Domain) #moderator model for domain
uanod<-subset(anodun, Domain == "U") #unilateral anodal updating studies
ianod<-subset(anodun,Domain == "I") #unilateral anodal inhibition studies
sanod<-subset(anodun,Domain == "S") #unilateral anodal shifting studies
ures<-rma(yi = uanod$g, vi = uanod$g.SE, data = uanod, method = "DL", slab = uanod$FirstAuthor) #subgroup analysis for updating AU studies
ires<- rma(yi = ianod$g, vi = ianod$g.SE, data = ianod, method = "DL", slab = ianod$FirstAuthor) #subgroup analysis for inhibition AU studies
sres<- rma(yi = sanod$g, vi = sanod$g.SE, data = sanod, method = "DL", slab = sanod$FirstAuthor) #subgroup analysis for shifting AU studies
#summary(aresdomains) #uncomment for domain moderator model
#forest(aresdomains) #uncomment for domain moderator plot
#title ("Moderator Analysis for Domain, Unilateral Anodal") #uncomment for domain moderator plot
#forest(ures) #uncomment for updating subgroup plot
#title("Anodal Unilateral, Updating Tasks") #uncomment for updating subgroup plot
#forest(ires) #uncomment for inhibition subgroup plot
#title("Anodal Unilateral, Inhibition Tasks") #uncomment for inhibition subgroup plot
#forest(sres) #uncomment for shifting subgroup plot
#title("Anodal Unilateral, Set Shifting Tasks") #uncomment for shifting subgroup plot



#moderator analysis for cathloc ####

#dummy code for cathloc, cranial is reference group
anodun<-anodun[anodun$CathLoc > 2,] #removes the few multiple cathode montages
#create a dummy code variable called c1
anodun$c1<-0
#make c1 equal to 1 for extracranial cathodes, stays 0 for cranial cathodes
anodun$c1[anodun$CathLoc==8]<-1 #8 meant deltoid
anodun$c1[anodun$CathLoc==10]<-1 #10 meant mastoid

arescathmod<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$c1) #cathode location moderator model
aub<-anodun[anodun$c1 == 0,] #cranial cathode studies
aunb<-anodun[anodun$c1 == 1,]#extracranial cathode studies
aubres<-rma(yi = aub$g, vi = aub$g.SE, data = aub, method = "DL", slab = aub$FirstAuthor) #cranial cathode subgroup model
aunbres<-rma(yi = aunb$g, vi = aunb$g.SE, data = aunb, method = "DL", slab = aunb$FirstAuthor) #extracranial cathode subgroup model
#summary(arescathmod) #uncomment for cathode location moderator model
#forest(arescathmod) #uncomment for cathode location moderator plot
#title("Anodal Unilateral, Cath Location as a moderator") #uncomment for cathode location moderator plot
#forest(aubres) #uncomment for cranial cathode subgroup plot
#title("Anodal Unilateral, Cranial Cathode") #uncomment for cranial cathode subgroup plot
#forest(aunbres) #uncomment for extracranial cathode subgroup plot
#title("Anodal Unilateral, Extracranial Cathode") #uncomment for extracranial cathode subgroup plot


#moderator analysis for anode size ####
areselec<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$FocalSize) #model for anode size moderator
alarge<-subset(anodun,FocalSize == 35) #large anode studies
asmall<-subset(anodun,FocalSize < 30) #small anode studies
alargeres<-rma(yi = alarge$g, vi = alarge$g.SE, data = alarge, method = "DL", slab = alarge$FirstAuthor) #large anode subgroup model
asmallres<-rma(yi = asmall$g, vi = asmall$g.SE, data = asmall, method = "DL", slab = asmall$FirstAuthor) #small anode subgorup model
#summary(areselec) #uncomment for anode size moderator model
#forest(alargeres) #uncomment for large anode subgroup plot
#title("Anodal Unilateral, Anode 35 cm sq") #uncomment for large anode subgroup plot
#forest(asmallres) #uncomment for small anode subgroup plot
#title("Anodal Unilateral, Anode 25 cm sq or less") #uncomment for small anode subgroup plot

#moderator analysis for stim intensity ####
aresint<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = factor(anodun$StimIntensity)) #stim intensity moderator model
#summary(aresint) #uncomment for stim intensity moderator summary

#moderator analysis for age ####
aresage<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$MeanAge) #age moderator model
#summary(aresage) #uncomment for age moderator summary


#moderator analysis for handedness  ####
ahand<-anodun[!is.na(anodun$Handedness),] #remove studes where handedness could not be determined
ahand$h1<-0 #dummy code for handedness, 0 is not controlled for and 1 is only right
ahand$h1[ahand$Handedness == "R"]<-1 #dummy code for handedness, 0 is not controlled for and 1 is only right
ahandres<- rma(yi = ahand$g, vi = ahand$g.SE, data = ahand, method = "DL", slab = ahand$FirstAuthor, mods = ahand$h1) #moderator model for handedness
#summary(ahandres) #uncomment for handedness model summary


#moderator analysis for gender ####
agen<-anodun[!is.na(anodun$Gender),]
agenres<- rma(yi = agen$g, vi = agen$g.SE, data = agen, method = "DL", slab = agen$FirstAuthor, mods = ~agen$Gender)
#summary(agenres) #uncomment for gender model summary

#moderator analysis for design (within or between) ####
anodun$d1<-0 #create dummy code column for design, 0 is between and 1 is within
anodun$d1[anodun$Design == "W"]<-1 #create dummy code column for design, 0 is between and 1 is within
adesres<- rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$d1)
#summary(adesres) #uncomment for design model summary

#moderator analysis for current density ####
anodun$density<-anodun$StimIntensity/anodun$FocalSize #create current density column
aresdens<- rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$density) #current density moderator model
#create rows to denote high, med and low density
anodun$hdens<-0
anodun$mdens<-1
anodun$ldens<-0
anodun$hdens[anodun$density >= quantile(anodun$density, 0.75)] <- 1
anodun$mdens[anodun$density >= quantile(anodun$density, 0.75)] <- 0
anodun$ldens[anodun$density <= quantile(anodun$density, 0.25)] <- 1
anodun$mdens[anodun$density <= quantile(anodun$density, 0.25)] <- 0
anodunld<-anodun[anodun$ldens == 1,] #low density studies
anodunmd<-anodun[anodun$mdens == 1,] #med density studies
anodunhd<-anodun[anodun$hdens == 1,] #high density studies
areshdens<- rma(yi = anodunhd$g, vi = anodunhd$g.SE, data = anodunhd, method = "DL", slab = anodunhd$FirstAuthor) #high density subgroup model
aresmdens<- rma(yi = anodunmd$g, vi = anodunmd$g.SE, data = anodunmd, method = "DL", slab = anodunmd$FirstAuthor) #med density subgroup model
aresldens<- rma(yi = anodunld$g, vi = anodunld$g.SE, data = anodunld, method = "DL", slab = anodunld$FirstAuthor) #low density subgroup model
#summary(aresdens) #uncomment for current density model summary
#forest(areshdens) #uncomment for high density subgroup plot
#title("Anodal Unilateral, High Density") #uncomment for high density subgroup plot
#forest(aresmdens) #uncomment for med density subgroup plot
#title("Anodal Unilateral, Medium Density") #uncomment for med density subgroup plot
#forest(aresldens) #uncomment for low density subgroup plot
#title("Anodal Unilateral, Low Density") #uncomment for low density subgroup plot

#Duration moderator analysis
aresdur<- rma(yi = aon$g, vi = aon$g.SE, data = aon, method = "DL", slab = aon$FirstAuthor, mods = aon$StimDur) #model for duration moderator
#summary(aresdur) #uncomment for duration moderator model summary

#anode location moderator analysis ####
anodun$a1<-0 #create dummy code column for left vs right dlpfc anode location
anodun$a1[anodun$AnLoc == 2]<-1 #create dummy code column for left vs right dlpfc anode location
aresloc<- rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$a1) #moderator model for anode location
#summary(aresloc) #uncomment for anode location moderator model summary


#moderator analysis for online vs offline ####
aoo<-anodun[anodun$Online.Offline<9,] #remove studies that did a mixture of online and offline, coded as 9
aoo$o1<-0 #dummy code column for online/offline, 1 is online and 0 is offline
aoo$o1[aoo$Online.Offline == 2]<-1 #dummy code column for online/offline, 1 is online and 0 is offline
aresoo<-rma(yi = aoo$g, vi = aoo$g.SE, data = aoo, method = "DL", slab = aoo$FirstAuthor, mods = aoo$o1) #online/offline moderator analysis model
#summary(aresoo) #uncomment for online/offline moderator analysis summary