#to view any forest plots or model summary information, uncomment the correspondingly labeled code
#to find them more easily search for "uncomment for", then uncomment only the beginning of the corresponding lines
#set working directory to the folder codesheet2.csv is in

require(ggplot2)
require(Rmisc)
require(meta)
require(metafor)

ds<-read.csv("codesheet2.csv", na.strings = c(999, "NA", ""))

#convert ES to number from factor
ds$ES<-as.numeric(ds$ES)


ds$FirstAuthor<-as.character(ds$FirstAuthor)
ds$Year<-as.character(ds$Year)

#change title
ds$FirstAuthor<- paste(ds$FirstAuthor, "et al. (", sep = " ")
ds$FirstAuthor<- paste(ds$FirstAuthor, ds$Year, sep = "")
ds$FirstAuthor<- paste(ds$FirstAuthor, ")", sep = "")
ds$Year<-NULL

########make study numbers
rown<-1
studyn<-1
for (rown in 1:length(ds$Link) ){
  if (rown == 1){
    ds$studyn[rown]<- studyn
    studyn <- studyn + 1
  }
  else if (ds$Link[rown] == ds$Link[rown - 1]){
    ds$studyn[rown] = studyn - 1
  }
  else {
    ds$studyn[rown]<- studyn
    studyn<- studyn + 1
  }
}
###################################################
##BEGIN ACTUAL MATH ANALYSES######################
##################################################

#calculate SE from CI if missing and possible
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
  if (!is.na(ds$ShamCIUpper[rown]) && is.na(ds$ShamSD[rown]) ){ #if we're missing SE in the row but have CI
    ds$ShamSE[rown] <- sqrt(ds$ShamN[rown])*( (ds$ShamCIUpper[rown] - ds$ShamMean[rown]) /1.96) #calculate Sham SE
    ds$ActiveStimSE[rown] <- sqrt(ds$ActiveStimN[rown])*( (ds$ActiveStimCIUpper[rown] - ds$ActiveStimMean[rown]) /1.96) #calc Stim SE
  } 
}


#calculate SD from SE if missing and possible
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
  if (!is.na(ds$ShamSE[rown]) && is.na(ds$ShamSD[rown]) ){ #if we're missing SD in the row but have SE
    ds$ShamSD[rown] <- ds$ShamSE[rown]*(sqrt(ds$ShamN[rown] - 1)) #calculate Sham SD
    ds$ActiveStimSD[rown] <- ds$ActiveStimSE[rown]*(sqrt(ds$ActiveStimN[rown])) #calc Stim SD
    ds$CalcSD[rown]<-TRUE
  } 
  else if (!is.na(ds$PooledSE[rown]) && is.na(ds$PooledSD[rown]) ){ #if we're missing pooled SD but have pooled SE
    ds$PooledSD[rown]<- ds$PooledSE[rown] / (sqrt((1/ds$ShamN[rown]) + (1/ds$ActiveStimN[rown] ) ) ) #calculate
    ds$CalcSD[rown]<-FALSE
  }
  else{
    ds$CalcSD[rown]<-FALSE
  }
}


#calculate t value from means and SDs if missing and possible (only for independent t tests)
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
  if (is.na(ds$t[rown] && !is.na(ds$ShamMean[rown]) && !is.na(ds$ActiveStimMean[rown]) 
      && !is.na(ds$ShamSD[rown]) && !is.na(ds$ActiveStimSD[rown]) ) ) {# if we don't have t but we have means and SDs
    if (ds$Design[rown] == "B"){ #for between subjects designs
      x1<-ds$ActiveStimMean[rown]
      x2<-ds$ShamMean[rown]
      var1<-(ds$ActiveStimSD[rown])^2
      var2<-(ds$ShamSD[rown])^2
      n1<-ds$ActiveStimN[rown]
      n2<-ds$ShamN[rown]
      ds$t[rown]<- (x1-x2)/sqrt( ( ((var1*(n1-1)) + (var2*(n2-1)) ) /(n1+n2-2) )*((1/n1) + (1/n2)) )
      ds$CalcTfromMean[rown]<-TRUE
    }
    else{
      ds$CalcTfromMean[rown]<-FALSE
    }
  }
  else{
    ds$CalcTfromMean[rown]<-FALSE
  }
}

#calculate t value from F value if missing and possible
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
  if (is.na(ds$t[rown]) && !is.na(ds$Fval[rown]) && ds$Fdf[rown]==1 ){ #if we have F, df is 1 and we don't have t
    ds$t[rown]<- ds$Fval[rown]^2
    ds$CalctFromF[rown]<-TRUE
  }
  else{
    ds$CalctFromF[rown]<-FALSE
  }
}

#calculate p value from t value if missing and possible
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
  if (is.na(ds$Exactp[rown]) && !is.na(ds$t[rown]) ){ #if we have t but not p
    #calculate df, different for between and within designs
    if (ds$Design[rown] == "B"){
      df = (ds$TotalN[rown] - 2)
    }
    else if (ds$Design[rown] == "W"){
      df = (ds$TotalN[rown] - 1)
    }
    #calculate p, two tailed t test
    ds$Exactp[rown]<- 2*pt(abs(ds$t[rown]), df, lower.tail = FALSE)
    ds$Calcp[rown]<- TRUE
  }
  else{
    ds$Calcp[rown]<- FALSE
  }
}


#calculate t value from p value if missing and possible
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
 if (!is.na(ds$Exactp[rown] && is.na(ds$t[rown]) ) ){ #if we have p but not t
   #calculate df depending on design, between or within
   if (ds$Design[rown] == "B"){
     df = ds$TotalN[rown] - 2
   }
   else if (ds$Design[rown] == "W"){
     df = ds$TotalN[rown] - 1
   }
   ds$t[rown]<-qt((ds$Exactp[rown]/2), df)
   ds$CalctFromp[rown]<-TRUE
 }
  else{
    ds$CalctFromp[rown]<-FALSE
  }
}  
   
##########################################################################################################
#############EFFECT SIZE CALCULATIONS#####################################################################
##########################################################################################################

ds$usable<-NA

#determine what is usable
for (rown in 1:length(ds$Link) ){ #for each row, where rown is row number
  if (ds$Design[rown]=="B" && 
      !is.na(ds$Exactp[rown])) { #if between and we have p, which means we have means and Sds or t
    ds$usable[rown]<-1 #can use it
  }
  else if (ds$Design[rown]=="W" && 
           !is.na(ds$ShamMean[rown]) && 
           !is.na(ds$ShamSD[rown])){#if within and we have means and SDs
    ds$usable[rown]<-1
  }
  else if (ds$Design[rown] == "W" &&
           !is.na(ds$t[rown])){#if within and we have t (debatable)
    ds$usable[rown]<-2 #maybe use?
  }
  else{ #if we don't have a p value
    ds$usable[rown]<-0 #cant use it
  }
} 

use<-ds[ds$usable != 0,] #this is the dataset with usable studies, including the "maybe use" (within subj for which all we have is F or t and p)
missing<-ds[ds$usable == 0,] #this is the dataset with missing numbers


#calculate pooled SD for between subj (for within we'll be using sham sd to calculate es)
for (rown in 1:length(use$Link) ){ #for each row, where rown is row number
  if (is.na(use$PooledSD[rown]) && !is.na(use$ShamSD[rown]) && !is.na(use$ActiveStimSD[rown] && use$Design[rown]=="B")){ #if we don't have it already but have both SDs
    use$PooledSD[rown]<- (sqrt((((use$ActiveStimN[rown]-1)*(use$ActiveStimSD[rown]^2)) + ((use$ShamN[rown]-1)*(use$ShamSD[rown]^2)))/(use$ActiveStimN[rown] + use$ShamN[rown] - 2)))
  }
}

#calculate mean diff if we don't have it and we can
for (rown in 1:length(use$Link) ){ #for each row, where rown is row number
 if (is.na(use$MeanDiff[rown]) && !is.na(use$ShamMean[rown]) && !is.na(use$ActiveStimMean[rown])){ #if we don't have meandiff already but we have means
   use$MeanDiff[rown]<-use$ActiveStimMean[rown] - use$ShamMean[rown]
 } 
}  

#calculate ES directly if possible
for (rown in 1:length(use$Link) ){ #for each row, where rown is row number
 if (use$Design[rown] == "B"){ #for between subjects studies
   if (is.na(use$ES[rown]) && !is.na(use$MeanDiff[rown]) && !is.na(use$PooledSD[rown]) ){ #if ES wasn't reported but we can calculate it directly
     use$ES[rown]<- use$MeanDiff[rown]/use$PooledSD[rown]
     use$calcESdirect[rown]<-TRUE
   }
   else{
     use$calcESdirect[rown]<-FALSE
   }
 }#end if between
  if (use$Design[rown] == "W"){#for within subjects studies
    if (is.na(use$ES[rown]) && !is.na(use$MeanDiff[rown]) && !is.na(use$ShamSD[rown]) ){ #if ES wasn't reported but we can calculate it directly
      use$ES[rown]<- use$MeanDiff[rown]/use$ShamSD[rown]
      use$calcESdirect[rown]<-TRUE
    }
    else{
      use$calcESdirect[rown]<-FALSE
    }
  }#end if within
}#end for loop


#calculate ES from t if possible and we couldn't calculate it directly (all F's have been converted to t's)
for (rown in 1:length(use$Link) ){ #for each row, where rown is row number
  
  #use line below if only for independent t tests
  #if (use$Design[rown] == "B" && is.na(use$ES[rown]) && !is.na(use$t[rown])){#if we don't have ES but we have t
  
  #comment out line below if only for independent t tests
  if (is.na(use$ES[rown]) && !is.na(use$t[rown])){#if we don't have ES but we have t
    use$ES[rown] <- use$t[rown]*sqrt((use$ShamN[rown] + use$ActiveStimN[rown])/(use$ShamN[rown]*use$ActiveStimN[rown]))
    use$calcESfromt[rown]<-TRUE
  }
  else{
    use$calcESfromt[rown]<-FALSE
  }
}



######################################################################################################
#########EFFECT SIZE ADJUSTMENTS######################################################################
######################################################################################################


#Hedges g correction
use$g<-use$ES*(1 - (3/(4*use$TotalN - 9)))


#calculate ES SE
use$g.SE<- sqrt( ( (use$ActiveStimN + use$ShamN) / (use$ActiveStimN*use$ShamN) ) + (use$ES^2/( 2*(use$ActiveStimN + use$ShamN) ) ) )


#account for direction
for (rown in 1:length(use$Link) ){ #for each row, where rown is row number
  if (use$ESDir[rown] == 1){#if ES should be positive
    if (use$g[rown] < 0){#and it's negative
      use$g[rown] <- use$g[rown] * -1
    }
  }
  else if (use$ESDir[rown] == -1){#if ES should be negative
    if (use$g[rown] > 0){#and it's positive
      use$g[rown] <- use$g[rown] * -1
    }
  }
}


###############################################################################################
#####FOREST PLOTS AND MODERATOR ANALYSES#######################################################
###############################################################################################


#sort by N for cumulative plot if desired
use<-use[order(-use$TotalN), ]

#Overall effect size models, split by stimulation type
anodun<-use[use$StimType == "A" & use$StimLaterality == "U",] #unilateral anodal studies
cathun<-use[use$StimType == "C" & use$StimLaterality == "U",] #unilateral cathodal studies
bilat<-use[use$StimLaterality == "B",] #bilateral studies
ares<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor) #total model for unilateral anodal studies
cres<-rma(yi = cathun$g, vi = cathun$g.SE, data = cathun, method = "DL", slab = cathun$FirstAuthor) #total model for unilateral cathodal studies
bres<-rma(yi = bilat$g, vi = bilat$g.SE, data = bilat, method = "DL", slab = bilat$FirstAuthor) #total model for bilateral studies
#forest(ares) #uncomment for anodal unilateral plot
#title("Anodal Unilateral Stimulation") #uncomment for anodal unilateral plot
#forest(cres) #uncomment for cathodal unilateral plot
#title("Cathodal Unilateral Stimulation") #uncomment for cathodal unilateral plot
#forest(bres) #uncomment for bilateral plot
#title("Bilateral Stimulation") #uncomment for bilateral plot


#moderator analysis for domain
adomains<-adomains[adomains$Domain != "G",] #removes the few "Gen EF" tasks that did not fit into a single domain
aresdomains<-rma(yi = adomains$g, vi = adomains$g.SE, data = adomains, method = "DL", slab = adomains$FirstAuthor, mods = ~ adomains$Domain) #moderator model for domain
uanod<-anodun[anodun$Domain == "U",] #unilateral anodal updating studies
ianod<-anodun[anodun$Domain == "I",] #unilateral anodal inhibition studies
sanod<-anodun[anodun$Domain == "S",] #unilateral anodal shifting studies
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



#moderator analysis for cathloc

#dummy code for cathloc, cranial is reference group
anodun<-anodun[anodun$CathLoc > 2,] #removes the few multiple cathode montages
anodun$c1<-0
for (rown in 1:length(anodun$Link) ){ #for each row, where rown is row number
  if (anodun$CathLoc[rown] == 8){#deltoid, 7 studies
    anodun$c1[rown]<- 1
  }
  else if (anodun$CathLoc[rown] == 10){#mastoid, 3 studies
    anodun$c1[rown]<- 1
  }
}
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


#moderator analysis for anode size 
areselec<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$FocalSize) #model for anode size moderator
alarge<-anodun[anodun$FocalSize == 35,] #large anode studies
asmall<-anodun[anodun$FocalSize < 30,] #small anode studies
alargeres<-rma(yi = alarge$g, vi = alarge$g.SE, data = alarge, method = "DL", slab = alarge$FirstAuthor) #large anode subgroup model
asmallres<-rma(yi = asmall$g, vi = asmall$g.SE, data = asmall, method = "DL", slab = asmall$FirstAuthor) #small anode subgorup model
#summary(areselec) #uncomment for anode size moderator model
#forest(alargeres) #uncomment for large anode subgroup plot
#title("Anodal Unilateral, Anode 35 cm sq") #uncomment for large anode subgroup plot
#forest(asmallres) #uncomment for small anode subgroup plot
#title("Anodal Unilateral, Anode 25 cm sq or less") #uncomment for small anode subgroup plot

#moderator analysis for stim intensity 
aresint<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = factor(anodun$StimIntensity)) #stim intensity moderator model
#summary(aresint) #uncomment for stim intensity moderator summary

#moderator analysis for age
aresage<-rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$MeanAge) #age moderator model
#summary(aresage) #uncomment for age moderator summary


#moderator analysis for handedness  
ahand<-anodun[!is.na(anodun$Handedness),] #remove studes where handedness could not be determined
ahand$h1<-0 #dummy code for handedness, 0 is not controlled for and 1 is only right
ahand$h1[ahand$Handedness == "R"]<-1 #dummy code for handedness, 0 is not controlled for and 1 is only right
ahandres<- rma(yi = ahand$g, vi = ahand$g.SE, data = ahand, method = "DL", slab = ahand$FirstAuthor, mods = ahand$h1) #moderator model for handedness
#summary(ahandres) #uncomment for handedness model summary


#moderator analysis for gender
agen<-anodun[!is.na(anodun$Gender),]
agenres<- rma(yi = agen$g, vi = agen$g.SE, data = agen, method = "DL", slab = agen$FirstAuthor, mods = ~agen$Gender)
#summary(agenres) #uncomment for gender model summary

#moderator analysis for design (within or between)
anodun$d1<-0 #create dummy code column for design, 0 is between and 1 is within
anodun$d1[anodun$Design == "W"]<-1 #create dummy code column for design, 0 is between and 1 is within
adesres<- rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$d1)
#summary(adesres) #uncomment for design model summary

#moderator analysis for current density
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

#anode location moderator analysis
anodun$a1<-0 #create dummy code column for left vs right dlpfc anode location
anodun$a1[anodun$AnLoc == 2]<-1 #create dummy code column for left vs right dlpfc anode location
aresloc<- rma(yi = anodun$g, vi = anodun$g.SE, data = anodun, method = "DL", slab = anodun$FirstAuthor, mods = anodun$a1) #moderator model for anode location
#summary(aresloc) #uncomment for anode location moderator model summary


#MODERATOR ANALYSIS FOR ONLINE VS OFFLINE
aoo<-anodun[anodun$Online.Offline<9,] #remove studies that did a mixture of online and offline, coded as 9
aoo$o1<-0 #dummy code column for online/offline, 1 is online and 0 is offline
aoo$o1[aoo$Online.Offline == 2]<-1 #dummy code column for online/offline, 1 is online and 0 is offline
aresoo<-rma(yi = aoo$g, vi = aoo$g.SE, data = aoo, method = "DL", slab = aoo$FirstAuthor, mods = aoo$o1) #online/offline moderator analysis model
#summary(aresoo) #uncomment for online/offline moderator analysis summary