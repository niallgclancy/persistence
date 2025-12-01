library(tidyverse)


#============================SECTION A: Net Species Changes to Native and Introduced species (if native anywhere in study area, restricted to native range)
#==================================================================================
#---------------------Net Species Changes to Only Native Populations
NET=read.csv("PersistenceMetrics.csv")
nogreen=read.csv("NOGREENmetrics.csv")
NOGREENsites=nogreen$RepetID


removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
NET=left_join(NET,removes,by="RepeatID")
NET=subset(NET,NET$GearMiss!="Y")
NET=subset(NET, NET$F_MAUG_<1000)
NET=subset(NET, NET$RepeatID!=234 & NET$RepeatID!=74& NET$RepeatID!=237& NET$RepeatID!=238)

NET=subset(NET,NET$RepeatID%in%NOGREENsites)

NET=NET%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")

####FROM HERE, MUST RUN "DatasetPrep.R" for section
NET=subset(NET,!is.na(NET$change))
NET$Species[which(NET$Species=="RMCOT" | NET$Species=="COLCOT")]="MOTCOT"
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
NET=left_join(NET,traits,by="Species")
NET$Status[which(NET$Species=="RDSH")]="Native"
NET$Status[which(NET$Species=="RDSH" & NET$RepeatID %in% rdshnn)]="Introduced"
NET$Status[which(NET$Species=="LKCH" & NET$RepeatID %in% lkchnn)]="Introduced"
NET$Status[which(NET$Species=="FHMN" & NET$RepeatID %in% FHMNnn)]="Introduced"
NET$Status[which(NET$Species=="LNDC" & NET$RepeatID %in% LNDCnn)]="Introduced"
NET$Status[which(NET$Species=="CRCH" & NET$RepeatID %in% CRCHnn)]="Introduced"
NET$Status[which(NET$Species=="BLBH")]="Native"
NET$Status[which(NET$Species=="BLBH" & NET$RepeatID %in% BLBHnn)]="Introduced"
NET$Status[which(NET$Species=="CCAT" & NET$RepeatID %in% CCATnn)]="Introduced"
NET$Status[which(NET$Species=="LING" & NET$RepeatID %in% LINGnn)]="Introduced"
NET$Status[which(NET$Species=="WSU" & NET$RepeatID %in% WSUnn)]="Introduced"
NET$Status[which(NET$Species=="RMCT" & NET$RepeatID %in% RMCTnn)]="Introduced"
NET$Status[which(NET$Species=="DRUM")]="Native"
NET$Status[which(NET$Species=="DRUM" & NET$RepeatID %in% drumnn)]="Introduced"
NET$Status[which(NET$Species=="PKF")]="Native"
NET$Status[which(NET$Species=="PKF" & NET$RepeatID %in% pkfnn)]="Introduced"
NET$Status[which(NET$Species=="PTMN")]="Native"
NET$Status[which(NET$Species=="PTMN" & NET$RepeatID %in% ptmnnn)]="Introduced"
NET$Status[which(NET$Species=="BRSB")]="Native"
NET$Status[which(NET$Species=="BRSB" & NET$RepeatID %in% brsbnn)]="Introduced"
invasive20s=c("LL","CARP","GSUN","RB","RSSH","EB","SMB","NP")
NET=subset(NET,NET$Status=="Native"|NET$Species%in%invasive20s)
nativeSites=NET%>%
  group_by(Species)%>%
  summarise(sites=length(RepeatID))
#write.csv(nativeSites,"nativesites.csv")

NET2=NET%>%group_by(Species)%>%summarise(net=mean(change), std.E = sd(change)/sqrt(length(change)), n=length(change))
NET2$CI90 = NA
NET2$CI90 = NET2$std.E*1.645
NETsub = subset(NET2, NET2$n>=20)
NETsub30 = subset(NET2, NET2$n>=50)




#---------------COMBINED GUILDS------------------------
traits$ThermalFlow=NA
traits$ThermalFlow=paste(traits$SimpleThermal,traits$LowFlowSens,sep = "")

traits$ThermalLifHis=NA
traits$ThermalLifHis=paste(traits$SimpleThermal,traits$LifeHist,sep = "")

traits$FlowLifHis=NA
traits$FlowLifHis=paste(traits$LowFlowSens,traits$LifeHist,sep = "")

#-----------------------------------------------------


traits2=traits[,c(1,3)]
NET=left_join(NETsub,traits,by="Species")
NETsub=subset(NET,NET$Species!="ONC" & NET$Species!="FMxWSU")
mean(NETsub$net)
write.csv(NETsub, "table1materialWITHconfidence.csv")

NETsub$CommonName[which(NETsub$CommonName=="Northern Redbelly Dace")]="Nor. Redbelly Dace"
NETsub$CommonName[which(NETsub$CommonName=="Rocky Mountain Cutthroat Trout")]="R.M. Cutthroat Trout"


NETsub$Direc=NA
NETsub$Direc[which(NETsub$net>=1)]="Pos"
NETsub$Direc[which(NETsub$net<1)]="Neg"
NETsub$Status=NA
NETsub$Status="Native"
NETsub$Status[which(NETsub$Species%in%invasive20s)]="Introduced"

library(ggbreak)
library(ggh4x)
changenative=NETsub%>%
  ggplot(aes(x=reorder(CommonName,-net),y=net, shape = Status, color=Status))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_segment( aes(x=CommonName, xend=CommonName, y=net-CI90, yend=net+CI90), color="black") +
  geom_point(size=3)+
  theme_light()+
  ylab(label = "Site Occupancy Trend")+
  xlab(label="")+
  scale_color_manual(values = c("#ff9999","#005555"))+
  scale_shape_manual(values=c(15,17))+
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1, by=0.2))+
  coord_flip()+
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,face = 2),
    legend.title = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=10, color = "black"),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))
changenative


ggsave(filename="NativeSpChange_CI_MISSOURIONLY.tiff",dpi = 400, width = 10, height = 8, units = "in")

NETsub$LifeHist[which(NETsub$LifeHist=="Equil")]="Equilibrium"
NETsub$LifeHist[which(NETsub$LifeHist=="Opp")]="Opportunistic"
NETsub$LifeHist[which(NETsub$LifeHist=="Per")]="Periodic"

graphA=NETsub%>%
  ggplot(aes(x = Status, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  xlab("Native vs. Introduced")+
  ggtitle("A. Native v. Introduced")+
  annotate("text", x = 2.3, y = -0.45, label = "p = 0.22")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))


graphB=NETsub%>%filter(SimpleThermal=="Cold"|SimpleThermal=="Cool"|SimpleThermal=="Warm")%>%
  ggplot(aes(x = SimpleThermal, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  xlab("Thermal Guild")+
  ggtitle("B. Thermal Guild")+
  annotate("text", x = 3, y = -0.48, label = "p = 0.18")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))

graphC=NETsub%>%filter(LowFlowSens=="Low"|LowFlowSens=="Moderate"|LowFlowSens=="High")%>%
  ggplot(aes(x = LowFlowSens, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  xlab("Low Flow Sensitivity Guild")+
  ggtitle("C. Flow Preference")+
  annotate("text", x = 2.3, y = -0.45, label = "p = 0.06")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))

graphD=NETsub%>%
  ggplot(aes(x = LifeHist, y= net))+
  geom_hline(yintercept = 0, color=alpha("black", alpha = 0.3), linewidth=1.1)+
  geom_boxplot()+
  theme_light()+
  ylab("Site Occupancy Trend")+
  ggtitle("D. Life History Strategy")+
  annotate("text", x = 3, y = -0.45, label = "p = 0.01")+
  annotate("text", x = 1.1, y = 0.4, label = "a")+
  annotate("text", x = 2.1, y = 0.15, label = "b")+
  annotate("text", x = 3.2, y = 0.3, label = "ab")+
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_rect(color = "black"),
    axis.text.y = element_text(size=12, color = "black"),
    axis.text.x = element_text(size=12, color = "black"),
    axis.title = element_blank(),
    legend.position = "none")+
  theme(ggh4x.axis.ticks.length.minor=rel(1))
graphD


library(ggpubr)


guild.net.graph=ggarrange(graphA, graphB, graphC, graphD,
          nrow = 1)
guild.net.graph

ggsave(filename = "guildNETgraph.tiff", dpi=400, height = 4, width = 15, units = "in")






#---------------ANOVAs & t.tests---------------
#NATIVE v INTRODUCED
NETsubnative=subset(NETsub, NETsub$Status=="Native")
NETsubintroduced=subset(NETsub, NETsub$Status!="Native")
t.test(NETsubnative$net,NETsubintroduced$net)
wilcox.test(NETsubnative$net,NETsubintroduced$net)
median(NETsubnative$net)
median(NETsubintroduced$net)

#THERMAL
NETsub2=NETsub%>%filter(SimpleThermal=="Cold"|SimpleThermal=="Cool"|SimpleThermal=="Warm")
therm.aov=aov(net~SimpleThermal, data=NETsub2)
summary(therm.aov) #p=0.248
#plot(therm.aov)
#kruskal.test(net~SimpleThermal, data=NETsub)
NETsub2%>%group_by(SimpleThermal)%>%summarise(mean(net), median(net))

#Flow
NETsubLow=subset(NETsub, NETsub$LowFlowSens=="Low")
NETsubMod=subset(NETsub, NETsub$LowFlowSens=="Moderate")
t.test(NETsubLow$net,NETsubMod$net)
mean(NETsubLow$net)
mean(NETsubMod$net)

#LifeHist
lifhist.aov=aov(net~LifeHist, data=NETsub)
summary(lifhist.aov) #p=0.014
TukeyHSD(lifhist.aov)
#plot(lifhist.aov)
NETsub%>%group_by(LifeHist)%>%summarise(mean(net), median(net))






###ADD NLCD mode type data in 5km buffer around each point
nlcd.dat=read.csv("nlcdcrop.csv")


#===================================================================================
#============================SECTION B: Missouri Basin ONLY ANALYSES & GRAPHS
#==================================================================================
#Manually removed all green river basin data from nssn$obs
nogreen=read.csv("NOGREENmetrics.csv")
library(lme4)
library(cv)
library(MuMIn)


nlcd.dat=nlcd.dat%>%rename("RepetID"="RepeatID")
nogreen=left_join(nogreen, nlcd.dat, by="RepetID")
colnames(nogreen)
####Variable Correlation
library(PerformanceAnalytics)
mdata=nogreen[,c(16,21,23,24,32,141,153)]
mdata=as.data.frame(mdata)
mdata$geometry=NULL
chart.Correlation(mdata, histogram=TRUE, pch=19)




nogreen%>%
  ggplot(aes(x=crop,y=Turnover))+
  geom_point()


#----LOO CV comparisions function-----
loocv_glm_subsets <- function(response, predictors, data, family = "binomial") {
  # ---- Interpret response & predictors ----
  if (is.character(response)) {
    resp_name <- response[1]
  } else {
    resp_name <- deparse(substitute(response))
  }
  
  if (!resp_name %in% names(data)) {
    stop("Response '", resp_name, "' is not a column in `data`.")
  }
  
  preds <- as.character(predictors)
  missing_cols <- setdiff(c(resp_name, preds), names(data))
  if (length(missing_cols) > 0) {
    stop("These variables are missing from `data`: ",
         paste(missing_cols, collapse = ", "))
  }
  
  p <- length(preds)
  if (p < 1L) stop("Need at least one predictor.")
  
  # ---- Handle family argument ----
  fam <- if (is.character(family)) {
    fam_fun <- get(family, mode = "function", envir = parent.frame())
    fam_fun()
  } else {
    family
  }
  
  # ---- Helper: coerce response to numeric for scoring ----
  coerce_y <- function(y) {
    if (is.matrix(y) && ncol(y) == 2) return(y[, 1] / rowSums(y))
    if (is.factor(y)) return(as.numeric(y) - 1L)
    as.numeric(y)
  }
  
  # ---- Helper: LOOCV for a single model formula ----
  loocv_single <- function(formula, dat, family) {
    mf_all <- model.frame(formula, data = dat)
    y_raw  <- model.response(mf_all)
    y      <- coerce_y(y_raw)
    n_obs  <- nrow(mf_all)
    
    preds_vec <- rep(NA_real_, n_obs)
    
    for (i in seq_len(n_obs)) {
      fit_i <- tryCatch(
        glm(formula,
            data   = mf_all[-i, , drop = FALSE],
            family = family),
        error = function(e) NULL
      )
      
      if (!is.null(fit_i)) {
        preds_vec[i] <- predict(
          fit_i,
          newdata = mf_all[i, , drop = FALSE],
          type = "response"
        )
      }
    }
    
    if (all(is.na(preds_vec))) return(NA_real_)
    mean((y - preds_vec)^2, na.rm = TRUE)   # Brier/MSE
  }
  
  # ---- Build all non-empty subsets ----
  subset_list <- unlist(
    lapply(1:p, function(m) combn(preds, m, simplify = FALSE)),
    recursive = FALSE
  )
  
  results <- list()
  
  ## ---- 1. Intercept-only model ----
  int_form_str <- paste(resp_name, "~ 1")
  int_form <- as.formula(int_form_str)
  
  # full-data fit for AIC
  int_fit <- tryCatch(
    glm(int_form, data = data, family = fam),
    error = function(e) NULL
  )
  int_aic <- if (is.null(int_fit)) NA_real_ else AIC(int_fit)
  
  cv_int <- loocv_single(int_form, data, fam)
  
  results[[1]] <- data.frame(
    formula      = int_form_str,
    predictors   = "(intercept only)",
    n_predictors = 0,
    loocv_score  = cv_int,
    AIC          = int_aic,
    stringsAsFactors = FALSE
  )
  
  ## ---- 2. All models with full interactions ----
  for (j in seq_along(subset_list)) {
    vars <- subset_list[[j]]
    
    rhs      <- paste(vars, collapse = " * ")
    form_str <- paste(resp_name, "~", rhs)
    sub_form <- as.formula(form_str)
    
    # full-data fit for AIC
    fit_full <- tryCatch(
      glm(sub_form, data = data, family = fam),
      error = function(e) NULL
    )
    aic_val <- if (is.null(fit_full)) NA_real_ else AIC(fit_full)
    
    cv_val <- loocv_single(sub_form, data, fam)
    
    results[[j + 1]] <- data.frame(
      formula      = form_str,
      predictors   = paste(vars, collapse = ","),
      n_predictors = length(vars),
      loocv_score  = cv_val,
      AIC          = aic_val,
      stringsAsFactors = FALSE
    )
  }
  
  out <- do.call(rbind, results)
  
  # sort by loocv_score (you could also sort by AIC if you prefer)
  out[order(out$loocv_score), ]
}











###-----TURNOVER-----
ngFULL=glm(Turnover~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(crop)*(Yrange), family="binomial", data=nogreen)
summary(ngFULL)

library(DescTools)
PseudoR2(ngFULL, which = "Nagelkerke")

ngFULL=glm(Turnover~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*scale(crop), family="binomial", data=nogreen,  na.action = na.fail)
summary(ngFULL)
PseudoR2(ngFULL, which = "Nagelkerke")
preds <- c("temp", "barrier", "length", "size", "prosper", "crop")
loo_tab <- loocv_glm_subsets(
  response   = "Turnover",     # <-- pass as string, safest
  predictors = preds,
  data       = nogreen,
  family     = "binomial"
)




#dredge.results=dredge(ngFULL)
#dredge.results=subset(dredge.results,dredge.results$delta<=2)

ngnointeractions=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=nogreen,  na.action = na.fail)
summary(ngnointeractions)
PseudoR2(ngnointeractions, which = "Nagelkerke")

ngTEMP=glm(Turnover~temp, data=nogreen, family = "binomial")
summary(ngTEMP)
PseudoR2(ngTEMP, which = "Nagelkerke")

ngTEMPCROP.turn=glm(Turnover~temp*crop, data=nogreen, family = "binomial")
summary(ngTEMPCROP.turn)
PseudoR2(ngTEMPCROP.turn, which = "Nagelkerke")





ngINTERCEPT=glm(Turnover~1, data = nogreen, family = "binomial")
summary(ngINTERCEPT)

cv(ngTEMP, k="loo")

cv(ngTEMPCROP, k="loo")
#cv(ngINTERCEPT)
cv(ngnointeractions,k="loo")
cv(ngFULL, k="loo")

#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture*"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

imp.turn=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Community Turnover")+
  xlab("") +
  coord_flip() +
  ylim(0,0.6)+
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
imp.turn

nogreen%>%filter(!is.na(pNat))%>%summarise(avgPNAT=mean(pNat))
nogreen%>%filter(!is.na(Turnover))%>%summarise(avgTurnover=mean(Turnover))


####
#----BASIN SPECIFIC MODELS----
nogreen$PU[which(nogreen$PU=="LTMO")]="CHEY"
ngUPMO=subset(nogreen,nogreen$PU=="UPMO")
ngPLAT=subset(nogreen,nogreen$PU=="PLAT")
ngCHEY=subset(nogreen,nogreen$PU=="CHEY")
ngYELL=subset(nogreen,nogreen$PU=="YELL")







#Upper Missouri
ngnointeractionsUPMO=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(ngnointeractionsUPMO)
PseudoR2(ngnointeractionsUPMO, which = "Nagelkerke")
#VAR IMP
coeff=as.data.frame(ngnointeractionsUPMO$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsUPMO$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

turn.upmo=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Upper Missouri Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))
turn.upmo


#Yellowstone
ngnointeractionsYELL=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngYELL,  na.action = na.fail)
summary(ngnointeractionsYELL)
PseudoR2(ngnointeractionsYELL, which = "Nagelkerke")


#VAR IMP
coeff=as.data.frame(ngnointeractionsYELL$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsYELL$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

turn.yell=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Yellowstone Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))
turn.yell



#Cheyenne-Little Missouri
ngnointeractionsCHEY=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngCHEY,  na.action = na.fail)
summary(ngnointeractionsCHEY)
PseudoR2(ngnointeractionsCHEY, which = "Nagelkerke")


#VAR IMP
coeff=as.data.frame(ngnointeractionsCHEY$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsCHEY$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

turn.chey=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Black Hills Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

turn.chey





#Platte
ngnointeractionsPLAT=glm(Turnover~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngPLAT,  na.action = na.fail)
summary(ngnointeractionsPLAT)
PseudoR2(ngnointeractionsPLAT, which = "Nagelkerke")


#VAR IMP
coeff=as.data.frame(ngnointeractionsPLAT$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsPLAT$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

turn.plat=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Platte Community Turnover")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

turn.plat

library(ggpubr)
ggarrange(turn.chey,turn.plat,turn.upmo,turn.yell, labels = "AUTO")
ggsave(filename = "BasinsTURNOVER.jpeg",dpi=400,width = 10, height = 10)

nogreen%>%group_by(PU)%>%
  summarise(avgT=mean(Turnover))















###-----NATIVE PERSISTENCE-----
ngFULL=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(size)*(Yrange), family="binomial", data=nogreen)
summary(ngFULL)

nogreenNATIVE=nogreen%>%filter(!is.na(pNat))
ngFULL=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*scale(pisc), family="binomial", data=nogreenNATIVE,  na.action = na.fail)
summary(ngFULL)
PseudoR2(ngFULL, which = "Nagelkerke")
loo_tab_pNAT <- loocv_glm_subsets(
  response   = "pNat",     # <-- pass as string, safest
  predictors = preds,
  data       = nogreenNATIVE,
  family     = "binomial"
)

ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(pisc)+scale(crop), family="binomial", data=nogreenNATIVE,  na.action = na.fail)
summary(ngnointeractions)

ngTEMP=glm(pNat~temp, data=nogreenNATIVE, family = "binomial")
summary(ngTEMP)
PseudoR2(ngTEMP, which = "Nagelkerke")


ngTEMPCROP.nat=glm(pNat~scale(temp)*scale(crop), data=nogreenNATIVE, family = "binomial")
summary(ngTEMPCROP.nat)
PseudoR2(ngTEMPCROP.nat, which = "Nagelkerke")

ngINTERCEPT=glm(pNat~1, data = nogreen, family = "binomial")
summary(ngINTERCEPT)

cv(ngTEMP, k="loo")
cv(ngTEMPPROSP, k="loo")
cv(ngnointeractions,k="loo")
cv(ngFULL, k="loo")

PseudoR2(ngTEMP, which = "Nagelkerke")
PseudoR2(ngFULL, which = "Nagelkerke")
PseudoR2(ngnointeractions, which = "Nagelkerke")

#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length*"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture*"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

imp.pnat=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylim(0,0.6)+
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12, color = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
imp.pnat


NOGREENsites=nogreen$RepetID





#----Basin Specific Models----
ngUPMO=subset(nogreenNATIVE,nogreenNATIVE$PU=="UPMO")
ngPLAT=subset(nogreenNATIVE,nogreenNATIVE$PU=="PLAT")
ngCHEY=subset(nogreenNATIVE,nogreenNATIVE$PU=="CHEY")
ngYELL=subset(nogreenNATIVE,nogreenNATIVE$PU=="YELL")

nguppmo=glm(pNat~scale(temp)*scale(length)*scale(size)*scale(prosper)*scale(pisc), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(nguppmo)
PseudoR2(nguppmo, which = "Nagelkerke")
#dredge.results.uppmo=dredge(nguppmo)
#dredge.results.uppmo=subset(dredge.results.uppmo,dredge.results.uppmo$delta<2)

#####Upper Missouri basin
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(pisc)+scale(crop), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

nper.upmo=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Up. Missouri Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

nper.upmo






#####Yellowstone basin
ngyell=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*scale(pisc), family="binomial", data=ngYELL,  na.action = na.fail)
summary(ngyell)
PseudoR2(ngyell, which = "Nagelkerke")
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(pisc)+scale(crop), family="binomial", data=ngYELL,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

nper.yell=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Yellowstone Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

nper.yell








#####Little Mo Chey basin
ngchey=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*scale(pisc), family="binomial", data=ngCHEY,  na.action = na.fail)
PseudoR2(ngchey, which = "Nagelkerke")
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(pisc)+scale(crop), family="binomial", data=ngCHEY,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

nper.chey=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Black Hills Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

nper.chey






#####Platte basin
ngplat=glm(pNat~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*scale(pisc), family="binomial", data=ngPLAT,  na.action = na.fail)
PseudoR2(ngplat, which = "Nagelkerke")
ngnointeractions=glm(pNat~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(pisc)+scale(crop), family="binomial", data=ngPLAT,  na.action = na.fail)
summary(ngnointeractions)


#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

nper.plat=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Platte Native Species Persistence")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

nper.plat


ggarrange(nper.chey,nper.plat,nper.upmo,nper.yell, labels = "AUTO")
ggsave(filename = "BasinsNatPersistence.jpeg",dpi=400,width = 10, height = 10)













#Create percent colonizing metrics
nogreen$numSpL = NA
nogreen$numSpL=nogreen$numSpE+nogreen$rchChng
nogreen$pCol = NA
nogreen$pCol = nogreen$nClnzSp/nogreen$numSpL

nogreenCOL=subset(nogreen,nogreen$numSpL>0)

###-----COLONIZATION-----
ngFULL=glm(pCol~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*(Yrange), family="binomial", data=nogreenCOL)
summary(ngFULL)

ngFULL=glm(pCol~scale(temp)*scale(barrier)*scale(length)*scale(size)*scale(prosper)*scale(crop), family="binomial", data=nogreenCOL,  na.action = na.fail)
summary(ngFULL)
PseudoR2(ngFULL, which = "Nagelkerke")
loo_tab_pCOL <- loocv_glm_subsets(
  response   = "pCol",     # <-- pass as string, safest
  predictors = preds,
  data       = nogreenCOL,
  family     = "binomial"
)


ngnointeractions=glm(pCol~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=nogreenCOL,  na.action = na.fail)
summary(ngnointeractions)

library(DescTools)
ngTEMP=glm(pCol~temp, data=nogreenCOL, family = "binomial")
summary(ngTEMP)
PseudoR2(ngTEMP, which="Nagelkerke")

ngTEMPCROP.col=glm(pCol~temp*crop, data=nogreenCOL, family = "binomial")
summary(ngTEMPCROP.col)
PseudoR2(ngTEMPCROP.col, which="Nagelkerke")


ngINTERCEPT=glm(pCol~1, data = nogreenCOL, family = "binomial")
summary(ngINTERCEPT)

cv(ngTEMP, k="loo")
cv(ngTEMPCROP, k="loo")
cv(ngnointeractions,k="loo")
cv(ngFULL, k="loo")
#cv(ngINTERCEPT, k="loo")

#VAR IMP
coeff=as.data.frame(ngnointeractions$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractions$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture*"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

imp.col=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Proportion Colonizing Species")+
  xlab("") +
  coord_flip() +
  ylim(0,0.6)+
  ylab("Variable Importance")+
  theme(axis.text.y=element_text(size=12, color = "black"),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=12, color = "black"))

imp.col



lm.pcol=nogreenCOL%>%
ggplot(aes(x=temp,y=pCol))+
  geom_point(color="lightgrey")+
  geom_smooth(method = "glm",se=F, color="black",method.args = list(family = "binomial"), size=2)+
  ylab(label="Colonization Proportion")+
  xlab(label="Temperature (Celsius)")+
  ylim(0,1)+
  annotate("text", x=13, y=0.9, label= "R^2=0.13", size=5) +
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"))
lm.pcol

lm.pnat=nogreenCOL%>%
  ggplot(aes(x=temp,y=pNat))+
  geom_point(color="lightgrey")+
  geom_smooth(method = "glm",se=F, color="black",method.args = list(family = "binomial"), size=2)+
  ylab(label="Native Species Persistence")+
  ylim(0,1)+
  annotate("text", x=13, y=0.05, label= "R^2=0.07", size=5) +
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
lm.pnat

lm.turn=nogreenCOL%>%
  ggplot(aes(x=temp,y=Turnover))+
  geom_point(color="lightgrey")+
  geom_smooth(method = "glm",se=F, color="black",method.args = list(family = "binomial"), size=2)+
  ylim(0,1)+
  ylab(label = "Community Turnover")+
  annotate("text", x=13, y=0.9, label= "R^2=0.15",size=5) +
  theme_classic()+
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank())
lm.turn



library(ggpubr)

ggarrange(imp.turn,lm.turn,
          imp.pnat,lm.pnat,
          imp.col,lm.pcol,
          ncol = 2, nrow = 3,labels = c("A","D","B","E","C","F"))

ggsave(filename = "importanceandtempgraph.tiff", dpi = 400, height = 10, width = 10, units = "in")


####
#----BASIN SPECIFIC MODELS----
nogreenCOL$PU[which(nogreenCOL$PU=="LTMO")]="CHEY"
ngUPMO=subset(nogreenCOL,nogreenCOL$PU=="UPMO")
ngPLAT=subset(nogreenCOL,nogreenCOL$PU=="PLAT")
ngCHEY=subset(nogreenCOL,nogreenCOL$PU=="CHEY")
ngYELL=subset(nogreenCOL,nogreenCOL$PU=="YELL")







#Upper Missouri
ngnointeractionsUPMO=glm(pCol~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngUPMO,  na.action = na.fail)
summary(ngnointeractionsUPMO)
cv(ngnointeractionsUPMO, k="loo")
PseudoR2(ngnointeractionsUPMO, which = "Nagelkerke")

#VAR IMP
coeff=as.data.frame(ngnointeractionsUPMO$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsUPMO$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

col.upmo=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Up. Missouri Prop. Colonized")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))
col.upmo


#Yellowstone
ngnointeractionsYELL=glm(pCol~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngYELL,  na.action = na.fail)
summary(ngnointeractionsYELL)
cv(ngnointeractionsYELL, k="loo")
PseudoR2(ngnointeractionsYELL, which = "Nagelkerke")

#VAR IMP
coeff=as.data.frame(ngnointeractionsYELL$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsYELL$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

col.yell=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Yellowstone Prop. Colonized")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

col.yell







#Cheyenne-Little Missouri
ngnointeractionsCHEY=glm(pCol~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngCHEY,  na.action = na.fail)
summary(ngnointeractionsCHEY)
cv(ngnointeractionsCHEY, k="loo")
PseudoR2(ngnointeractionsCHEY, which = "Nagelkerke")


#VAR IMP
coeff=as.data.frame(ngnointeractionsCHEY$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsCHEY$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature*"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

col.chey=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Black Hills Prop. Colonized")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))

col.chey





#Platte
ngnointeractionsPLAT=glm(pCol~scale(temp)+scale(barrier)+scale(length)+scale(size)+scale(prosper)+scale(crop), family="binomial", data=ngPLAT,  na.action = na.fail)
summary(ngnointeractionsPLAT)
cv(ngnointeractionsPLAT, k="loo")
PseudoR2(ngnointeractionsPLAT, which = "Nagelkerke")


#VAR IMP
coeff=as.data.frame(ngnointeractionsPLAT$coefficients)
coeff$Variable=NA
coeff$Variable=rownames(coeff)
coeff=coeff%>%rename("Score"="ngnointeractionsPLAT$coefficients")
coeff=coeff[-1,]
coeff$Score2=NA
coeff$Score2=abs(coeff$Score)
coeff$Var2=NA
coeff$Var2[which(coeff$Variable=="scale(temp)")]="Temperature"
coeff$Var2[which(coeff$Variable=="scale(size)")]="Stream Size"
coeff$Var2[which(coeff$Variable=="scale(barrier)")]="Barriers"
coeff$Var2[which(coeff$Variable=="scale(prosper)")]="Flow Permanence"
coeff$Var2[which(coeff$Variable=="scale(length)")]="Fragment Length"
coeff$Var2[which(coeff$Variable=="scale(pisc)")]="Piscivores"
coeff$Var2[which(coeff$Variable=="scale(crop)")]="Crop/Pasture"
coeff$Sign=NA
coeff$Sign[which(coeff$Score>0)]="Pos"
coeff$Sign[which(coeff$Score<0)]="Neg"

col.plat=coeff%>%arrange(Score2) %>%
  mutate(Var2=factor(Var2, levels=Var2)) %>%
  ggplot(aes(x=Var2, y=Score2, colour = Sign)) +
  geom_segment( aes(xend=Var2, y=0,yend=Score2)) +
  geom_point(size=4) +
  theme_light() +
  scale_color_manual(values = c("#ff9999","#018080"))+
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  ) +
  ggtitle(label="Platte Prop. Colonized")+
  xlab("") +
  coord_flip() +
  ylab("Variable Importance")+
  theme(axis.text.y = element_text(size=12))
col.plat



ggarrange(col.chey,col.plat,col.upmo,col.yell, labels = "AUTO")
ggsave(filename = "BasinsColonization.jpeg",dpi=400,width = 10, height = 10)




nogreenCOL%>%group_by(PU)%>%
  summarise(avgPCOL=mean(pCol), basinavgTemp=mean(temp))
nogreenCOL%>%filter(!is.na(pCol))%>%summarise(avgpCol=mean(pCol))

























#####-----INDIVIDUAL SPECIES PERSISTENCE MODELS-----
library(dplyr)
library(ggplot2)

species_cols <- NETsub30$Species
predictors_all <- c("temp", "barrier", "length", "size","prosper", "pisc","crop")

# species that should NOT include "pisc"
no_pisc_species <- c("LL", "GSUN", "SMB", "NP")

models <- setNames(vector("list", length(species_cols)), species_cols)

for (i in species_cols) {
  cat("\n--------------------------------------------\n")
  cat("Fitting:", i, "\n")
  
  # Choose predictors
  predictors <- if (i %in% no_pisc_species) {
    setdiff(predictors_all, "pisc")
  } else {
    predictors_all
  }
  
  # Keep rows with no NAs in the response or predictors
  keep <- complete.cases(nogreen[, c(i, predictors)])
  i.sub <- nogreen[keep, , drop = FALSE]
  
  # Build formula and fit model
  terms_scaled <- paste0("scale(", predictors, ")")
  fml <- reformulate(terms_scaled, response = i)
  i.model <- glm(fml, family = binomial, data = i.sub)
  models[[i]] <- i.model
  
  # Print model summary
  print(summary(i.model))
  
  # --- Variable importance plot ---
  coeff <- as.data.frame(i.model$coefficients)
  coeff$Variable <- rownames(coeff)
  names(coeff)[1] <- "Score"
  coeff <- coeff[-1, , drop = FALSE]  # remove intercept
  coeff$Score2 <- abs(coeff$Score)
  coeff$Sign <- ifelse(coeff$Score > 0, "Pos", "Neg")
  
  coeff$Var2 <- NA
  coeff$Var2[coeff$Variable == "scale(temp)"]    <- "Temperature"
  coeff$Var2[coeff$Variable == "scale(barrier)"] <- "Barriers"
  coeff$Var2[coeff$Variable == "scale(length)"]  <- "Fragment Length"
  coeff$Var2[coeff$Variable == "scale(prosper)"] <- "Flow Permanence"
  coeff$Var2[coeff$Variable == "scale(pisc)"]    <- "Piscivores"
  coeff$Var2[coeff$Variable == "scale(size)"] <- "Stream Size"
  coeff$Var2[coeff$Variable == "scale(crop)"] <- "Crop/Pasture"
  coeff <- coeff %>%
    arrange(Score2) %>%
    mutate(Var2 = factor(Var2, levels = Var2))
  
  p <- ggplot(coeff, aes(x = Var2, y = Score2, colour = Sign)) +
    geom_segment(aes(xend = Var2, y = 0, yend = Score2)) +
    geom_point(size = 4) +
    theme_light() +
    scale_color_manual(values = c("#ff9999", "#018080")) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(label = paste(i, "- Variable Importance")) +
    xlab("") +
    coord_flip() +
    ylab("Absolute Coefficient")
  
  print(p)
  
  # Optional: save plot
  # ggsave(filename = paste0("varimp_", i, ".png"), plot = p, width = 6, height = 4)
}





# --- Build a species x predictors table (signed version) ---

predictors_all <- c("temp", "barrier", "length", "prosper", "size","pisc","crop")

term_name <- function(var) paste0("scale(", var, ")")

# choose scaling method: "sum" or "max"
scale_method <- "sum"  # or "max"

make_species_row <- function(sp, fit) {
  sm <- summary(fit)
  coefs <- sm$coefficients
  
  # Pull coefficients and p-values for our predictors
  est <- sapply(predictors_all, function(v) {
    rn <- term_name(v)
    if (rn %in% rownames(coefs)) coefs[rn, "Estimate"] else NA_real_
  })
  pvals <- sapply(predictors_all, function(v) {
    rn <- term_name(v)
    if (rn %in% rownames(coefs)) coefs[rn, "Pr(>|z|)"] else NA_real_
  })
  
  # Scale coefficients by sum or max of absolute values
  denom <- if (scale_method == "sum") sum(abs(est), na.rm = TRUE) else max(abs(est), na.rm = TRUE)
  if (is.finite(denom) && denom != 0) {
    scaled_est <- est / denom
  } else {
    scaled_est <- est
  }
  
  # Format with asterisk for p < 0.1 (but keep sign)
  fmt <- function(val, p) {
    if (is.na(val)) return(NA_character_)
    paste0(format(round(val, 3), nsmall = 3),
           ifelse(!is.na(p) && p < 0.1, "*", ""))
  }
  cells <- mapply(fmt, scaled_est, pvals, USE.NAMES = FALSE)
  names(cells) <- predictors_all
  
  # AIC
  aic_val <- AIC(fit)
  
  # Return one row
  data.frame(
    species = sp,
    as.list(cells),
    AIC = round(aic_val, 2),
    check.names = FALSE
  )
}

# Build the full table
species_vec <- names(models)
rows <- lapply(species_vec, function(sp) make_species_row(sp, models[[sp]]))
varimp_table_signed <- do.call(rbind, rows)
varimp_table_signed <- varimp_table_signed[, c("species", predictors_all, "AIC")]

# Preview
print(varimp_table_signed)

write.csv(varimp_table_signed, "species_var_importance_signed.csv", row.names = FALSE)


















#####INDIVIDUAL SPECIES COLONIZATION-PERSISTENCE MODELS
NET=read.csv("PersistenceMetrics.csv")
nogreen=read.csv("NOGREENmetrics.csv")
NOGREENsites=nogreen$RepetID


removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
NET=left_join(NET,removes,by="RepeatID")
NET=subset(NET,NET$GearMiss!="Y")
NET=subset(NET, NET$F_MAUG_<1000)
NET=subset(NET, NET$RepeatID!=234 & NET$RepeatID!=74& NET$RepeatID!=237& NET$RepeatID!=238)

NET=subset(NET,NET$RepeatID%in%NOGREENsites)

NET=NET%>%pivot_longer(cols = c(28:121),names_to = "Species",values_to = "change")

####FROM HERE, MUST RUN "DatasetPrep.R" for section
NET=subset(NET,!is.na(NET$change))
NET$Species[which(NET$Species=="RMCOT" | NET$Species=="COLCOT")]="MOTCOT"
traits=read.csv("traits.csv")
traits=traits%>%rename("Species"="Code")
NET=left_join(NET,traits,by="Species")
NET$Status[which(NET$Species=="RDSH")]="Native"
NET$Status[which(NET$Species=="RDSH" & NET$RepeatID %in% rdshnn)]="Introduced"
NET$Status[which(NET$Species=="LKCH" & NET$RepeatID %in% lkchnn)]="Introduced"
NET$Status[which(NET$Species=="FHMN" & NET$RepeatID %in% FHMNnn)]="Introduced"
NET$Status[which(NET$Species=="LNDC" & NET$RepeatID %in% LNDCnn)]="Introduced"
NET$Status[which(NET$Species=="CRCH" & NET$RepeatID %in% CRCHnn)]="Introduced"
NET$Status[which(NET$Species=="BLBH")]="Native"
NET$Status[which(NET$Species=="BLBH" & NET$RepeatID %in% BLBHnn)]="Introduced"
NET$Status[which(NET$Species=="CCAT" & NET$RepeatID %in% CCATnn)]="Introduced"
NET$Status[which(NET$Species=="LING" & NET$RepeatID %in% LINGnn)]="Introduced"
NET$Status[which(NET$Species=="WSU" & NET$RepeatID %in% WSUnn)]="Introduced"
NET$Status[which(NET$Species=="RMCT" & NET$RepeatID %in% RMCTnn)]="Introduced"
NET$Status[which(NET$Species=="DRUM")]="Native"
NET$Status[which(NET$Species=="DRUM" & NET$RepeatID %in% drumnn)]="Introduced"
NET$Status[which(NET$Species=="PKF")]="Native"
NET$Status[which(NET$Species=="PKF" & NET$RepeatID %in% pkfnn)]="Introduced"
NET$Status[which(NET$Species=="PTMN")]="Native"
NET$Status[which(NET$Species=="PTMN" & NET$RepeatID %in% ptmnnn)]="Introduced"
NET$Status[which(NET$Species=="BRSB")]="Native"
NET$Status[which(NET$Species=="BRSB" & NET$RepeatID %in% brsbnn)]="Introduced"
invasive20s=c("LL","CARP","GSUN","RB","RSSH","EB","SMB","NP")
NET=subset(NET,NET$Status=="Native"|NET$Species%in%invasive20s)
nativeSites=NET%>%
  group_by(Species)%>%
  summarise(sites=length(RepeatID))
#write.csv(nativeSites,"nativesites.csv")

NET2=NET%>%group_by(Species)%>%summarise(net=mean(change), std.E = sd(change)/sqrt(length(change)), n=length(change))
NET2$CI90 = NA
NET2$CI90 = NET2$std.E*1.645
NETsub = subset(NET2, NET2$n>=20)
NETsub30 = subset(NET2, NET2$n>=30)





NETcol=NET
NETcol$change[which(NETcol$change==0)]=1
NETcol$change[which(NETcol$change==-1)]=0

prosperdat=nogreen[,c(19,141)]
prosperdat=prosperdat%>%rename("RepeatID"="RepetID")
NETcol=left_join(NETcol, prosperdat, by="RepeatID")

NETcol2=NETcol%>%group_by(Species)%>%summarise(n=length(CommonName))
NETcol30=subset(NETcol2,NETcol2$n>=30)
NETcol30sp=NETcol30$Species
nlcd.dat=read.csv("nlcdcrop.csv")
NETcol=subset(NETcol,NETcol$Species%in%NETcol30sp)
NETcol=left_join(NETcol, nlcd.dat, by="RepeatID")
NETcol=NETcol%>%rename("temp"="S1_93_11","barrier"="DSbarrier",
                       "length"="Length_km", "size"="F_MAUG_HIS",
                       "pisc"="nnPreds")

NETcol3=NETcol[,c(2,3,15,18,19,22,24,25,26,33,36,37,55,56)]
NETcol3=NETcol3%>%pivot_wider(names_from = "Species", values_from = "change")
#write.csv(NETcol3,"NETcol3.csv")

#####-----INDIVIDUAL SPECIES COLONIZATION-PERSISTENCE MODELS-----
library(dplyr)
library(ggplot2)

species_cols <- NETcol30sp
predictors_all <- c("temp", "barrier", "length", "size","prosper", "pisc","crop")

# species that should NOT include "pisc"
no_pisc_species <- c("LL", "GSUN", "SMB", "NP")

models <- setNames(vector("list", length(species_cols)), species_cols)

for (i in species_cols) {
  cat("\n--------------------------------------------\n")
  cat("Fitting:", i, "\n")
  
  # Choose predictors
  predictors <- if (i %in% no_pisc_species) {
    setdiff(predictors_all, "pisc")
  } else {
    predictors_all
  }
  
  # Keep rows with no NAs in the response or predictors
  keep <- complete.cases(NETcol3[, c(i, predictors)])
  i.sub <- NETcol3[keep, , drop = FALSE]
  
  # Build formula and fit model
  terms_scaled <- paste0("scale(", predictors, ")")
  fml <- reformulate(terms_scaled, response = i)
  i.model <- glm(fml, family = binomial, data = i.sub)
  models[[i]] <- i.model
  
  # Print model summary
  print(summary(i.model))
  
  # --- Variable importance plot ---
  coeff <- as.data.frame(i.model$coefficients)
  coeff$Variable <- rownames(coeff)
  names(coeff)[1] <- "Score"
  coeff <- coeff[-1, , drop = FALSE]  # remove intercept
  coeff$Score2 <- abs(coeff$Score)
  coeff$Sign <- ifelse(coeff$Score > 0, "Pos", "Neg")
  
  coeff$Var2 <- NA
  coeff$Var2[coeff$Variable == "scale(temp)"]    <- "Temperature"
  coeff$Var2[coeff$Variable == "scale(barrier)"] <- "Barriers"
  coeff$Var2[coeff$Variable == "scale(length)"]  <- "Fragment Length"
  coeff$Var2[coeff$Variable == "scale(prosper)"] <- "Flow Permanence"
  coeff$Var2[coeff$Variable == "scale(pisc)"]    <- "Piscivores"
  coeff$Var2[coeff$Variable == "scale(size)"] <- "Stream Size"
  
  coeff <- coeff %>%
    arrange(Score2) %>%
    mutate(Var2 = factor(Var2, levels = Var2))
  
  p <- ggplot(coeff, aes(x = Var2, y = Score2, colour = Sign)) +
    geom_segment(aes(xend = Var2, y = 0, yend = Score2)) +
    geom_point(size = 4) +
    theme_light() +
    scale_color_manual(values = c("#ff9999", "#018080")) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(label = paste(i, "- Variable Importance")) +
    xlab("") +
    coord_flip() +
    ylab("Absolute Coefficient")
  
  print(p)
  
  # Optional: save plot
  # ggsave(filename = paste0("varimp_", i, ".png"), plot = p, width = 6, height = 4)
}



library(dplyr)
library(ggplot2)



# species column names (make sure this is a character vector of names)
species_cols <- NETcol30sp          # e.g. names(NETcol3)[13:32]

predictors_all <- c("temp", "barrier", "length", "size", "prosper", "pisc","crop")

# species that should NOT include "pisc"
no_pisc_species <- c("LL", "GSUN", "SMB", "NP")

models <- setNames(vector("list", length(species_cols)), species_cols)

for (i in species_cols) {
  cat("\n--------------------------------------------\n")
  cat("Fitting:", i, "\n")
  
  # Choose predictors
  predictors <- if (i %in% no_pisc_species) {
    setdiff(predictors_all, "pisc")
  } else {
    predictors_all
  }
  
  # Keep rows with no NAs in the response or predictors
  keep  <- complete.cases(NETcol3[, c(i, predictors)])
  i.sub <- NETcol3[keep, , drop = FALSE]
  
  # Build formula and fit model
  terms_scaled <- paste0("scale(", predictors, ")")
  fml <- reformulate(terms_scaled, response = i)
  i.model <- glm(fml, family = binomial, data = i.sub)
  models[[i]] <- i.model
  
  # Print model summary
  print(summary(i.model))
  
  # --- Variable importance plot ---
  coeff <- as.data.frame(i.model$coefficients)
  coeff$Variable <- rownames(coeff)
  names(coeff)[1] <- "Score"
  coeff <- coeff[-1, , drop = FALSE]  # remove intercept
  coeff$Score2 <- abs(coeff$Score)
  coeff$Sign <- ifelse(coeff$Score > 0, "Pos", "Neg")
  
  coeff$Var2 <- NA
  coeff$Var2[coeff$Variable == "scale(temp)"]    <- "Temperature"
  coeff$Var2[coeff$Variable == "scale(barrier)"] <- "Barriers"
  coeff$Var2[coeff$Variable == "scale(length)"]  <- "Fragment Length"
  coeff$Var2[coeff$Variable == "scale(prosper)"] <- "Flow Permanence"
  coeff$Var2[coeff$Variable == "scale(pisc)"]    <- "Piscivores"
  coeff$Var2[coeff$Variable == "scale(size)"]    <- "Stream Size"
  
  coeff <- coeff %>%
    arrange(Score2) %>%
    mutate(Var2 = factor(Var2, levels = Var2))
  
  p <- ggplot(coeff, aes(x = Var2, y = Score2, colour = Sign)) +
    geom_segment(aes(xend = Var2, y = 0, yend = Score2)) +
    geom_point(size = 4) +
    theme_light() +
    scale_color_manual(values = c("#ff9999", "#018080")) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(label = paste(i, "- Variable Importance")) +
    xlab("") +
    coord_flip() +
    ylab("Absolute Coefficient")
  
  print(p)
}







# --- Build a species x predictors table (signed version) ---

predictors_all <- c("temp", "barrier", "length", "prosper", "size","pisc","crop")

term_name <- function(var) paste0("scale(", var, ")")

# choose scaling method: "sum" or "max"
scale_method <- "sum"  # or "max"

make_species_row <- function(sp, fit) {
  sm <- summary(fit)
  coefs <- sm$coefficients
  
  # Pull coefficients and p-values for our predictors
  est <- sapply(predictors_all, function(v) {
    rn <- term_name(v)
    if (rn %in% rownames(coefs)) coefs[rn, "Estimate"] else NA_real_
  })
  pvals <- sapply(predictors_all, function(v) {
    rn <- term_name(v)
    if (rn %in% rownames(coefs)) coefs[rn, "Pr(>|z|)"] else NA_real_
  })
  
  # Scale coefficients by sum or max of absolute values
  denom <- if (scale_method == "sum") sum(abs(est), na.rm = TRUE) else max(abs(est), na.rm = TRUE)
  if (is.finite(denom) && denom != 0) {
    scaled_est <- est / denom
  } else {
    scaled_est <- est
  }
  
  # Format with asterisk for p < 0.1 (but keep sign)
  fmt <- function(val, p) {
    if (is.na(val)) return(NA_character_)
    paste0(format(round(val, 3), nsmall = 3),
           ifelse(!is.na(p) && p < 0.1, "*", ""))
  }
  cells <- mapply(fmt, scaled_est, pvals, USE.NAMES = FALSE)
  names(cells) <- predictors_all
  
  # AIC
  aic_val <- AIC(fit)
  
  # Return one row
  data.frame(
    species = sp,
    as.list(cells),
    AIC = round(aic_val, 2),
    check.names = FALSE
  )
}

# Build the full table
species_vec <- names(models)
rows <- lapply(species_vec, function(sp) make_species_row(sp, models[[sp]]))
varimp_table_signed <- do.call(rbind, rows)
varimp_table_signed <- varimp_table_signed[, c("species", predictors_all, "AIC")]

# Preview
print(varimp_table_signed)

write.csv(varimp_table_signed, "species_var_importance_COLandPERS.csv", row.names = FALSE)








#####################################################################################################
#============================SECTION C: ORDINATION
#####################################################################################################
#============================Dataset preparation

##Load repeated site data
wc=read.csv("S2DR_v_1_1_WILDCARD.csv")

wc=wc[-240,]#remove sample that was added twice
#wcincluded=subset(wc, wc$RepeatID%in%included)
wc_nn=wc
wc_nn=wc_nn[,c(29,137)]
wc_nn=wc_nn%>%group_by(RepeatID)%>%summarise(pu=first(PU))
wc_nn=wc_nn%>%rename("PU"="pu")
write.csv(wc_nn,file = "processingunits.csv")

############Correct for samples that differentially classified Hybognathus sp.
wchybogchecks=subset(wc, wc$HYBOG==1 | wc$WS_PLMN==1)
wchybogchecks2=as.list(wchybogchecks$RepeatID)
wchybogchecks=subset(wc,wc$RepeatID%in%wchybogchecks2)
wchybogchecks=wchybogchecks[,c(29,30,74:77,134,136)]
wc$PLMN[which(wc$RepeatID==23 & wc$TIME=="LATE")]=0
wc$WS_PLMN[which(wc$RepeatID==23 & wc$TIME=="LATE")]=1
wc$PLMN[which(wc$RepeatID==167 & wc$TIME=="LATE")]=0
wc$WS_PLMN[which(wc$RepeatID==167 & wc$TIME=="LATE")]=1
wc$PLMN[which(wc$RepeatID==266 & wc$TIME=="LATE")]=0
wc$WS_PLMN[which(wc$RepeatID==266 & wc$TIME=="LATE")]=1
wc$PLMN[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=0
wc$WSMN[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=0
wc$HYBOG[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=1
wc$numSp[which(wc$RepeatID==341 & wc$TIME=="EARLY")]=6
wc$PLMN[which(wc$RepeatID==387 & wc$TIME=="EARLY")]=0
wc$HYBOG[which(wc$RepeatID==387 & wc$TIME=="EARLY")]=1
wc=wc%>%rename("temp"="S1_93_11", "size"="F_MAUG_HIS", "barrier"="DSbarrier","length"="Length_km")


removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
wc=left_join(wc,removes,by="RepeatID")
wc=subset(wc,wc$GearMiss!="Y")
wc=subset(wc, wc$size<1000)
wc=subset(wc, wc$RepeatID!=234 & wc$RepeatID!=74& wc$RepeatID!=237& wc$RepeatID!=238)
#write.csv(NET, "PersistenceMetrics_filtered.csv")
wc=subset(wc,wc$PU!="GREEN")
wc=subset(wc,wc$RepeatID!=368)
wc=subset(wc,wc$RepeatID!=431)
wc=subset(wc,wc$RepeatID!=358)
wc=subset(wc,wc$RepeatID!=356)
wc=subset(wc,wc$RepeatID!=397)
wc=subset(wc,wc$RepeatID!=329)

wc$Site_ID=NA
wc$Site_ID=paste(wc$RepeatID, wc$TIME, sep = "_")
wc.comm=wc[,c(139,42:135)]
colSums(wc.comm[2:95])
wc.comm <- wc.comm[, c(1, which(colSums(wc.comm[, 2:95]) > 0) + 1)]

wc.meta=wc[,c(29,30,137,26,36,39,41)]


#####################################################
#---------------------NMDS--------------------------
####################################################
library(vegan)



rowSums(wc.comm[,-1])
# Make first column row names and convert to matrix
wc.mat <- as.matrix(wc.comm[,-1])
rownames(wc.mat) <- wc.comm[,1]
# Example:
nmds <- metaMDS(wc.mat, distance = "jaccard", binary=T, trymax = 20,autotransform = F)

sc <- as.data.frame(scores(nmds, display = "sites"))
sc=cbind(sc, wc.meta)
sc$NMDS1=as.numeric(sc$NMDS1)
sc$NMDS2=as.numeric(sc$NMDS2)
# Ensure TIME and PU are proper factors
sc$TIME <- factor(sc$TIME, levels = c("EARLY", "LATE"))

sc$temp_cat <- NA
sc$temp_cat[sc$temp >= 20] <- "WARM"
sc$temp_cat[sc$temp <  20] <- "COOL"
sc$temp_cat <- factor(sc$temp_cat, levels = c("COOL","WARM"))  # nice legend order

# Colors for temperature categories
temp_cols <- setNames(c("steelblue", "tomato"), levels(sc$temp_cat))

# Base empty plot (no points)
plot(sc$NMDS1, sc$NMDS2, type = "n", xlab = "NMDS1", ylab = "NMDS2")

# Arrows per site, colored by temp category at that site; ensure EARLY -> LATE
by(sc, sc$RepeatID, function(df) {
  if (nrow(df) == 2) {
    df <- df[order(df$TIME), ]  # EARLY (row 1) -> LATE (row 2)
    col_here <- temp_cols[as.character(df$temp_cat[1])]  # same site temp_cat
    arrows(df$NMDS1[1], df$NMDS2[1],
           df$NMDS1[2], df$NMDS2[2],
           length = 0.06, angle = 20, lwd = 2,
           col = col_here)
  }
})

# Legend for temperature categories
legend("bottomleft", lwd = 2, col = temp_cols,
       legend = names(temp_cols), title = "Temperature", bty = "n")
title("Community Trajectories by Stream Temperature")

pl=recordPlot()
pl


#RECREATE IN GGPLOT
library(dplyr)
library(ggplot2)
library(grid)   # for unit()


# Make a data frame of EARLY -> LATE arrows for each RepeatID
arrows_df <- sc %>%
  group_by(RepeatID) %>%
  filter(n() == 2) %>%        # keep only pairs
  arrange(TIME) %>%           # EARLY then LATE (given factor levels)
  summarise(
    NMDS1_start = first(NMDS1),
    NMDS2_start = first(NMDS2),
    NMDS1_end   = last(NMDS1),
    NMDS2_end   = last(NMDS2),
    temp_cat    = first(temp_cat),
    .groups = "drop"
  )

# Recreate your color mapping (or reuse existing temp_cols)
temp_cols <- setNames(c("steelblue", "tomato"),
                      c("COOL", "WARM"))

# ggplot version of the NMDS arrows
plarrow=ggplot(arrows_df) +
  geom_segment(
    aes(x = NMDS1_start, y = NMDS2_start,
        xend = NMDS1_end, yend = NMDS2_end,
        colour = temp_cat),
    arrow     = arrow(length = unit(0.2, "cm"), angle = 20),
    linewidth = 0.8
  ) +
  scale_colour_manual(values = temp_cols, name = "Temperature") +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "Community Trajectories by Temperature"
  ) +
  coord_equal() +
  theme_bw()+
  theme(
    legend.position = c(0.15, 0.15),   # (x, y) in relative plot coordinates
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )



plarrow






library(dplyr)

# Make sure TIME is ordered
sc$TIME <- factor(sc$TIME, levels = c("EARLY", "LATE"))
sc$PU[which(sc$PU=="CHEY")]="Black Hills"
sc$PU[which(sc$PU=="LTMO")]="Black Hills"
sc$PU[which(sc$PU=="PLAT")]="Platte"
sc$PU[which(sc$PU=="UPMO")]="U. Missouri"
sc$PU[which(sc$PU=="YELL")]="Yellowstone"

sc$PU=as.factor(sc$PU)

# Calculate arrow lengths (distance between EARLY and LATE)
arrow_lengths <- sc %>%
  arrange(RepeatID, TIME) %>%
  group_by(RepeatID) %>%
  filter(n() == 2) %>%                    # only include sites with both samples
  summarize(
    PU        = first(PU),
    temp_cat  = first(temp_cat),
    length    = sqrt((NMDS1[2] - NMDS1[1])^2 + (NMDS2[2] - NMDS2[1])^2)
  )

# Inspect results
head(arrow_lengths)

# Average arrow length by temperature category
arrows=arrow_lengths %>%
  group_by(temp_cat) %>%
  summarize(
    mean_length = mean(length, na.rm = TRUE),
    sd_length   = sd(length, na.rm = TRUE),
    se_length = sd(length)/sqrt(length(length)),
    n           = n()
  )

arrows

meanarrow=
arrows%>%
ggplot(aes(x=temp_cat, y= mean_length, color=temp_cat))+
  scale_color_manual(values = c("steelblue", "tomato"))+
  ylim(limits=c(0.2,0.8))+
  ylab("NMDS Vector Length")+
  xlab("Stream Temperature")+
  ggtitle("Magnitude of Change by Temperature")+
  geom_pointrange(ymin=arrows$mean_length-arrows$se_length, 
                ymax=arrows$mean_length+arrows$se_length,
                linewidth=1.5, size=1)+
  annotate("text",label="a", x=1.05, y=0.5, size=5)+
  annotate("text",label="b", x=2.05, y=0.72, size=5)+
  theme_classic()+
  theme(axis.title = element_text(size=12, face = "bold"),
        axis.text = element_text(size=12, color = "black"),
        plot.title = element_text(size=14, color = "black"),
        legend.position = "none")
meanarrow



sc$PU   <- factor(sc$PU)



# 2. Build your community matrix and distance matrix
# (Assuming wc.mat or comm is your sites × species presence/absence matrix)
dist_jac <- vegdist(wc.mat, method = "jaccard", binary = TRUE)

# 3. Combine metadata (must match row order of wc.mat)
env <- sc[, c("RepeatID", "PU", "TIME", "temp_cat")]

# 4. PERMDISP setup
# Test whether within-basin dispersions differ between TIME periods
# (First, group by PU × TIME)
group <- interaction(env$PU, env$TIME)

disp <- betadisper(dist_jac, group = group)

# 5. Test for differences in dispersion between time periods within each basin
anova(disp)              # overall test, p-value <0.01
perm_res=permutest(disp, pairwise = TRUE)  # permutation test (pairwise comparisons)

perm_res$pairwise$permuted


# 6. Extract mean distances to centroid for plotting
centroid_df <- as.data.frame(disp$group)
centroid_df$distance <- disp$distances
centroid_df <- cbind(centroid_df, env)

aggregate(distance ~ PU + TIME, data = centroid_df, FUN = mean)





#BIOTIC HOMOGENIZATION TEST WITH TEMPERATURE
sc$PU   <- factor(sc$PU)
sc$temp_cat = factor(sc$temp_cat)
# 2. Build your community matrix and distance matrix
# (Assuming wc.mat or comm is your sites × species presence/absence matrix)
dist_jac <- vegdist(wc.mat, method = "jaccard", binary = TRUE)

# 3. Combine metadata (must match row order of wc.mat)
env <- sc[, c("RepeatID", "PU", "TIME", "temp_cat")]

# 4. PERMDISP setup
# Test whether within-basin dispersions differ between TIME periods
# (First, group by PU × TIME)
group <- interaction(env$PU, env$temp_cat,env$TIME)

disp <- betadisper(dist_jac, group = group)

# 5. Test for differences in dispersion between time periods within each basin
anova(disp)              # overall test, p-value <0.01
perm_res=permutest(disp, pairwise = TRUE)  # permutation test (pairwise comparisons)

perm_res$pairwise$permuted


# 6. Extract mean distances to centroid for plotting
centroid_df <- as.data.frame(disp$group)
centroid_df$distance <- disp$distances
centroid_df <- cbind(centroid_df, env)

homog.test=aggregate(distance ~ PU + temp_cat+TIME, data = centroid_df, FUN = mean)
homog.test=homog.test%>%pivot_wider(names_from = "TIME", values_from = "distance")
homog.test$change=NA
homog.test$change=homog.test$LATE-homog.test$EARLY
homog.test%>%
  ggplot(aes(x=temp_cat, y=change, shape = PU, color=temp_cat))+
  scale_shape_manual(values = 15:18)+
  scale_color_manual(values = c("steelblue", "tomato"))+
  geom_point(size=4)+
  theme_classic()





#------Trajectories by Basin---------
# Ensure TIME and PU are proper factors
sc$TIME <- factor(sc$TIME, levels = c("EARLY", "LATE"))
sc$PU <- factor(sc$PU)
unique(sc$PU)
sc$temp_cat <- NA
sc$temp_cat[sc$temp >= 20] <- "WARM"
sc$temp_cat[sc$temp <  20] <- "COOL"
sc$temp_cat <- factor(sc$temp_cat, levels = c("COOL","WARM"))  # nice legend order

# Colors for basins
basin_cols <- setNames(c("#A38DBA","#F88B78","#ADDEFF","#FEE658"), levels(sc$PU))

# Base empty plot (no points)
plot(sc$NMDS1, sc$NMDS2, type = "n", xlab = "NMDS1", ylab = "NMDS2")

# Arrows per site, colored by temp category at that site; ensure EARLY -> LATE
by(sc, sc$RepeatID, function(df) {
  if (nrow(df) == 2) {
    df <- df[order(df$TIME), ]  # EARLY (row 1) -> LATE (row 2)
    col_here <- basin_cols[as.character(df$PU[1])]  # same site temp_cat
    arrows(df$NMDS1[1], df$NMDS2[1],
           df$NMDS1[2], df$NMDS2[2],
           length = 0.06, angle = 20, lwd = 2,
           col = col_here)
  }
})

# Legend for temperature categories
legend("bottomleft", lwd = 2, col = basin_cols,
       legend = names(basin_cols), title = "Basin", bty = "n")
title("Community Trajectories by Basin")



pl2=recordPlot()
pl2




###----GGPLOT NMDS----
library(dplyr)
library(ggplot2)
library(grid)   # arrow units

# --- Recode basins ---

# --- Basin color palette ---
basin_cols <- setNames(
  c("#A38DBA", "#F88B78", "#ADDEFF", "#FEE658"),
  levels(sc$PU)
)

# --- Build arrows: EARLY -> LATE for each RepeatID ---
arrows_basin_df <- sc %>%
  group_by(RepeatID) %>%
  filter(n() == 2) %>%
  arrange(TIME) %>%
  summarise(
    NMDS1_start = first(NMDS1),
    NMDS2_start = first(NMDS2),
    NMDS1_end   = last(NMDS1),
    NMDS2_end   = last(NMDS2),
    PU          = first(PU),
    .groups = "drop"
  )

# --- ggplot with legend inside ---
arrows_basin=ggplot(arrows_basin_df) +
  geom_segment(
    aes(x = NMDS1_start, y = NMDS2_start,
        xend = NMDS1_end, yend = NMDS2_end,
        colour = PU),
    arrow     = arrow(length = unit(0.2, "cm"), angle = 20),
    linewidth = 1
  ) +
  scale_colour_manual(values = basin_cols, name = "Basin") +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    title = "Community Trajectories by Basin"
  ) +
  coord_equal() +
  theme_bw() +
  theme(
    legend.position = c(0.2, 0.2),   # (x, y) in relative plot coordinates
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

arrows_basin















# Make sure TIME is ordered
sc$TIME <- factor(sc$TIME, levels = c("EARLY", "LATE"))

# Calculate arrow lengths (distance between EARLY and LATE)
arrow_lengths <- sc %>%
  arrange(RepeatID, TIME) %>%
  group_by(RepeatID) %>%
  filter(n() == 2) %>%                    # only include sites with both samples
  summarize(
    PU        = first(PU),
    temp_cat  = first(temp_cat),
    length    = sqrt((NMDS1[2] - NMDS1[1])^2 + (NMDS2[2] - NMDS2[1])^2)
  )

# Inspect results
head(arrow_lengths)



######STATISTICAL TESTS: temp
arrows_warm=subset(arrow_lengths,arrow_lengths$temp_cat=="WARM")
arrows_cool=subset(arrow_lengths,arrow_lengths$temp_cat=="COOL")

t.test(arrows_cool$length,arrows_warm$length) #p=0.02
wilcox.test(arrows_cool$length,arrows_warm$length)


######STATISTICAL TESTS:basin
PUaov=aov(arrow_lengths$length~arrow_lengths$PU)
summary(PUaov)
TukeyHSD(PUaov) #Upper MO sig dif than others



# Average arrow length by temperature category
arrows=arrow_lengths %>%
  group_by(PU) %>%
  summarize(
    mean_length = mean(length, na.rm = TRUE),
    sd_length   = sd(length, na.rm = TRUE),
    se_length = sd(length)/sqrt(length(length)),
    n           = n()
  )

arrows
arrows$PU=as.character(arrows$PU)
arrows$PU[which(arrows$PU=="CHEY")]="Black Hills"
arrows$PU[which(arrows$PU=="PLAT")]="Platte"
arrows$PU[which(arrows$PU=="UPMO")]="U. Missouri"
arrows$PU[which(arrows$PU=="YELL")]="Yellowstone"
arrows$PU=as.factor(arrows$PU)
meanarrowPU=
  arrows%>%
  ggplot(aes(x=PU, y= mean_length, color=PU))+
  ylab("NMDS Vector Length")+
  xlab("Major Basin")+
  ylim(0,1)+
  ggtitle("Magnitude of Change by Basin")+
  geom_pointrange(ymin=arrows$mean_length-arrows$se_length, 
                  ymax=arrows$mean_length+arrows$se_length,
                  linewidth=1.5, size=1)+
  scale_colour_manual(values = basin_cols, name = "Basin")+
  
  theme_classic()+
  annotate("text",label="a", x=1.07, y=0.45, size=5)+
  annotate("text",label="a", x=2.07, y=0.4, size=5)+
  annotate("text",label="b", x=3.09, y=0.95, size=5)+
  annotate("text",label="a", x=4.07, y=0.56, size=5)+
  theme(axis.title = element_text(size=12, face = "bold"),
        axis.text = element_text(size=12, color = "black"),
        plot.title = element_text(size=14, color = "black"),
        legend.position = "none")

meanarrowPU
library(ggpubr)
ggarrange(plarrow, arrows_basin,meanarrow, meanarrowPU, ncol = 2, nrow = 2, labels = "AUTO")


ggsave(filename = "ORDINATIONrestuls.jpeg", width = 10, height = 10, units = "in", dpi = 400)












#-----Temp and Crop Interaction-------
nogreen$TempCat <- NA
nogreen$TempCat[nogreen$temp >= 20] <- "WARM"
nogreen$TempCat[nogreen$temp <  20] <- "COOL"

nogreenCOL$TempCat <- NA
nogreenCOL$TempCat[nogreenCOL$temp >= 20] <- "WARM"
nogreenCOL$TempCat[nogreenCOL$temp <  20] <- "COOL"

summary(ngTEMPCROP.turn)
PseudoR2(ngTEMPCROP.turn, which = "Nagelkerke")

crop.turngraph=nogreen%>%
  ggplot(aes(x=crop, y=Turnover, color=TempCat))+
  geom_point()+
  geom_smooth(method = "glm",se=F, method.args = list(family = "binomial"),size=2)+
  scale_color_manual(values = c("steelblue", "tomato"))+
  labs(x="Crop/Pasture Proportion Within 5km", y="Community Turnover", color="Temperature")+
  annotate("text",x=0.75, y=0.15, label="R^2 = 0.24", size=5)+
  theme_classic()+
  theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color="white"),
        legend.position = "none")

crop.turngraph

PseudoR2(ngTEMPCROP.nat, which = "Nagelkerke")
crop.natgraph=nogreen%>%
  ggplot(aes(x=crop, y=pNat, color=TempCat))+
  geom_point()+
  geom_smooth(method = "glm",se=F, method.args = list(family = "binomial"), size=2)+
  scale_color_manual(values = c("steelblue", "tomato"))+
  theme_classic()+
  labs(x="Crop/Pasture Proportion Within 5km", y="Native Species Persistence", color="Temperature")+
  annotate("text",x=0.75, y=0.15, label="R^2 = 0.15", size=5)+
  theme_classic()+
  theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color="white"),
        legend.position = "none")
crop.natgraph

PseudoR2(ngTEMPCROP.col, which = "Nagelkerke")
crop.colgraph=nogreenCOL%>%
  ggplot(aes(x=crop, y=pCol, color=TempCat))+
  geom_point()+
  geom_smooth(method = "glm",se=F, method.args = list(family = "binomial"), size=2)+
  scale_color_manual(values = c("steelblue", "tomato"))+
  theme_classic()+
  labs(x="Crop/Pasture Proportion Within 5km", y="Colonization Proportion", color="Temperature")+
  annotate("text",x=0.75, y=0.15, label="R^2 = 0.17", size=5)+
  theme_classic()+
  theme(axis.title = element_text(size=16),
        axis.text.x = element_text(size=14, color = "black"),
        axis.title.y=element_blank(),
        axis.text.y=element_text(color="white"),
        legend.position = "none")
crop.colgraph

cropgraphs=ggarrange(crop.turngraph, crop.natgraph, crop.colgraph, nrow = 3, common.legend = T,legend = "right", labels = "AUTO")
annotate_figure(cropgraphs,
                top = text_grob("Community Metrics by Agricultural Land Use and Temperature", color = "black", face = "bold", size = 13))
ggsave(filename = "croptempgraphs.jpeg", dpi=400, width = 7, height=9, units = "in")

ggarrange(imp.turn,lm.turn,crop.turngraph,
          imp.pnat,lm.pnat,crop.natgraph,
          imp.col,lm.pcol,crop.colgraph,
          ncol = 3, nrow = 3,labels = c("A","D","G",
                                        "B","E","H", 
                                        "C","F","I"))

ggsave(filename = "importanceandtempgraph.tiff", dpi = 400, height = 10, width = 13, units = "in")

