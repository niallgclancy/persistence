#---------------------Net Species Changes to Only Native Populations------------------------------------------
NET=read.csv("PersistenceMetrics.csv")
removes=read.csv("REMOVES.csv")
removes=removes[,c(1,3)]
#removes=removes%>%rename("RepetID"="RepeatID")
NET=left_join(NET,removes,by="RepeatID")
NET=subset(NET,NET$GearMiss!="Y")
NET=subset(NET, NET$F_MAUG_<1000)
NET=subset(NET, NET$RepeatID!=234 & NET$RepeatID!=74& NET$RepeatID!=237& NET$RepeatID!=238)

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

traits2=traits[,c(1,3)]
NET=left_join(NETsub,traits2,by="Species")
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
  ylab(label = "Occupancy Change")+
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


ggsave(filename="NativeSpChange.tiff",dpi = 400, width = 10, height = 8, units = "in")

