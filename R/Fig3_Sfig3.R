
indir.prof<-"E:/OneDrive - BGI Tech Solutions (Hongkong) Co., Ltd/女性生殖道微生态/MRKH/4.数据产出/Profile/"
indir.meta<-"E:/OneDrive - BGI Tech Solutions (Hongkong) Co., Ltd/女性生殖道微生态/MRKH/4.数据产出/metadata/"
outdir.prefix <- "E:/OneDrive - BGI Tech Solutions (Hongkong) Co., Ltd/女性生殖道微生态/MRKH/5.分析结果/"
## phen
phen    <- read_excel(paste0(indir.meta,"/enroll.507.xlsx"))
phen$SeqID <- gsub("-",".",phen$SeqID)
colnames(phen)[c(2,4)]<-c("Time_guess","Sites")
## 提供 lac信息
phen2    <- read_excel(paste0(indir.meta,"/Sub_group_U_1.xlsx"))
phen2$SeqID <- phen2$SeqID_1
## prof.mp3
prof.raw <- read.table(paste0(indir.prof,"/mp3/metaphlan3_vir.profile.merge.817.dedup.txt"),header=T)
prof.raw  <- column_to_rownames(prof.raw,"ID")
prof.sp  <- prof.split(prof.raw,"s")
colnames(prof.sp) <- gsub(".y","",colnames(prof.sp))
####
phen.va.pre <- phen2[phen2$Time_guess=="1D"&phen2$group=="peritoneal_graft",]
phen.va <- phen %>% subset(Sites=="V")%>%mutate(L.class.pre = phen.va.pre$L.Class[match(SubjectID,phen.va.pre$SubjectID)])
####Sfig3 ########################################################################################################
order <-
  phen.va %>% select(SeqID, SubjectID, Time_guess,L.class.pre
  )%>% 
  mutate(timepoints = factor(Time_guess,levels = c("preop","14D", "90D", "follow-up")),
        SubjectID=as.factor(SubjectID)
  ) %>%
  merge(data.frame(xtabs( ~ SubjectID, .))) %>% mutate(sub.lab = paste0(SubjectID, " (", Freq, ")")) %>%
  arrange(SubjectID, timepoints)
seq <- c()
for (i in 1:length(levels(order$SubjectID))) {
  j <- table(order$SubjectID)[i]
  seq <- c(seq, seq(1, j))
}
order <- order %>% mutate(timepoint = seq#,
                          #SubjectID = factor(order$SubjectID, levels = label.x)
                          ) %>%
  arrange(L.class.pre,-Freq, timepoints) %>% ###reorder by the cluster
  mutate(sub.lab = factor(sub.lab, levels = as.vector(sub.lab[!duplicated(SubjectID)])))

####
select.tax <-c(
      "Lactobacillus_crispatus",
      "Lactobacillus_iners",
      "Prevotella_bivia",
      "Prevotella_timonensis",
      "Prevotella_disiens",
      "Prevotella_buccalis",
      "Prevotella_sp_S7_1_8",
      "Gardnerella_vaginalis",
      "Atopobium_vaginae",
      "Ureaplasma_parvum",
      "Ureaplasma_urealyticum",
      "Mycoplasma_hominis"  )
re.select.tax<- c(
  "Lactobacillus crispatus",
  "Lactobacillus iners",
  "Prevotella bivia",
  "Prevotella timonensis",
  "Prevotella disiens",
  "Prevotella buccalis",
  "Prevotella sp. S7_1_8",
  "Gardnerella vaginalis",
  "Atopobium vaginae",
  "Ureaplasma parvum",
  "Ureaplasma urealyticum",
  "Mycoplasma hominis"  )
dat.draw <-
  prof.sp[select.tax, ]
temp <- data.frame(1 - colSums(dat.draw))
dat.draw <- rbind(dat.draw, t(temp))
rownames(dat.draw)[nrow(dat.draw)] <- "other"
dat.draw$tax <- rownames(dat.draw)
dat.m  <- melt(dat.draw, id.vars = "tax", variable.name = "SeqID")
##
dat.m <- merge(order, dat.m, all.x = T)
#tax.order <-  arrange(data.frame(aggregate(dat.m$value, by = list(dat.m$tax), mean)), -x)
dat.m$tax <- gsub("_"," ",dat.m$tax)
dat.m$tax[dat.m$tax=="Prevotella sp S7 1 8"]<-"Prevotella sp. S7_1_8"
dat.m$tax <- 
  factor(dat.m$tax, levels = rev(c(re.select.tax, "other")))
  #factor(dat.m$tax, levels = c(tax.order$Group.1[-which(tax.order$Group.1 == "other")], "other"))
dat.m$timepoints <- recode_factor(dat.m$timepoints,preop="Pre",`14D`="P14D",`90D`="P90D",`follow-up`="PFU")
label.timepoints <- dat.m %>% group_by(timepoints, SubjectID) %>% filter(row_number() == 1) %>% ungroup()
##
cbPal1<-c(brewer.pal(8,"YlGn")[c(7,5)],
          brewer.pal(9,"PuBu")[rev(c(3,4,5,7,8))],
          brewer.pal(9,"RdPu")[rev(c(4,7))],
          brewer.pal(9,"YlOrRd")[rev(c(2,4,6))],"grey90")%>%rev(.)

p1<- ggplot(dat.m%>%subset(L.class.pre=="Lac"), aes(x = timepoint, y = value, fill = tax)) +
  geom_area(position = "stack") +
  #geom_point(aes(x = timepoint, y = -0.05,color=timepoints), size = 2) +
  geom_text(data=label.timepoints%>%subset(L.class.pre=="Lac"),aes(x = timepoint, y = -0.05,color=timepoints,label=timepoints),hjust=0.35, size = 3) +
  scale_fill_manual(values = cbPal1, guide = guide_legend(ncol = 1,reverse = T,label.theme=element_text(face = "italic"))) +
  scale_color_manual(values = brewer.pal(8, "Set1"), guide = guide_legend(ncol = 1)) +
  facet_wrap( ~ sub.lab, scales = "fixed", ncol = 7) +
  labs(fill = "Taxa",color="Timepionts", size = 10) + labs(x = "", y = "Relative Abundance") +
  theme_light() + theme(
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 18),
    panel.grid = element_blank(),strip.background = element_rect(color = "grey")
  )
p2<-ggplot(dat.m%>%subset(L.class.pre=="no-Lac"), aes(x = timepoint, y = value, fill = tax)) +
  geom_area(position = "stack") +
  #geom_point(aes(x = timepoint, y = -0.05,color=timepoints), size = 2) +
  geom_text(data=label.timepoints%>%subset(L.class.pre=="no-Lac"),aes(x = timepoint, y = -0.05,color=timepoints,label=timepoints),hjust=0.35, size = 3) +
  scale_fill_manual(values = cbPal1, guide = guide_legend(ncol = 1,reverse = T,label.theme=element_text(face = "italic"))) +
  scale_color_manual(values = brewer.pal(8, "Set1"), guide = guide_legend(ncol = 1, )) +
  facet_wrap( ~ sub.lab, scales = "fixed", ncol = 7) +
  labs(fill = "Taxa", color="Timepionts",  size= 10) + labs(x = "", y = "Relative Abundance") +
  theme_light() + theme(
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    text = element_text(size = 18),
    panel.grid = element_blank(),strip.background = element_rect(color = "grey")
  )
ggarrange(p1,p2,common.legend = T,align = "hv",ncol = 1,nrow = 2,legend = "right",labels = c("a", "b") )



#### fig3b-h ################################################
fig3.select.tax<- c(
  "Ureaplasma parvum",
  "Lactobacillus crispatus",
  "Lactobacillus iners",
  "Gardnerella vaginalis",
  "Prevotella bivia",
  "Prevotella timonensis",
  "Atopobium vaginae" )
re.fig3.select.tax<- c(
  "U. parvum",
  "L. crispatus",
  "L. iners",
  "G. vaginalis",
  "P. bivia",
  "P. timonensis",
  "A. vaginae" )
pick.list<- list()
res.data <- c()
for (i in 1:length(fig3.select.tax)){
data<-dat.m%>%subset(tax==fig3.select.tax[i])%>% mutate(reshape=ifelse(value>0,"Presence","Absence"))
## pre abs/pre
id.abs<-(data%>%subset(timepoints=="Pre"&reshape=="Absence"))$SubjectID%>%as.vector()
id.pre<-(data%>%subset(timepoints=="Pre"&reshape=="Presence"))$SubjectID%>%as.vector()
data$sub <- ifelse(data$SubjectID%in%id.abs,"Absence_sub","Presence_sub")%>%factor(.,levels=c("Presence_sub","Absence_sub"))
label.x <- merge(data%>%group_by(timepoints,sub)%>%summarise(n=n()),
      data%>%subset(value>0)%>%group_by(timepoints,sub)%>%
   summarise(n_4=n()),by=c("timepoints","sub"),all = T)%>%mutate(n_x=paste0("(",n_4,"/",n,")"))%>%arrange(sub,timepoints)
label.x$n_x<-gsub("NA","0",label.x$n_x)
data <- merge(data,label.x,by=c("timepoints","sub"),all = T)


p<-ggplot(data,aes(timepoints,value,color=SubjectID))+
  geom_point(size=1)+
  geom_line(aes(group=SubjectID))+
  scale_color_manual(values = cbPa.39,breaks = levels(dat.m$SubjectID),guide=guide_legend(ncol = 2))+
  facet_wrap(~sub,scales = "fixed")+
  labs(x="",y=paste0("Relative Abundance"),title = re.fig3.select.tax[i],caption = paste(label.x$n_x,collapse = " "),color="Patient ID")+
  theme_light()+theme(text = element_text(size=7),
                      axis.text = element_text(size=7),panel.grid = element_blank(),plot.title = element_text(face="italic"),
                      strip.background = element_blank(),strip.text = element_text(color = "black",size = 6))
pick.list[[i]] <- p
res.data <- rbind(res.data,data%>%select(tax,sub,timepoints,value,reshape,SubjectID,n:n_x))
}
ggarrange(NULL,pick.list[[1]],pick.list[[2]],pick.list[[3]],pick.list[[4]],pick.list[[5]],pick.list[[6]],pick.list[[7]],
          ncol = 2,nrow = 4,common.legend = T,legend = "right",align = "hv",
          labels = letters,font.label=list(color="black",size=7,face="bold"))
ggsave(paste0(outdir.prefix, "Fig3.220908.pdf"),
       width = 180,
       height = 210,units = c("mm"))
