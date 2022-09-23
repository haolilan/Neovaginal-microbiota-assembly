###bray-distance####
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
##############################
## prof.mp3_sz295
prof.raw <- read.table(paste0(indir.prof,"/mp3/metaphlan3_vir.profile.merge.reads.763.peritoneum541_sz295.txt"),header=T)
prof.raw  <- column_to_rownames(prof.raw,"ID")
prof.sp2  <- prof.split(prof.raw,"s")
colnames(prof.sp2) <- gsub(".y","",colnames(prof.sp2))

###
phen2.s <- phen2[,c("SeqID","group.L")]
phen2.s <- phen2.s%>%subset(SeqID %in% intersect(phen2.s$SeqID,colnames(prof.sp2)))%>%arrange(group.L)
prof.v  <- prof.sp2[,phen2.s$SeqID]

### bray #################################################
groups <- dplyr::select(phen2.s,SeqID,group.L);names(groups) <- c("V1","V2")
dist   <- vegdist(t(prof.v) ,method = "bray") %>% as.matrix(.) %>%data.frame()

## bray β多样性
dist.lower <- dist;
dist.lower[upper.tri(dist.lower)]<-NA 
dist.melt.raw  <- dist.lower %>%rownames_to_column("source")%>%melt(id = "source",variable.name="variable")%>%subset(!(value%in%c(0,NA)))

dist.melt.raw <- dist.melt.raw %>% ##","g1","g2","s1","s2"
  mutate(g1=phen2$group.L[match(source,phen2$SeqID)],
         g2=phen2$group.L[match(variable,phen2$SeqID)],
         s1=phen2$SubjectID[match(source,phen2$SeqID)],
         s2=phen2$SubjectID[match(variable,phen2$SeqID)])%>%
  mutate(across.indi=ifelse(s1==s2,"intra","inter"))%>%
  mutate(group=paste(g1,g2))

### fig 2d ##########################
## 留下所有和"H_nL/H_L"
dist.melt <- dist.melt.raw%>%
  subset(group %in% c("H_L 1D","H_nL 1D" ,"H_L 14D","H_nL 14D",
         "H_L 90D","H_nL 90D","H_L follow-up","H_nL follow-up") )%>%
  mutate(g2 = recode_factor(g2,`1D`="Pre",`14D`="P14D",`90D`="P90D",`follow-up`="PFU"))

data.dis.draw <- dist.melt %>%
  mutate(groupV2=factor(g2,levels = c("Pre","P14D","P90D","PFU")),
         groupV3=factor(g1,levels=c("H_L","H_nL")))

stat.mean <- data.dis.draw %>% 
  group_by(groupV2,groupV3) %>% 
  summarise(value = mean(value))  # 按continent统计平均值
## 
compare_means(value~groupV2,data.dis.draw,group.by = "groupV3",method = "kruskal.test",paired = T)
##
ggplot(data.dis.draw,aes(x=groupV2,y=value,fill=groupV3))+
  geom_boxplot(alpha=0.8,width=0.7,outlier.size = 0.1,size=0.3)+
  geom_point(data=stat.mean,aes(group=groupV3,color=groupV3),size =1.5)+
  geom_line(data=stat.mean,aes(group=groupV3,color=groupV3),size =1,linetype= 1)+
  stat_compare_means(aes(label = paste0( ..p.signif..)),method = "wilcox.test",size=2)+
  labs(x="",y="Bray–Curtis dissimilarity")+theme_light()+theme(axis.text = element_text(size=10),panel.grid = element_blank())


### fig 2g ##################
unique(dist.melt.raw$group)
dist.melt <- dist.melt.raw%>%subset(group%in%(unique(dist.melt.raw$group)%>%.[grep("H",.,invert = T)]) & g1!=g2)
unique(dist.melt$group)
dist.melt <- dist.melt%>%
  mutate(group=recode_factor(group, 
                             `1D 14D`="Pre/P14D",  `14D 1D` ="Pre/P14D",      
                             `90D 14D`="P14D/P90D",  `14D 90D` ="P14D/P90D", 
                             `follow-up 90D`="P90D/PFU",  `90D follow-up`="P90D/PFU",
                             `90D 1D`="Pre/P90D",         `1D 90D`  ="Pre/P90D",   
                             `follow-up 1D`="Pre/PFU",   `1D follow-up`  ="Pre/PFU",  
                             `follow-up 14D`="P14D/PFU", `14D follow-up`="P14D/PFU"))
data.dis.draw <- dist.melt %>%
  mutate(groupV2=group,
         groupV3=factor(across.indi,levels=c("intra","inter"))) %>% 
  subset(groupV2%in% c("Pre/P14D",  "P14D/P90D", "P90D/PFU") )
table(data.dis.draw$groupV3,data.dis.draw$groupV2)
stat.mean <- data.dis.draw%>%
  group_by(groupV2,groupV3) %>% 
  summarise(value = mean(value))  # 按continent统计平均值
## 
compare_means(value~groupV2,data.dis.draw,group.by = "groupV3",method = "kruskal.test",paired = T)
##
ggplot(data.dis.draw,aes(x=groupV2,y=value,fill=groupV3))+
  geom_boxplot(alpha=0.8,width=0.5,outlier.size = 0.1,size=0.2)+
  geom_point(data=stat.mean,aes(group=groupV3,color=groupV3),size =1.5)+
  geom_line(data=stat.mean,aes(group=groupV3,color=groupV3),size =1,linetype= 1)+
  stat_compare_means(aes(label = paste0( ..p.signif..)),method = "wilcox.test",size=2)+
  labs(x="",y="Bray–Curtis dissimilarity")+
  theme_light()+theme(axis.text = element_text(size=10),panel.grid = element_blank())

