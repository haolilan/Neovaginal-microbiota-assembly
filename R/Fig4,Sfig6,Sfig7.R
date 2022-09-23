library(xtable)
library(kableExtra)
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

## prof.mp3_sz295
prof.raw <- read.table(paste0(indir.prof,"/mp3/metaphlan3_vir.profile.merge.reads.763.peritoneum541_sz295.txt"),header=T)
prof.raw  <- column_to_rownames(prof.raw,"ID")
prof.sp2  <- prof.split(prof.raw,"s")
colnames(prof.sp2) <- gsub(".y","",colnames(prof.sp2))

##
setwd(paste0(outdir.prefix,"Transmission_Peritoneal/20220705/"))
outdir <- paste0(outdir.prefix,"Transmission_Peritoneal/20220705/")

### Fig4 POINTPLOT##########################
dist.sum.raw <- read.delim(paste0(paste0(outdir.prefix,"Transmission_Peritoneal/20220705/"),"/dist.tree.sum.raw.txt"))
## 个体内dist 
dist.gv <- dist.sum.raw %>% subset((site1 == "GUT"&site2=="V")|(site1 == "V"&site2=="GUT")|
                                     (site1 == "V"&site2=="V")|
                                     (site1 == "T"&site2=="V")|(site1 == "V"&site2=="T")|
                                     (site1 == "GUT"&site2=="T")|(site1 == "T"&site2=="GUT")) %>% subset(s1 == s2)
dist.gv <- dist.gv %>%mutate(shift = ifelse(value>=0.1,"No","Yes"),site=paste(site1,site2))
dist.gv$site<-gsub("GUT V","G V",gsub("V GUT","G V",dist.gv$site))%>%
  gsub("V T","T V",.)
dist.gv$site<-gsub("GUT T","T G",gsub("T GUT","T G",dist.gv$site))
dist.gv$group1<-paste(dist.gv$site1,dist.gv$g1)
dist.gv$group2<-paste(dist.gv$site2,dist.gv$g2)
##
dist.gv$timepoint1 <- ifelse(dist.gv$groups%in%c("Preop to Preop"),"Preop",NA)
dist.gv$timepoint2 <- ifelse(dist.gv$groups%in%c("Preop to 14D","14D to 14D"),"14D",NA) 
dist.gv$timepoint3 <- ifelse(dist.gv$groups%in%c("Preop to 90D","14D to 90D","90D to 90D"),"90D",NA)
dist.gv$timepoint4 <- ifelse(dist.gv$groups%in%c("Preop to Follow-up","14D to Follow-up",
                                                 "90D to Follow-up","Follow-up to Follow-up"),"Follow-up",NA)
write.csv(dist.gv,paste0(outdir,"dist.GV.VV.samesub.csv"),quote = F)

dist.reshape<-rbind(aggregate(dist.gv$value,by=list(dist.gv$tax,dist.gv$site,dist.gv$s1,dist.gv$timepoint1),min),
                    aggregate(dist.gv$value,by=list(dist.gv$tax,dist.gv$site,dist.gv$s1,dist.gv$timepoint2),min),
                    aggregate(dist.gv$value,by=list(dist.gv$tax,dist.gv$site,dist.gv$s1,dist.gv$timepoint3),min),
                    aggregate(dist.gv$value,by=list(dist.gv$tax,dist.gv$site,dist.gv$s1,dist.gv$timepoint4),min))
colnames(dist.reshape)<-c("tax","site","s1","timepoint","value")
dist.reshape <- dist.reshape %>%mutate(shift = ifelse(value>=0.1,"No","Yes"),y=paste(s1,site))

## 
marker.diffgrs.intrasub <- unique(dist.reshape$tax)
stat.marker.diffgrs.intrasub <- data.frame(xtabs(~tax,dist.reshape))%>% arrange(-Freq)
picked.markers <- stat.marker.diffgrs.intrasub$tax%>%as.vector(.)

## pick.main ##
picked.main <- c("Lactobacillus_crispatus","Lactobacillus_iners","Gardnerella_vaginalis","Atopobium_vaginae",
                 "Prevotella_bivia","Prevotella_timonensis","Prevotella_disiens","Ureaplasma_parvum")
dist.gv.picked <- dist.reshape %>% subset(tax %in% picked.main) %>% mutate(tax = factor(tax,levels=picked.main))
dist.gv.picked$timepoint<-recode_factor(dist.gv.picked$timepoint,"Pre"="Preop","P14D"="14D","P90D"="90D","PFU"="Follow-up")
dist.gv.picked$site%>%unique()
## 
pal.site <- c("#3182bd","#fd8d3c")#pal_d3("category20c")(20)
ggplot(dist.gv.picked, aes(x=timepoint,y=y,color=site))+
  geom_point(aes(shape=shift),size=3)+facet_wrap(~tax,scales = "free_y",nrow = 2)+ 
  theme_pubclean()+ theme(strip.background = element_rect (fill="grey90"),legend.position = "right")+
  scale_color_manual(values = pal.site,guide = guide_legend(reverse=TRUE))+
  scale_shape_manual(values=c(1,16))+labs(x="",y="")#0,15


## Sfig7： pick.sup##################################
picked.sup <- setdiff(picked.markers,picked.main)
#picked.sup <- picked.main
dist.gv.picked <- dist.reshape %>% subset(tax %in% picked.sup) %>% mutate(tax = factor(tax,levels=picked.sup))
dist.gv.picked$timepoint<-recode_factor(dist.gv.picked$timepoint,"Preop"="Pre","14D"="P14D","90D"="P90D","Follow-up"="PFU")
dist.gv.picked$site     <-factor(dist.gv.picked$site,levels = c("V V","G V","T V","T G"))
## order 
count_y = dist.gv.picked%>%group_by(tax,y)%>%summarise(n=n())%>%group_by(tax)%>%summarise(n=n())%>%arrange(-n)
dist.gv.picked <- dist.gv.picked%>%mutate(tax = factor(tax,levels=count_y$tax))
## 
pal.site <- c("#fd8d3c","#3182bd","green","pink")
##
theme_set( 
  theme_pubclean()+theme(text = element_text(size=5),panel.grid = element_blank(),strip.background = element_blank(),strip.text = element_text(color = "black",face="italic",size = 5),legend.position = "right")
)
p1<- ggplot(dist.gv.picked%>%subset(tax%in%count_y$tax[1:5]), aes(x=timepoint,y=y,color=site))+
  geom_point(aes(shape=shift),size=1)+facet_wrap(~tax,ncol = 5,scales = "free_y")+
  scale_color_manual(values = pal.site)+
  scale_shape_manual(values=c(1,16))+labs(x="",y="")
p2<- ggplot(dist.gv.picked%>%subset(tax%in%count_y$tax[6:15]), aes(x=timepoint,y=y,color=site))+
  geom_point(aes(shape=shift),size=1)+facet_wrap(~tax,ncol = 5,scales = "free_y")+ 
  scale_color_manual(values = pal.site)+
  scale_shape_manual(values=c(1,16))+labs(x="",y="")
p3<- ggplot(dist.gv.picked%>%subset(tax%in%count_y$tax[16:25]), aes(x=timepoint,y=y,color=site))+
  geom_point(aes(shape=shift),size=1)+facet_wrap(~tax,ncol = 5,scales = "free_y")+ 
  scale_color_manual(values = pal.site)+
  scale_shape_manual(values=c(1,16))+labs(x="",y="")
p4<- ggplot(dist.gv.picked%>%subset(tax%in%count_y$tax[26:length(count_y$tax)]), aes(x=timepoint,y=y,color=site))+
  geom_point(aes(shape=shift),size=1)+facet_wrap(~tax,ncol = 5,scales = "free_y")+ 
  scale_color_manual(values = pal.site)+
  scale_shape_manual(values=c(1,16))+labs(x="",y="")
ggarrange(p1,p2,p3,p4,
          heights = c(1.2,1.3,1,2.6),ncol=1,nrow = 4,common.legend = T,legend = "bottom")




### Sfig6#######################
dist.sum.raw <- read.delim(paste0(paste0(outdir.prefix,"Transmission_Peritoneal/20220705/"),"/dist.tree.sum.raw.txt"))
## vagina samples
dist.vagina <- dist.sum.raw %>% 
  subset(source %in% pairs.samples & variable %in% pairs.samples) 
## 
stat.intrasub <- data.frame(xtabs(~tax,dist.vagina%>% subset(s1 == s2) ))%>% arrange(-Freq)
wilcox.vv <- compare_means(value~ss,group.by = "tax",dist.vagina,method = "wilcox.test")%>%subset(p<0.05)
picked.markers <- stat.intrasub$tax%>%subset(.%in%wilcox.vv$tax)%>%as.vector()
vv.pm <- picked.markers
##
dist.vagina.picked <- dist.vagina %>% subset(tax %in% picked.markers) %>% mutate(tax = factor(tax,levels=picked.markers))
##
counts <- dist.vagina.picked%>%group_by(tax,ss)%>%summarise(count = n())
dist.sum.draw <- dist.vagina.picked%>% mutate(tax = factor(tax,levels=rev(picked.markers)));
dist.sum.draw$value[dist.sum.draw$value==0] <- 0.01*min(dist.sum.draw$value[dist.sum.draw$value!=0]) 
p1<- ggplot(dist.sum.draw,aes(x=tax,y=value,fill=ss))+geom_boxplot(size=0.3,width = 0.7,alpha=0.8,outlier.size =  0.2)+scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,30),labels=c(0.001,0.01,0.1,1,10,30))+
  geom_text(data=counts,aes(x=tax,y=35,color=ss,label=paste("n",count,sep="=")),size=2,position = position_dodge(width = 0.7))+
  theme_light()+theme(legend.position = "right")+coord_flip()+scale_color_uchicago(palette = "dark")+scale_fill_uchicago()+
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.y = log10(10))+
  labs(x="",y="Normalized phylogenetic distance",fill="")
p1

## GV  
dist.gv <- dist.sum.raw %>% subset((site1 == "GUT"&site2=="V")|(site1 == "V"&site2=="GUT"))
## 
stat.intrasub <- data.frame(xtabs(~tax,dist.gv%>% subset(s1 == s2) ))%>% arrange(-Freq)
wilcox.gv <- compare_means(value~ss,group.by = "tax",dist.gv,method = "wilcox.test")%>%subset(p<0.05)
picked.markers <-  stat.intrasub$tax%>%subset(.%in%wilcox.gv$tax)%>%as.vector()
vg.pm <- picked.markers
##
dist.gv.picked <- dist.gv %>% subset(tax %in% picked.markers) %>% mutate(tax = factor(tax,levels=picked.markers))
## 
counts <- dist.gv.picked%>%group_by(tax,ss)%>%summarise(count = n())
dist.sum.draw <- dist.gv.picked%>% mutate(tax = factor(tax,levels=rev(picked.markers)));
dist.sum.draw$value[dist.sum.draw$value==0] <- 0.01*min(dist.sum.draw$value[dist.sum.draw$value!=0]) 
p2<-ggplot(dist.sum.draw,aes(x=tax,y=value,fill=ss))+geom_boxplot(size=0.3,width = 0.7,alpha=0.8,outlier.size =  0.2)+scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,30),labels=c(0.001,0.01,0.1,1,10,30))+
  geom_text(data=counts,aes(x=tax,y=35,color=ss,label=paste("n",count,sep="=")),size=2,position = position_dodge(width = 0.7))+
  theme_light()+theme(legend.position = "right")+coord_flip()+scale_color_uchicago(palette = "dark")+scale_fill_uchicago()+
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.y = log10(10))+
  labs(x="",y="Normalized phylogenetic distance",fill="")
p2


## VT  
dist.vt <- dist.sum.raw %>% subset((site1 == "T"&site2=="V")|(site1 == "V"&site2=="T"))
stat.intrasub <- data.frame(xtabs(~tax,dist.vt%>% subset(s1 == s2) ))%>% arrange(-Freq)
wilcox.vt <- compare_means(value~ss,group.by = "tax",dist.vt,method = "wilcox.test")%>%subset(p<0.05)
picked.markers <- stat.intrasub$tax%>%subset(.%in%c(vv.pm,vg.pm,wilcox.vt$tax))%>%as.vector()
##
dist.vt.picked <- dist.vt %>% subset(tax %in% picked.markers) 
## 
counts <- dist.vt.picked%>%group_by(tax,ss)%>%summarise(count = n())
dist.sum.draw <- dist.vt.picked%>% mutate(tax = factor(tax,levels=rev(picked.markers)));
dist.sum.draw$value[dist.sum.draw$value==0] <- 0.01*min(dist.sum.draw$value[dist.sum.draw$value!=0]) 
p3<- ggplot(dist.sum.draw,aes(x=tax,y=value,fill=ss))+geom_boxplot(size=0.3,width = 0.7,alpha=0.8,outlier.size =  0.2)+scale_y_log10(breaks=c(0.001,0.01,0.1,1,10,30),labels=c(0.001,0.01,0.1,1,10,30))+
  geom_text(data=counts,aes(x=tax,y=35,color=ss,label=paste("n",count,sep="=")),size=2,position = position_dodge(width = 0.7))+
  theme_light()+theme(legend.position = "right")+coord_flip()+scale_color_uchicago(palette = "dark")+scale_fill_uchicago()+
  stat_compare_means(method = "wilcox.test",label = "p.signif",label.y = log10(10))+
  labs(x="",y="Normalized phylogenetic distance",fill="")
p3


### FIG4 Marker tree#########################################
indir_tree <- paste0(indir.prof,"tree/20220705/")
outdir     <- paste0(outdir,"/tree/")
picked.markers <- stat.marker.diffgrs.intrasub$tax%>%as.vector(.)
phen.meta <- read_excel(paste0(indir.meta,"/enroll.507.xlsx"))
phen.meta$Sites[phen.meta$enroll=="#N/A"]<-"OTHERS"
phen.meta$SubjectID[phen.meta$enroll=="#N/A"]<-"NA"
phen.meta$Sub_site <- paste(phen.meta$SubjectID,phen.meta$Site)
phen.meta$Sub_site[phen.meta$enroll=="#N/A"]<-"NA"
phen.meta <- phen.meta[phen.meta$enroll!="#N/A",]

#### 重新排 ###
picked.markers <-c("Lactobacillus_crispatus","Ureaplasma_parvum","Atopobium_vaginae","Lactobacillus_iners",
                   "Gardnerella_vaginalis","Prevotella_timonensis","Prevotella_bivia","Prevotella_disiens")
cbPa.39 <- c(brewer.pal(8,"Dark2"),brewer.pal(8,"Accent")[c(1:3,5:8)],brewer.pal(8,"Set2"),
             brewer.pal(10,"Paired"),brewer.pal(8,"Pastel2"))

for(sp.names in picked.markers){
  ##
  tree.file <- paste0(indir_tree,"RAxML_bestTree.s__",sp.names,".StrainPhlAn3.tre")   
  outfiles  <- paste0(outdir,sp.names,"_tree.pdf")   
  e.sir.tre <- read.tree(tree.file)
  e.sir.gg <- ggtree(e.sir.tre,layout="rectangular")#branch.length='none'
  e.sir.meta <- phen.meta
  intersect(rownames(e.sir.meta),e.sir.gg$data$label)
  e.sir.gg <- e.sir.gg %<+% e.sir.meta
  e.sir.gg$data$Sites[is.na(e.sir.gg$data$Sites)]<-"OTHERS"
  e.sir.gg$data$Sites<- factor(e.sir.gg$data$Sites,levels = c("preop","14D","90D","follow-up" ,"OTHERS"))
  ###
  p <- e.sir.gg +  
    geom_tippoint(size = 2.5, aes(color = SubjectID,shape=Sites )) +
    aes(branch.length = 'length' ) +
    geom_tiplab(aes(label=Sub_site,col=SubjectID ),hjust = -0.2,size = 2, align=F) +
    scale_shape_manual(values = c(15,16,17,18,1),breaks = c("preop","14D","90D","follow-up" ,"OTHERS"))+
    scale_color_manual(values = cbPa.39,breaks = unique(phen.meta$SubjectID))+
    theme_tree2() + theme(legend.position="bottom")+labs(title = sp.names)
  ggsave(p,paste0(sp.names,".pdf"))
}
