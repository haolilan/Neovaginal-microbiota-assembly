library(vegan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(tibble)

### inputs ####
setwd("E:/NeovaginaAssembly 2023/")
indir  <- paste0(getwd(), "/data/")
outdir <- paste0(getwd(), "/output/")

##
phen.raw <- read.csv(paste0(indir, "MRKH_study_472_456_202310.csv"))
profile <-
  merge(
    read.table(paste0(indir, "prof.count.s__MRKH.txt"), header = T) %>% rownames_to_column("ID"),
    read.table(paste0(indir, "prof.count.s__v.health.txt"), header = T) %>%
      rownames_to_column("ID"),
    by = "ID",
    all = T
  ) %>% column_to_rownames("ID")
profile[is.na(profile)] <- 0
##
prefix <- "DMM.allsample.39"
phen.input <-
  phen.raw %>% subset(
    SeqID %in% intersect(colnames(profile), as.vector(phen.raw$SeqID)) &
      Sites == "V" &
      Cohort == "Study_cohort" & Group != "PFU"
  )
dim(phen.input)
profile <-
  profile %>% .[, as.vector(phen.input$SeqID)] %>% .[rowSums(.) != 0, ]  ## 39

# DMM fit ####
library(parallel)
library(DirichletMultinomial)
set.seed(123)
dataT <- t(profile)
##
#HMP.fit.file <- paste0(indir,prefix,".fit")
#cl <- makeCluster(6)
#clusterEvalQ(cl,library(DirichletMultinomial))
#clusterExport(cl,c("dataT"))
#HMP.fit <- parLapply(cl, 1:6, dmn, count=dataT,verbose=TRUE)
#stopCluster(cl)
#save(HMP.fit, file=HMP.fit.file)

# Fit figure
load(paste0(indir, prefix, ".fit"))
HMP.lplc <- sapply(HMP.fit, laplace)
#pdf(paste0(outdir,prefix,".nclusters.pdf"), width = 5,height = 4)
plot(HMP.lplc,
     type = "b",
     xlab = "Number of Dirichlet Components",
     ylab = "Model Fit")
#dev.off()

# Best fit
which.fit <- HMP.fit # Choose the fit model to apply
which.model <- which.min(HMP.lplc)  # k
best <- which.fit[[which.model]]## best model
coll <- mixture(best)
col <- apply(coll, 1, which.max)  ##cluster

# output
res.cluster <- data.frame(V1 = names(col), V2 = as.factor(col))
#write.table(res.cluster ,paste0(outdir,prefix,".k=",which.model,"cluster",".txt"),quote = F,sep = "\t")

## Top contributing species to each cluster ####
taxa.top.driver <- c()
top.n <- 10
for (k in seq(which.model)) {
  d = melt(fitted(best))
  colnames(d) <- c("taxa", "cluster", "value")
  #write.table(d,paste0(outdir,prefix,".k=",which.model,".Feature_component_estimates.txt"),quote=F,row.name=F)
  d <- subset(d, cluster == k) %>%
    arrange(value) %>%
    mutate(taxa = factor(taxa, levels = unique(taxa))) %>%
    filter(abs(value) > quantile(abs(value), (length(unique(
      d$taxa
    )) - top.n) / length(unique(d$taxa))))
  p <- ggplot(d, aes(taxa, value)) + geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("The top contribution to cluster", k))
  #ggsave(paste0(outdir,prefix,".k=",which.model,".Top_contribution",".cluster_",k,".pdf"),p,width = 5,height = 4)
  #print(p)
  taxa.top.driver <- c(taxa.top.driver, as.vector(d$taxa)) %>% unique()
}

## DMM Heatmap ####
profile_relab <- decostand(profile, method = "total", 2)

## heatmap
annotation <- res.cluster %>% dplyr::select(V2)
colnames(annotation) <- c("Cluster")
anno_colors = list(Cluster = brewer.pal(8, "Set1")[1:length(levels(annotation$Cluster))])
names(anno_colors$Cluster) <- levels(annotation$Cluster)

## pre-cluster
dat.draw <-
  profile_relab[order(apply(profile_relab, 1, mean), decreasing = TRUE)[1:30], ]
dat.draw[dat.draw == 0] <- min(dat.draw[dat.draw != 0]) * 0.5

order.clust <- c()
for (k in seq(which.model)) {
  g <-
    pheatmap(log10(dat.draw[, rownames(annotation)[annotation$Cluster == k]]), clustering_method = "ward.D")
  order.clust.temp <-
    data.frame(ID = g$tree_col$labels,
               order.all = g$tree_col$order) %>%
    mutate(cluster = annotation$Cluster[match(ID, rownames(annotation))])
  order.clust <- rbind(order.clust, order.clust.temp)
}
order.clust <- order.clust %>% arrange(cluster, order.all)

## draw
#dev.off()
#pdf(paste(outdir,prefix,".k=",which.model,".heatmap.pdf",sep=""),width = 16, height = 8 , onefile = FALSE)
dat.draw <-
  profile_relab[order(apply(profile_relab, 1, mean), decreasing = TRUE)[1:30],
                order.clust$ID]
dat.draw[dat.draw == 0] <- min(dat.draw[dat.draw != 0]) * 0.5

pheatmap(
  log10(dat.draw),
  cluster_cols = F,
  cluster_rows = T,
  cellwidth = 2,
  cellheight = 12,
  fontsize = 13,
  show_rownames = T,
  show_colnames = F,
  color = colorRampPalette(brewer.pal(9, "YlGnBu"))(100),
  border_color = "white",
  annotation_col = annotation,
  annotation_colors = anno_colors,
  main = prefix
)
#dev.off()

## DMM type transition ####
phen.input$cluster  <-
  res.cluster$V2[match(phen.input$SeqID, res.cluster$V1)]

### 4 timepoints #####
data.point <- data.frame(xtabs( ~ cluster + Group, phen.input))
data <-
  data.frame(xtabs( ~ cluster + Group + SubjectID, phen.input)) %>% subset(Freq !=
                                                                             0)
data.m <-
  merge(data, data, by = "SubjectID", suffixes = c("", ".V2")) %>% subset(Group != Group.V2)
data.l <-
  data.m %>% group_by(cluster, Group, cluster.V2, Group.V2) %>% summarise(count =
                                                                            n())
data.seg <- data.l %>% subset(
  Group == "Pre" & Group.V2 == "P14D" |
    Group == "P14D" & Group.V2 == "P90D" |
    Group == "P90D" &
    Group.V2 == "P3Y"
) %>%
  arrange(Group)
data.seg$Group <-
  factor(data.seg$Group, levels = c("Pre", "P14D", "P90D", "P3Y"))
data.seg$cluster <- as.factor(data.seg$cluster)
ggplot(data.point, aes(Group, cluster)) +
  geom_segment(
    data = data.seg,
    aes(
      x = Group,
      y = cluster,
      xend = Group.V2,
      yend = cluster.V2,
      size = count,
      color = count
    )
  ) +
  geom_point(aes(size = Freq), color = brewer.pal(8, "Set2")[6]) +
  geom_text(
    data = data.point %>% group_by(Group) %>% summarise(SamNum = paste0("N=", sum(Freq))),
    aes(label = SamNum),
    y = 0.5
  ) +
  geom_text(aes(label = Freq)) +
  scale_x_discrete(limits = c("Pre", "P14D", "P90D", "P3Y")) +
  scale_color_gradient2() +
  scale_size_area(max_size = 20) +
  theme_pubclean()
#ggsave(paste0(outdir,prefix,".k=",which.model,".trajectory.pdf"), width = 8.5,height = 6)

### pre-3Y #####
data.point <-
  data.frame(xtabs(
    ~ cluster + Group,
    phen.input %>% subset(SubjectID %in% id.3Y &
                            Group %in% c("Pre", "P3Y"))
  ))
data.seg <- data.l %>% subset(Group == "Pre" & Group.V2 == "P3Y") %>%
  arrange(Group)
data.seg$Group <- factor(data.seg$Group, levels = c("Pre", "P3Y"))
data.seg$cluster <- as.factor(data.seg$cluster)
ggplot(data.point, aes(Group, cluster)) +
  geom_segment(
    data = data.seg,
    aes(
      x = Group,
      y = cluster,
      xend = Group.V2,
      yend = cluster.V2,
      size = count,
      color = count
    )
  ) +
  geom_point(aes(size = Freq), color = brewer.pal(8, "Set2")[6]) +
  geom_text(
    data = data.point %>% group_by(Group) %>% summarise(SamNum = paste0("N=", sum(Freq))),
    aes(label = SamNum),
    y = 0.5
  ) +
  geom_text(aes(label = Freq)) +
  scale_x_discrete(limits = c("Pre", "P3Y")) +
  scale_color_gradient2() +
  scale_size_area(max_size = 10) +
  theme_pubclean()
#ggsave(paste0(outdir,prefix,".k=",which.model,".trajectory.Pre_3Y.pdf"), width = 8,height = 6)

## Shifts on pcoa-log10  ####
#PCoA showing the shifts in the vaginal microbiota from PRE to P2/4Y.
id.3Y <- phen.input$SubjectID[phen.input$Group == "P3Y"]
groups <-
  phen.input %>% subset(Group %in% c("Pre", "P3Y") &
                          SubjectID %in% id.3Y) %>% select(SeqID, Group)
colnames(groups) <- c("V1", "V2")
prof.input <- profile_relab[, groups$V1] %>% .[rowSums(.) != 0, ]

## LOG10-PCOA
data.log <- log10(prof.input)
data.log[data.log == "-Inf"] <- 0
data.log <- log10(prof.input) + 8
data.log[data.log == "-Inf"] <- 0
dist <- vegdist(t(data.log), method = "bray")
#
pcoa <- pcoa(dist, correction = "none", rn = NULL)
##
## plotdata
Axis.1 = pcoa$vectors[, 1]
Axis.2 = pcoa$vectors[, 2]
plotdata <-
  data.frame(rownames(pcoa$vectors), Axis.1, Axis.2, groups$V2)
colnames(plotdata) <- c("sample", "Axis.1", "Axis.2", "group")
plotdata <- arrange(plotdata, group)
plotdata$SubjectID <-
  phen.input$SubjectID[match(plotdata$sample, phen.input$SeqID)]
plotdata$Cluster   <-
  res.cluster$V2[match(plotdata$sample, res.cluster$V1)]
plotdata$Cluster   <-  as.factor(plotdata$Cluster)
## plot vectors
pich = rep(c(19, 17, 15), 20) #pich=rep(c(21:25,1:5),5)
cbPalette <- brewer.pal(8, "Set1")#col_Group
##
pc1 <- floor(pcoa$values$Relative_eig[1] * 100)
pc2 <- floor(pcoa$values$Relative_eig[2] * 100)

## set segment
plotdata <- plotdata %>% arrange(group) ##order

data.line <- t(combn(plotdata$sample, 2)) %>% data.frame()
colnames(data.line) <- c("sample1", "sample2")
data.line <- mutate(
  data.line,
  s1 = plotdata$SubjectID[match(data.line$sample1, plotdata$sample)],
  s2 = plotdata$SubjectID[match(data.line$sample2, plotdata$sample)],
  t1 = plotdata$group[match(data.line$sample1, plotdata$sample)],
  t2 = plotdata$group[match(data.line$sample2, plotdata$sample)],
  Type11 = plotdata$Cluster[match(data.line$sample1, plotdata$sample)],
  Type12 = plotdata$Cluster[match(data.line$sample2, plotdata$sample)]
) %>%
  mutate(t12 = paste(t1, "-", t2),
         Type112 = paste(Type11, "-", Type12)) %>%
  subset(sample1 != sample2 & s1 == s2)
table(data.line$t12)
xtabs( ~ Type112 + t12, data.line)

## 添加坐标
data.line <- mutate(
  data.line,
  x1 = plotdata$Axis.1[match(data.line$sample1, plotdata$sample)],
  y1 = plotdata$Axis.2[match(data.line$sample1, plotdata$sample)],
  x2 = plotdata$Axis.1[match(data.line$sample2, plotdata$sample)],
  y2 = plotdata$Axis.2[match(data.line$sample2, plotdata$sample)]
)
data.line$SubjectID <- data.line$s1
## average
data.line.avg <- data.line %>% group_by(Type112, t12) %>%
  summarize(
    x1 = mean(x1),
    y1 = mean(y1),
    x2 = mean(x2),
    y2 = mean(y2),
    count = n()
  )

names(pich) <- as.vector(unique(plotdata$group))

### cluster.pre
plotdata$Type11 <-
  data.line$Type11[match(plotdata$SubjectID, data.line$s1)]
ggplot(plotdata, aes(Axis.1, Axis.2)) +
  geom_point(aes(
    colour = Cluster,
    fill = Cluster,
    shape = group
  ),
  alpha = 0.8,
  size = 4) +
  geom_segment(
    data = data.line,
    aes(
      x = x1,
      y = y1,
      xend = x2,
      yend = y2,
      linetype = t12
    ),
    alpha = 0.8,
    arrow = arrow(length = unit(0.5, "lines"))
  ) +
  scale_shape_manual(values = pich) +
  scale_colour_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  xlab(paste("PC1 ( ", pc1, "%", " )", sep = "")) +
  ylab(paste("PC2 ( ", pc2, "%", " )", sep = "")) +
  theme(text = element_text(size = 10)) +
  theme(
    panel.background = element_rect(fill = 'white', colour = 'black'),
    panel.grid = element_blank(),
    axis.title = element_text(color = 'black', size = 10),
    axis.ticks.length = unit(0.4, "lines"),
    axis.ticks = element_line(color = 'black'),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(colour = 'black', size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.key = element_blank(),
    legend.key.height = unit(0.8, "cm")
  ) +
  theme(plot.title = element_text(
    size = 20,
    colour = "black",
    hjust = 0.5,
    face = "bold"
  )) +
  facet_wrap( ~ Type11)
#ggsave(paste0(outdir,prefix,".k=",which.model,".sfig.pcoa.cluster.pre.cluster.pdf"),width = 10,height = 4)
