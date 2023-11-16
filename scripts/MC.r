library(vegan)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(ggpubr)
library(tibble)
library(stringr)

## inputs ####
setwd("E:/NeovaginaAssembly 2023/")
indir  <- paste0(getwd(), "/data/")
outdir <- paste0(getwd(), "/output/")
#
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
profile_relab <- decostand(profile, method = "total", 2)
##
phen.input <-
  phen.raw %>% subset(SeqID %in% intersect(colnames(profile_relab), as.vector(phen.raw$SeqID))) %>%
  subset(Sites == "V" & Group != "PFU") ## exclude clude PFU
prof.input <-
  profile_relab %>% .[, as.vector(phen.input$SeqID)] %>% .[rowSums(.) != 0, ]
# ORDER
phen.input$Group <-
  factor(phen.input$Group,
         levels = c("Pre", "P14D", "P90D", "P3Y", "Hpre", "Hpost"))

##Set  color####
levels(phen.input$Group)
col_Group <-
  c("#FDEF35",
    "#BE8FBE",
    "#FEB24C",
    "#8F1A40",
    "#227E40",
    "#88CDA4")#,"brown1"
names(col_Group) <- levels(phen.input$Group)
#
col_tax_top <-
  c(
    "#238443",
    "#78C679",
    "#045A8D",
    "#0570B0",
    "#74A9CF",
    "#A6BDDB",
    "#D0D1E6",
    "#AE017E",
    "#FA9FB5",
    "#FC4E2A",
    "#FEB24C",
    "#FFEDA0",
    "grey42",
    "#88CDA4",
    "#6FC05D",
    "#00A08A",
    "#4AC0A5",
    "#92C5DE",
    "#6BCCE6",
    "#3B9AB2",
    "#D64E4C",
    "#F98400",
    "#EA3324",
    "#F5BC5B",
    "#EDEB64",
    "#9F267E",
    "#F37D53",
    "#F5E79B",
    "#F7BC85",
    "#E7AEB9",
    "#99B0B1",
    "#F4A582",
    "#FDDBC7",
    "#aa8736",
    "#FF0000",
    "#F2AD00",
    "#E64B35FF"
  )
names(col_tax_top) <-
  c(
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
    "Mycoplasma_hominis",
    "Others",
    "Lactobacillus_jensenii",
    "Lactobacillus_paragasseri",
    "Lactobacillus_gasseri",
    "Lactobacillus_delbrueckii",
    "Prevotella_amnii",
    "Prevotella_corporis",
    "Escherichia_coli",
    "Streptococcus_anginosus_group",
    "Streptococcus_agalactiae",
    "Atopobium_minutum",
    "Enterococcus_faecalis",
    "Sneathia_amnii",
    "Porphyromonas_sp_HMSC065F10",
    "Citrobacter_koseri",
    "Anaerococcus_tetradius",
    "Peptostreptococcus_anaerobius",
    "Haemophilus_parainfluenzae",
    "Klebsiella_aerogenes",
    "Porphyromonas_uenonis",
    "Veillonella_montpellierensis",
    "Finegoldia_magna",
    "Stenotrophomonas_maltophilia",
    "Veillonellaceae_bacterium_DNF00626",
    "Tepidiphilus_thermophilus"
  )

### VENN ####
#### VENN MRKH(PRE,3Y)_Health(MR) ##
venn.groups <- function(prof, group, filenames, cutoff) {
  ## prof : tax & sample
  ## group : colnames(V1,V2) :sampleID, group1
  library(dplyr)
  g <- unique(group$V2)
  Gnum <- length(g)
  number <- group %>% group_by(V2) %>% count(V2)
  gf <- c()
  k <- 1
  for (i in 1:(Gnum)) {
    gf[k] <- sum(number[1:i, 2])
    k <- k + 1
  }
  
  dg <- rownames(prof)
  k <- 1
  for (i in 1:(Gnum)) {
    a <- as.data.frame(apply(prof[, c(k:gf[i])], 1, mean))
    dg <- cbind(dg, a)
    k <- gf[i] + 1
  }
  dg <- dg[, -1]
  colnames(dg) <- g
  dg.lst <- dg
  write.table(dg.lst,
              paste0(outdir, filenames, ".txt"),
              quote = F,
              sep = "\t")
  
  dg <- ifelse(dg > cutoff, rownames(dg), NA)
  dg <- as.data.frame(dg)
  df <- list()
  for (i in 1:length(dg)) {
    a <- na.omit(dg[, i])
    df <- c(df, list(a))
  }
  names(df) <- colnames(dg)
  
  library(VennDiagram)
  library(RColorBrewer)
  venn.diagram(
    df,
    filename = paste0(outdir, filenames),
    height = 6,
    width = 6,
    resolution = 600,
    imagetype = "svg",
    units = "px",
    lwd = 1,
    lty = 1,
    col = "grey90",
    fill = (c(
      "deepskyblue3", "tomato", brewer.pal(12, 'Paired')[c(5, 6, 1, 2)]
    ))[1:length(dg)],
    cat.col = (c(
      "deepskyblue3", "tomato", brewer.pal(12, 'Paired')[c(5, 6, 1, 2)]
    ))[1:length(dg)],
    cex = 1.3,
    cat.cex = 2,
    alpha = 0.8,
    margin = 0.05,
    fontface = 2,
    cat.fontface = 2,
    print.mode = c("raw")
  ) #,"percent"
  
  return(dg.lst)
}
phen.input$group <- recode_factor(
  phen.input$Group,
  "Pre" = "PRE",
  "P14D" = "PM",
  "P90D" = "PM",
  "P3Y" = "P2/4Y",
  "Hpre" = "H_R",
  "Hpost" = "H_M"
) %>%
  factor(., levels = c("PRE", "PM", "P2/4Y", "H_R", "H_M"))
prefix <- "MRKH_PRE_PM_3Y_Health_MR.svg"
groups <-
  dplyr::select(phen.input, SeqID, group) %>% arrange(group)
names(groups) <- c("V1", "V2")
prof.venn.input <- prof.input[, as.vector(groups$V1)]
venn.lst <- venn.groups(prof.venn.input, groups, prefix, 0)

## species accounting with >1‰ mean relative abundance in each of the different stages of MRKH patients and healthy references
prefix <- "MRKH_PRE_PM_3Y_Health_MR_any0.001.svg"
taxa <-
  rownames(venn.lst)[apply(venn.lst[, 1:5], 1, function(x) {
    any(x > 0.001)
  })]
prof.venn.input <- prof.input[taxa, as.vector(groups$V1)]
venn.groups(prof.venn.input, groups, prefix, 0)

###Alpha ####
alpha.stat <- function(profile) {
  library(vegan)
  dataT <- t(profile) ###rm the tax &tax.all
  shannon.wiener = data.frame(values = diversity(dataT, "shannon"))
  simpson = data.frame(values = diversity(dataT, "simpson"))
  s = data.frame(values = specnumber(dataT))
  pielou = data.frame(values = shannon.wiener / log(s))
  ##
  alpha <- data.frame(
    Richness = s$values,
    Shannon = shannon.wiener$values,
    Pielou = pielou$values,
    SeqID = rownames(dataT)
  )
  return(alpha)
}
alpha <- alpha.stat(prof.input)
phen.input$Shannon  <-
  alpha$Shannon [match(phen.input$SeqID, alpha$SeqID)]
phen.input$Richness <-
  alpha$Richness[match(phen.input$SeqID, alpha$SeqID)]
####  Shannon ####
mat <- combn(unique(as.vector(phen.input$Group)), 2)
my_com <- split(mat, col(mat))
temp <-
  compare_means(Shannon ~ Group, phen.input, method = "wilcox.test")
ggplot(phen.input, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(
    alpha = 0.8,
    width = 0.5,
    outlier.size = 0.1,
    size = 0.2
  ) +
  scale_fill_manual(values = col_Group) +
  scale_x_discrete(labels = c("PRE", "P14D", "P90D", "P2/4Y", "H_R", "H_M")) +
  stat_compare_means(method = "kruskal.test", size = 2) +
  stat_compare_means(
    aes(label = ..p.signif..),
    comparisons =  my_com,
    method = "wilcox.test",
    tip.length = 0.01
  ) +
  labs(x = "", y = "Shannon's H") +
  theme_light() + theme(axis.text = element_text(size = 10), panel.grid = element_blank())
ggsave(paste0(outdir, "Sfig.Shannon_group.pdf"),
       width = 6,
       height = 6)

#### Richness ####
mat <- combn(unique(as.vector(phen.input$Group)), 2)
my_com <- split(mat, col(mat))
ggplot(phen.input, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot(
    alpha = 0.8,
    width = 0.5,
    outlier.size = 0.1,
    size = 0.2
  ) +
  scale_fill_manual(values = col_Group) +
  scale_x_discrete(labels = c("PRE", "P14D", "P90D", "P2/4Y", "H_R", "H_M")) +
  stat_compare_means(method = "kruskal.test", size = 2) +
  stat_compare_means(
    aes(label = ..p.signif..),
    comparisons =  my_com,
    method = "wilcox.test",
    tip.length = 0.01
  ) +
  labs(x = "", y = "Observed Species") +
  theme_light() + theme(axis.text = element_text(size = 10), panel.grid = element_blank())
ggsave(paste0(outdir, "Sfig.Richness_group.pdf"),
       width = 6,
       height = 6)

### Beta ####
groups <-
  dplyr::select(phen.input, SeqID, Group)
names(groups) <- c("V1", "V2")
dist   <-
  vegdist(t(prof.input) , method = "bray") %>% as.matrix(.) %>% data.frame()

## bray
dist.lower <- dist

dist.lower[upper.tri(dist.lower)] <- NA
dist.melt.raw  <-
  dist.lower %>% rownames_to_column("source") %>% melt(id = "source", variable.name =
                                                         "variable") %>% subset(!(value %in% c(0, NA)))
dist.melt.raw <- dist.melt.raw %>% ##","g1","g2","s1","s2"
  mutate(
    g1 = phen.input$Group[match(source, phen.input$SeqID)],
    g2 = phen.input$Group[match(variable, phen.input$SeqID)],
    s1 = phen.input$SubjectID[match(source, phen.input$SeqID)],
    s2 = phen.input$SubjectID[match(variable, phen.input$SeqID)]
  ) %>%
  mutate(across.indi = ifelse(s1 == s2, "intra", "inter")) %>%
  mutate(group = paste(g1, g2))

#### Group ####
dist.melt.intra.group <- dist.melt.raw %>% subset(g1 == g2)
mat <- combn(unique(as.vector(dist.melt.intra.group$g1)), 2)
my_com <- split(mat, col(mat))
ggplot(dist.melt.intra.group, aes(x = g1, y = value, fill = g1)) +
  geom_boxplot(
    alpha = 0.8,
    width = 0.5,
    outlier.size = 0.1,
    size = 0.2
  ) +
  scale_fill_manual(values = col_Group, guide = guide_legend(ncol = 1, title = "Group")) +
  scale_x_discrete(labels = c("PRE", "P14D", "P90D", "P2/4Y", "H_R", "H_M")) +
  stat_compare_means(method = "kruskal.test", size = 2) +
  stat_compare_means(
    aes(label = ..p.signif..),
    comparisons =  my_com,
    method = "wilcox.test",
    tip.length = 0.01
  ) +
  labs(x = "", y = "Bray–Curtis dissimilarity") +
  theme_light() + theme(axis.text = element_text(size = 10), panel.grid = element_blank())
ggsave(paste0(outdir, "sfig.bray_intra_group.pdf"),
       width = 6,
       height = 6)

#### Adjacent timepoints  ####
unique(dist.melt.raw$group)
dist.melt <-
  dist.melt.raw %>% subset(group %in% (unique(dist.melt.raw$group) %>% .[grep("Hpre|Hpost", ., invert = T)]) &
                             g1 != g2)
unique(dist.melt$group)
dist.melt <- dist.melt %>%
  mutate(
    group = recode_factor(
      group,
      `Pre P14D` = "Pre/P14D",
      `P14D Pre` = "Pre/P14D",
      `P90D P14D` = "P14D/P90D",
      `P14D P90D` = "P14D/P90D",
      `P90D Pre` = "Pre/P90D",
      `Pre P90D` = "Pre/P90D",
      `P3Y P14D` = "P14D/P3Y",
      `P3Y P90D` = "P90D/P3Y",
      `P3Y Pre`  = "Pre/P3Y"
    )
  )
table(dist.melt$across.indi, dist.melt$group)
dist.melt.test <-
  compare_means(value ~ group,
                dist.melt,
                group.by = "across.indi",
                method = "wilcox.test")
write.table(
  dist.melt.test,
  paste0(outdir, "dist.melt.test.wilcox", ".txt"),
  quote = F,
  sep = "\t",
  row.names = F
)

##
data.dis.draw <- dist.melt %>%
  mutate(groupV2 = group,
         groupV3 = factor(across.indi, levels = c("intra", "inter"))) %>%
  subset(groupV2 %in% c("Pre/P14D",  "P14D/P90D", "P90D/P3Y"))
table(data.dis.draw$groupV3, data.dis.draw$groupV2)
stat.mean <- data.dis.draw %>%
  group_by(groupV2, groupV3) %>%
  summarise(value = mean(value))   
##
compare_means(value ~ groupV2,
              data.dis.draw,
              group.by = "groupV3",
              method = "kruskal.test")
compare_means(value ~ groupV2,
              data.dis.draw,
              group.by = "groupV3",
              method = "wilcox.test")

##
ggplot(data.dis.draw, aes(x = groupV2, y = value, fill = groupV3)) +
  geom_boxplot(
    alpha = 0.8,
    width = 0.5,
    outlier.size = 0.1,
    size = 0.2
  ) +
  geom_point(data = stat.mean,
             aes(group = groupV3, color = groupV3),
             size = 1.5) +
  geom_line(
    data = stat.mean,
    aes(group = groupV3, color = groupV3),
    size = 1,
    linetype = 1
  ) +
  scale_x_discrete(labels = c("PRE vs P14D", "P14D vs P90D", "P90D vs P2/4Y")) +
  stat_compare_means(aes(label = paste0(..p.format..)), method = "wilcox.test", size =
                       2) +
  labs(x = "", y = "Bray–Curtis dissimilarity") +
  theme_light() + theme(axis.text = element_text(size = 10), panel.grid = element_blank())
ggsave(paste0(outdir, "sfig.bray-adjacent.pdf"),
       width = 6,
       height = 3)


#### Pre vs post-surgery####
data.dis.draw <- dist.melt %>%
  mutate(groupV2 = group,
         groupV3 = factor(across.indi, levels = c("intra", "inter"))) %>%
  subset(groupV2 %in% c("Pre/P14D", "Pre/P90D", "Pre/P3Y"))
table(data.dis.draw$groupV3, data.dis.draw$groupV2)
stat.mean <- data.dis.draw %>%
  group_by(groupV2, groupV3) %>%
  summarise(value = mean(value))  
##
compare_means(value ~ groupV2,
              data.dis.draw,
              group.by = "groupV3",
              method = "kruskal.test")
compare_means(value ~ groupV2,
              data.dis.draw,
              group.by = "groupV3",
              method = "wilcox.test")
##
ggplot(data.dis.draw, aes(x = groupV2, y = value, fill = groupV3)) +
  geom_boxplot(
    alpha = 0.8,
    width = 0.5,
    outlier.size = 0.1,
    size = 0.2
  ) +
  geom_point(data = stat.mean,
             aes(group = groupV3, color = groupV3),
             size = 1.5) +
  geom_line(
    data = stat.mean,
    aes(group = groupV3, color = groupV3),
    size = 1,
    linetype = 1
  ) +
  scale_x_discrete(labels = c("PRE vs P14D", "PRE vs P90D", "PRE vs P2/4Y")) +
  stat_compare_means(aes(label = paste0(..p.format..)), method = "wilcox.test", size =
                       2) +
  labs(x = "", y = "Bray–Curtis dissimilarity") +
  theme_light() + theme(axis.text = element_text(size = 10), panel.grid = element_blank())
ggsave(paste0(outdir, "sfig.bray-Pre.pdf"),
       width = 6,
       height = 3)


### Pairwise.adonis ####
groups <-
  dplyr::select(phen.input, SeqID, Group)
names(groups) <- c("V1", "V2")
prof.adonis = adonis(dist ~ V2, data = groups, permutations = 999)
prof.adonis$aov.tab
write.table(prof.adonis$aov.tab,
            paste0(outdir, "fig2c.pcoa.adonis.Group.txt"),
            quote = F)
##
#library(ranacapa) ## R 4.0.2
#groups <- dplyr::select(phen.input,SeqID,Group);names(groups) <- c("V1","V2")
#groups <- arrange(groups,V2)
##prof.input <- t(prof.input[,groups$V1])
#pairwise.adonis <- pairwise_adonis(x=prof.input, factors=groups$V2,sim_method = "bray",
#                                   p_adjust_m = "fdr",reduce = NULL )
#write.table(pairwise.adonis,paste0(outdir,"adonis.paired.txt"),quote = F)

# draw
pairwise.adonis <-
  read.table(paste0(outdir, "adonis.paired.txt"), sep = " ")
pairwise.adonis$x <-
  str_split_fixed(pairwise.adonis$pairs, ".vs.", 2)[, 1]
pairwise.adonis$y <-
  str_split_fixed(pairwise.adonis$pairs, ".vs.", 2)[, 2]
pairwise.adonis$x <-
  factor(pairwise.adonis$x,
         levels = c("Pre", "P14D", "P90D", "PFU", "P3Y", "Hpre", "Hpost"))
pairwise.adonis$y <-
  factor(pairwise.adonis$y,
         levels = c("Pre", "P14D", "P90D", "PFU", "P3Y", "Hpre", "Hpost"))
pairwise.adonis$R2 <- as.numeric(pairwise.adonis$R2)
pairwise.adonis$p.adjusted <- as.numeric(pairwise.adonis$p.adjusted)
pairwise.adonis$p.adj <-
  pairwise.adonis$p.adjusted %>% cut(
    .,
    breaks = c(0, 0.001, 0.05, 1),
    labels = c("<0.001", "<0.05", "ns"),
    right = F
  ) %>% as.factor()
##
ggplot(merge(pairwise.adonis, pairwise.adonis %>% rename(., x = y, y = x), all = T),
       aes(y, x, size = R2 * 100)) +
  geom_point(aes(color = p.adj), alpha = 0.7) +
  scale_y_discrete(labels = c("PRE", "P14D", "P90D", "P2/4Y", "H_R", "H_M")) +
  scale_x_discrete(labels = c("PRE", "P14D", "P90D", "P2/4Y", "H_R", "H_M")) +
  theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  scale_color_manual(values = brewer.pal(8, "Set1")) +
  labs(x = "", y = "")
ggsave(paste0(outdir, "adonis.paired.pdf"),
       width = 4,
       height = 3)

## Wilcox test for Timepoints ####
tax.top.MRKH <-
  rowSums(prof.input) %>% .[order(., decreasing = T)] %>% .[1:20] %>% names

## species ## add health: 只放tax.top.MRKH
dat.draw <- prof.input[tax.top.MRKH, ]
dat.m  <-
  melt(
    dat.draw %>% rownames_to_column("Taxonomy"),
    id.vars = "Taxonomy",
    variable.name = "SeqID",
    value.name = "rel_ab"
  )
dat.m  <-
  merge(dat.m,
        phen.input %>% dplyr::select(SeqID, SubjectID, Timepoints, Group),
        all.x = T)
dat.m$Timepoints <-
  ifelse(is.na(dat.m$Timepoints),
         as.vector(dat.m$Group),
         as.vector(dat.m$Timepoints))
t0 <-
  compare_means(
    data = dat.m,
    rel_ab ~ Timepoints,
    method = "wilcox.test",
    group.by = c("Taxonomy"),
    p.adjust.method = "fdr"
  )
dat.m.avg  <-
  group_by(dat.m, Taxonomy, Timepoints) %>% summarize(rel_ab_avg = mean(rel_ab), .groups =
                                                        "drop_last")
t0 <-
  merge(t0,
        dat.m.avg %>% rename(group2 = Timepoints, mean2 = rel_ab_avg),
        all.x = T) %>% merge(dat.m.avg %>% rename(group1 = Timepoints, mean1 = rel_ab_avg),
                             all.x = T)
t0$group.direct <-
  ifelse(t0$mean1 > t0$mean2,
         as.vector(t0$group1),
         ifelse(t0$mean1 < t0$mean2, as.vector(t0$group2), NA))
write.table(
  t0,
  paste0(outdir, "TIMEPOINTS.sp.wilcox", ".txt"),
  quote = F,
  sep = "\t",
  row.names = F
)

## Vagina type 30% ####
cst_max <- function(profile, cutoff) {
  cst <- data.frame(
    SeqID = colnames(profile),
    cst = ifelse (apply(profile, 2,  max) > cutoff , rownames(profile) [apply(profile, 2, which.max)], "Diverse"),
    max = apply(profile, 2,  max)
  )
  return(cst)
}
cst_s <- cst_max(prof.input, 0.3)
phen.v.MRKH <-
  phen.input %>% subset(Cohort == "Study_cohort") ## exclude PFU
phen.v.MRKH$cluster <-
  cst_s$cst[match(phen.v.MRKH$SeqID, cst_s$SeqID)]
##
data.point <-
  data.frame(xtabs( ~ cluster + Group, phen.v.MRKH)) %>% subset(Freq != 0)
## order
cluster.order <-
  data.point %>% arrange(Group, -Freq) %>% .[, "cluster"] %>% unique() %>% as.vector()

ggplot(data.point, aes(Group, cluster)) +
  geom_point(aes(size = Freq), color = brewer.pal(8, "Set2")[3]) +
  geom_text(aes(label = Freq), hjust = -1.5) +
  geom_text(
    data = data.point %>% group_by(Group) %>% summarise(SamNum = paste0("N=", sum(Freq))),
    aes(label = SamNum),
    y = 0.7
  ) +
  scale_x_discrete(labels = c("PRE", "P14D", "P90D", "P2/4Y")) +
  scale_y_discrete(limits = rev(cluster.order)) +
  scale_color_gradient2() +
  scale_size_area(max_size = 10) +
  labs(x = "", y = "") +
  theme_pubclean()
ggsave(paste0(outdir, "CST_30_reorder.pdf"),
       width = 8,
       height = 8)

## Individual barplots ####
groups <- phen.v.MRKH %>% select(SeqID, Group, SubjectID)
colnames(groups) <- c("V1", "V2", "V3")
#cluster pre
mp_clust = hclust(vegdist(t(prof.input[, groups$V1[groups$V2 == "Pre"]]), method = "bray"), "ward.D2")
label.x = mp_clust$labels[mp_clust$order]
label.sub = (groups$V3[match(label.x, groups$V1)])
##
prof.v.input <-
  prof.input[, phen.v.MRKH$SeqID] ## MRKH V , EXCLUDE PFU
tax.top.MRKH <-
  rowSums(prof.v.input) %>% .[order(., decreasing = T)] %>% .[1:20] %>% names  ## exclude PFU
#
profile.top <- prof.v.input[tax.top.MRKH, phen.v.MRKH$SeqID]
profile.top["Others", ] <-
  colSums(prof.v.input) - colSums(profile.top)
profile.top$Taxonomy <- rownames(profile.top)
##
mp.dat = melt(
  profile.top,
  id.vars =  "Taxonomy",
  variable.name = "SeqID",
  value.name = "rel_ab"
)
mp.dat$Taxonomy = factor(mp.dat$Taxonomy, levels = c(tax.top.MRKH, "Others"))

mp.dat$Group = phen.v.MRKH$Group[match(mp.dat$SeqID, phen.v.MRKH$SeqID)]
mp.dat$Group = mp.dat$Group %>% recode(., "P90D" = "P90D", "P3Y" = "P2/4Y")
mp.dat$SubjectID = phen.v.MRKH$SubjectID[match(mp.dat$SeqID, phen.v.MRKH$SeqID)]
mp.dat$SubjectID = factor(mp.dat$SubjectID, levels = label.sub)
ggplot(mp.dat, aes(x = SubjectID, y = rel_ab, fill = Taxonomy)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(" ") + ylab("Relative Abundance") +
  scale_fill_manual(values = col_tax_top[c(tax.top.MRKH, "Others")],
                    guide = guide_legend(ncol = 6, title = "Taxa")) +
  facet_grid(Group ~ .) +
  theme_test() + theme(
    legend.position = 'bottom',
    axis.text.x = element_text(angle = 90),
    axis.ticks.x = element_blank()
  )
ggsave(paste0(outdir, "Sfig.individual.pdf"),
       width = 12,
       height = 8)


## Mean among body sites ####
#The mean relative abundances of the 10 most abundant neovaginal species
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
profile_relab <- decostand(profile, method = "total", 2)
##
phen.input <-
  phen.raw %>% subset(SeqID %in% intersect(colnames(profile_relab), as.vector(phen.raw$SeqID)) &
                        Cohort == "Study_cohort")
prof.all.input <-
  profile_relab %>% .[, as.vector(phen.input$SeqID)] %>% .[rowSums(.) != 0, ]  ## 39
##
phen.input$group1 <-
  ifelse(phen.input$Sites == "V",
         as.vector(phen.input$Group),
         as.vector(phen.input$Sites))
phen.input$group1 = factor(
  phen.input$group1,
  levels = c(
    "Pre",
    "P14D",
    "P90D",
    "PFU",
    "P3Y",
    "GUT",
    "SK",
    "T",
    "PF",
    "Saline"
  )
)
tax.top.neo2 <-
  rowSums(prof.all.input[, phen.input$SeqID[phen.input$Sites == "V" &
                                              phen.input$Group %in% c("P14D", "P90D", "PFU", "P3Y")]]) %>%
  .[order(., decreasing = T)] %>% .[1:10] %>% names  ## include PFU
#
mp.dat = melt(
  prof.all.input %>% rownames_to_column("Taxonomy"),
  id = "Taxonomy",
  variable.name = "SeqID",
  value.name = "rel_ab"
)

## Pre P14D P90D PFU P3Y Hpre Hpost GUT SK T PF Saline
groups <- phen.input %>% select(SeqID, group1)
colnames(groups) <- c("V1", "V2")
unique(groups$V2)
mp.dat$group1   <- groups$V2[match(mp.dat$SeqID, groups$V1)]
mp.dat.avg  <-
  group_by(mp.dat, Taxonomy, group1) %>% summarize(rel_ab_avg = mean(rel_ab), .groups =
                                                     "drop_last")
profile.avg <-
  dcast(mp.dat.avg, Taxonomy ~ group1) %>% column_to_rownames("Taxonomy")
profile.avg <- profile.avg[rowSums(profile.avg) != 0, ]
#write.table(profile.avg ,paste0(outdir,"prof.s.all.AVERAGE",".txt"),quote = F,sep = "\t")

## mean bar: top 10 of tax.top.neo2 (include PFU) ##
draw.dat <- mp.dat.avg %>% subset(Taxonomy %in% tax.top.neo2)
draw.dat$Taxonomy = factor(draw.dat$Taxonomy, levels = tax.top.neo2)
##
ggplot(draw.dat, aes(x = group1, y = rel_ab_avg, fill = Taxonomy)) +
  geom_bar(stat = "identity", position = "stack") +
  xlab(" ") + ylab("Relative Abundance") + ylim(0, 1) +
  scale_x_discrete(
    labels = c(
      "PRE",
      "P14D",
      "P90D",
      "P6/12M",
      "P2/4Y",
      "Stool",
      "Skin",
      "Tongue",
      "Peritoneal fluid",
      "Saline"
    )
  ) +
  scale_fill_manual(values = col_tax_top[tax.top.neo2],
                    guide = guide_legend(ncol = 1, title = "")) +
  theme_test() + theme(
    legend.position = 'right',
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 90)
  )
ggsave(
  paste0(outdir, "Multi-sites-relative-abundance-neovagina.pdf"),
  width = 8.5,
  height = 6.5
)
