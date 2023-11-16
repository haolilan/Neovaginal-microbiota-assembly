

library(dplyr)
library(readxl)
library(vegan)
library(reshape2)
library(ggplot2)
library(tibble)
library(ggpubr)

## inputs ################################################
setwd("E:/NeovaginaAssembly 2023/")
indir  <- paste0(getwd(), "/data/")
outdir <- paste0(getwd(), "/output/")
#
phen.raw <- read.csv(paste0(indir, "MRKH_study_472_456_202310.csv"))
phen.test <-
  read_excel(paste0(indir, "/mrkh_phen_cc20230719_ALL.xlsx"))
profile <-
  read.table(paste0(indir, "prof.count.s__MRKH.txt"), header = T)
profile_relab <- decostand(profile, method = "total", 2)
##
phen.input <- phen.raw %>% subset(Cohort == "Study_cohort")
prof.s.MRKH <-
  profile_relab %>% .[, intersect(colnames(profile_relab),
                                  (subset(phen.input, Sites == "V" &
                                            Group != "PFU"))$SeqID)] %>%
  .[rowSums(.) != 0,]
#
phen.input.3Y  <- phen.input %>% subset(Timepoints == "P3Y")
phen.input.14D <- phen.input %>% subset(Timepoints == "P14D")
phen.test <-
  phen.test %>% dplyr::select(
    ID,
    `MRKH type`,
    `Mould types within a month at P3Y`:`Lubricating oil use with one month at P3Y`,
    `Sexual activity`:`Sexual pain at P3Y`,
    BV_nugent,
    `Nugent score at P3Y`
  ) %>%
  mutate(SeqID = phen.input.3Y$SeqID[match(ID, phen.input.3Y$SubjectID)],
         SeqID.14 = phen.input.14D$SeqID[match(ID, phen.input.14D$SubjectID)])
phen.test$`Mould types within a month at P3Y` <-
  phen.test$`Mould types within a month at P3Y` %>% recode_factor(., "N" =
                                                                    "0", "P" = "1", "R" = "1")
phen.test$`Condom use with one month at P3Y` <-
  as.factor(phen.test$`Condom use with one month at P3Y`)
phen.test$`Lubricating oil use with one month at P3Y` <-
  as.factor(phen.test$`Lubricating oil use with one month at P3Y`)
phen.test$`Sexual activity` <-
  as.factor(phen.test$`Sexual activity`)
phen.test$BV_nugent <- as.factor(phen.test$BV_nugent)
phen.test$`MRKH type` <- as.factor(phen.test$`MRKH type`)


## PERMANOVA for variable metadata####
##Permutational Multivariate Analysis of Variance (PERMANOVA) was used to assess the statistically significant differences between entire microbial community at P2/4Y and metadata including MRKH type and variable metadata at P2/4Y
metadat <-
  phen.test %>% dplyr::select(-ID, -SeqID.14) %>% subset(!is.na(SeqID))
rownames(metadat) <- NULL
metadat <- metadat %>% column_to_rownames("SeqID")
prof.input <-
  prof.s.MRKH %>% dplyr::select(rownames(metadat)) %>% .[rowSums(.) != 0, ]
prof.bray  <- vegdist(t(prof.input), method = "bray") %>% as.matrix()
individual.perm.fun <-
  function (metadat,
            bray,
            file.prefix,
            n.perm = 999,
            n.proc = 2) {
    print(dim(metadat))
    print(dim(bray))
    mute.low.freq <- function(dat) {
      cutoff <- round(nrow(dat) * 0.01)
      for (i in 1:ncol(dat)) {
        if (length(unique(dat[, i])) < 8 &
            sum(!is.na(dat[, i])) > 0) {
          # categorical data type
          frq <- as.data.frame(table(dat[, i]))
          # print(paste0(i,'the factor before mute low freq: '))
          # print(table(dat[,i]))
          entries <- frq[frq[, 2] < cutoff, 1]
          dat[dat[, i] %in% entries, i] <- NA
          # print(paste0(i,'th factor after mute low freq: '))
          # print(table(dat[,i]))
        }
      }
      return(dat)
    }
    adonis.plot <- function(pe.adonis, file.prefix) {
      pe.adonis.sig <-
        pe.adonis[pe.adonis$Padjust < 0.05 & !is.na(pe.adonis$Padjust), ]
      pe.adonis.sig$R2 <-
        round(as.numeric(pe.adonis.sig$R2), digits = 3)
      pe.adonis.sig <-
        pe.adonis.sig[order(pe.adonis.sig$R2, decreasing = F), ]
      pe.adonis.sig$factor <-
        factor(pe.adonis.sig$factor, levels = pe.adonis.sig$factor)
      f <-
        ggplot(pe.adonis.sig, aes(x = factor, y = R2)) + geom_bar(stat = 'identity') +
        theme_bw() + coord_flip()
      ggsave(paste0(file.prefix, '.pdf'), f)
      f
    }
    metadat <- mute.low.freq(metadat)
    fac <- colnames(metadat)
    
    library(vegan)
    
    out <- matrix(NA, nrow = ncol(metadat), ncol = 10)
    colnames(out) <-
      list(
        "factor",
        "SampleNum",
        "distribution",
        "Df",
        "SumsOfSqs",
        "R2",
        "F.Model",
        "Pr(>F)",
        "R2.adjust",
        "Padjust"
      )
    out <- as.data.frame(out)
    for (i in 1:ncol(metadat)) {
      idx <- !is.na(metadat[, i])
      # valid data length==0 or ( for factor only one valid level)
      print(fac[i])
      if (is.factor(metadat[, i])) {
        levels(metadat[, i]) -> t1
        as.integer(table(metadat[, i])) -> t2
        paste(paste(t1, t2, sep = ':'), collapse = ' - ') -> dis
      } else {
        if (sum(idx) > 0)
          dis <-
            paste(round(summary(metadat[, i]), digits = 2), collapse = ' - ')
      }
      
      if (sum(idx) < 2 * nrow(metadat) * 0.01 |
          length(unique(metadat[idx, i])) == 1) {
        if (sum(idx) > 0) {
          out[i, ] <- c(fac[i], sum(idx), dis, 0, NA, NA, NA, NA, NA, NA, NA)
        } else {
          out[i, ] <- c(fac[i], sum(idx), dis, NA, NA, NA, NA, NA, NA, NA, NA)
        }
      } else {
        samples <- as.character(rownames(metadat)[idx])
        fac.dat <- metadat[idx, i]
        fac.bray <- bray[samples, samples]
        res <-
          adonis2(fac.bray ~ fac.dat,
                  permutations = n.perm,
                  parallel = n.proc)
        tmp <- res[1, ]
        #res <- adonis(fac.bray ~ fac.dat, permutations = n.perm,parallel = n.proc)
        #tmp <- res$aov.tab[1,]
        # print(tmp)
        out[i, ] <-
          c(fac[i],
            sum(idx),
            dis,
            as.character(tmp),
            RsquareAdj(tmp$R2, sum(idx), tmp$Df),
            NA)
      }
      print(paste('done with ', fac[i], sep = ''))
    }
    out$Padjust <- p.adjust(out[, 8], method = 'BH')
    out <- out[order(as.numeric(out$R2), decreasing = T), ]
    write.csv(out,
              paste0(file.prefix, '.csv'),
              quote = F,
              row.names = F)
    f <- adonis.plot(out, file.prefix)
    return(list(f, out))
  }
perm <-
  individual.perm.fun(
    metadat,
    prof.bray,
    paste0(outdir, "/permanova.perm"),
    n.perm = 999,
    n.proc = 2
  )
perm[2]

## species correlated metadata ####
dat.draw <-
  prof.s.MRKH[, na.omit(phen.test$SeqID)] %>% .[rowSums(.) != 0, ]
dat.draw <-
  dat.draw[(apply(dat.draw, 1, mean) > 1e-3)  &
             (apply(dat.draw, 1, function(x) {
               sum(x > 0) / length(x)
             }) > 0.1) ,]
dat.m  <-
  melt(
    dat.draw %>% rownames_to_column("Taxonomy"),
    id.vars = "Taxonomy",
    variable.name = "SeqID",
    value.name = "rel_ab"
  ) %>%
  subset(Taxonomy %in% rownames(dat.draw))
dat.m  <-
  merge(
    dat.m,
    phen.test %>% select(SeqID, `MRKH type`:`Sexual activity`, BV_nugent),
    by = "SeqID",
    all.x = T
  )
colnames(dat.m) <- gsub(" ", "_", colnames(dat.m))

## wilcox
w1 <-
  compare_means(
    data = dat.m %>% subset(!is.na(MRKH_type)),
    rel_ab ~ MRKH_type,
    method = "wilcox.test",
    group.by = c("Taxonomy"),
    p.adjust.method = "fdr"
  ) %>% subset(p.adj < 0.1)
w1
compare_means(
  data = dat.m %>% subset(!is.na(Mould_types_within_a_month_at_P3Y)),
  rel_ab ~ Mould_types_within_a_month_at_P3Y,
  method = "wilcox.test",
  group.by = c("Taxonomy"),
  p.adjust.method = "fdr"
) %>% subset(p.adj < 0.05)
compare_means(
  data = dat.m %>% subset(!is.na(Condom_use_with_one_month_at_P3Y)),
  rel_ab ~ Condom_use_with_one_month_at_P3Y,
  method = "wilcox.test",
  group.by = c("Taxonomy"),
  p.adjust.method = "fdr"
) %>% subset(p.adj < 0.05)
compare_means(
  data = dat.m %>% subset(!is.na(
    Lubricating_oil_use_with_one_month_at_P3Y
  )),
  rel_ab ~ Lubricating_oil_use_with_one_month_at_P3Y,
  method = "wilcox.test",
  group.by = c("Taxonomy"),
  p.adjust.method = "fdr"
) %>% subset(p.adj < 0.05)
compare_means(
  data = dat.m %>% subset(!is.na(Sexual_activity)),
  rel_ab ~ Sexual_activity,
  method = "wilcox.test",
  group.by = c("Taxonomy"),
  p.adjust.method = "fdr"
) %>% subset(p.adj < 0.05)
compare_means(
  data = dat.m %>% subset(!is.na(BV_nugent)),
  rel_ab ~ BV_nugent,
  method = "wilcox.test",
  group.by = c("Taxonomy"),
  p.adjust.method = "fdr"
) %>% subset(p.adj < 0.05)
#write.table(w1,paste0(outdir,"wilcox.phen.q_0.1.txt"),quote = F,sep = "\t")
ggplot(
  dat.m %>% subset(!is.na(MRKH_type) & Taxonomy %in% w1$Taxonomy),
  aes(x = MRKH_type, y = rel_ab, fill = MRKH_type)
) +
  geom_boxplot(alpha = 0.6) + geom_point(position = position_jitterdodge()) +
  facet_wrap( ~ Taxonomy, nrow = 1) + theme_test() + scale_y_log10()
#ggsave(paste0(outdir,"wilcox.phen.q_0.1.pdf"),width = 12,height = 4)

## lm
metadat <-
  phen.test %>% dplyr::select(SeqID,
                              `Female Sexual Function Index at P3Y`:`Sexual pain at P3Y`,
                              `Nugent score at P3Y`)
rownames(metadat) <- NULL
res.lm.phen <- c()
for (j in colnames(metadat)[2:9]) {
  ##
  phen.input <-
    metadat %>% .[, c("SeqID", j)] %>% na.omit(.) %>% column_to_rownames("SeqID")
  prof.input <-
    prof.s.MRKH %>% .[rownames(dat.draw), rownames(phen.input)] %>% .[rowSums(.) !=
                                                                          0, ] %>% t()
  
  res.test <- c()
  pseu.min <- 10 ^ -6
  
  for (i in 1:ncol(prof.input)) {
    data.test <-
      data.frame(grp = phen.input[, j],
                 rel_ab = log10(prof.input[, i] + pseu.min))
    #
    formula <- paste(c("rel_ab ~ ", "grp"), collapse = "")
    m.exx <- lm(as.formula(formula), data = data.test)
    res.test = rbind(res.test,
                     data.frame(
                       taxa = colnames(prof.input)[i],
                       groups = j,
                       #rownames(coef(summary(m.exx)))[2],
                       Estimate = coef(summary(m.exx))[2, 1],
                       p = coef(summary(m.exx))[2, 4]
                     ))
  }
  res.test$q <- res.test[, 4] %>% p.adjust(., method = "fdr")
  write.table(res.test,
              paste0(outdir, "lm.", j, ".txt"),
              quote = F,
              sep = "\t")
  res.lm.phen <- rbind(res.lm.phen, res.test)
}
#write.table(res.lm.phen%>%subset(q<0.1),paste0(outdir,"lm.phen.q_0.1.txt"),quote = F,sep = "\t")
res.lm.phen %>% subset(q < 0.1)
#
dat.lm.m <-
  melt(
    dat.draw %>% rownames_to_column("Taxonomy"),
    id.vars = "Taxonomy",
    variable.name = "SeqID",
    value.name = "rel_ab"
  ) %>%
  subset(Taxonomy %in% rownames(dat.draw))
dat.lm.m  <- merge(dat.lm.m, metadat, by = "SeqID", all.x = T)
ggplot(
  dat.lm.m %>% subset(
    !is.na(`sexual Desire at P3Y`) &
      Taxonomy %in% (res.lm.phen %>% subset(q < 0.1))$taxa
  ),
  aes(x = `sexual Desire at P3Y`, y = rel_ab)
) +
  xlim(1, 5) +
  geom_smooth(method = 'lm') +
  geom_point() + facet_wrap( ~ Taxonomy, nrow = 1) + theme_test() + scale_y_log10()
#ggsave(paste0(outdir,"res.lm.phen.q_0.1.pdf"),width = 4.5,height = 4)
