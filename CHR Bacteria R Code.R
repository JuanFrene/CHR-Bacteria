rfPermuteTutorial()


packages <- c('ggcorrplot','ggthemes','dplyr', "ape", "ShortRead", "Biostrings",
              "phyloseq", "corncob","DESeq2", "microbiome", "DECIPHER", "phangorn",
              "tibble", "lme4", "lmerTest", "ggplot2", "vegan", "car", "rcompanion", 'microbiomeSeq',
              "emmeans", "RVAideMemoire",'gplots','plotly','tidyr','VennDiagram','venneuler')
sapply(packages, require, character.only = TRUE)              
                   
# Set working directory
setwd("G:/My Drive/labs/NCSU/Projecto Shaneka 2/Bacteria/")

# Assign variables for imported data
sharedfile = "bacteria.asv.ASV.subsample.shared"
taxfile = "bacteria.asv.ASV.cons.taxonomy"

# Import mothur data
mothur_data <- import_mothur(mothur_shared_file = sharedfile,
                             mothur_constaxonomy_file = taxfile)

colnames(mothur_data@otu_table)
meta <- read.csv("Metadata.csv", header = TRUE, row.names = 1)
sample_dat = sample_data(meta)
sample_names(sample_dat) = rownames(meta)

meta$Replicates=as.factor(meta$Replicates)

physeq1 = merge_phyloseq(mothur_data, sample_dat)
sample_sums(physeq1)
ps.1 = prune_taxa(taxa_sums(physeq1) >= 5, physeq1)
ps.perc <- transform_sample_counts(ps.1, function(x) x / sum(x) * 100)
ps.perc.1 = prune_taxa(taxa_sums(ps.perc) >= 0.01, ps.perc)

# 1. NMDS plot####
paleta_alive <- c("#C200FF",'#FF0000','#8B7500','#00008B',"#FFB919",'#FF7F50',"#00CC7A")

plot_ordination(ps.perc, ordinate(ps.perc, "PCoA", "bray"), color = "State", shape = 'Growth') + theme_few() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  #facet_grid(.~Species, )+
  geom_point(size = 4, aes(shape = Species)) +
  scale_color_manual(values = paleta_alive)

ps.perc.RO = subset_samples(ps.perc, Species ==  "NRO")
plot_ordination(ps.perc.RO, ordinate(ps.perc.RO, "PCoA", "bray"), color = "State", shape = 'Growth') + theme_few() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  geom_point(size = 4) +
  scale_color_manual(values = paleta_alive)

ps.perc.BW = subset_samples(ps.perc, Species !=  "NRO")
plot_ordination(ps.perc.BW, ordinate(ps.perc.BW, "PCoA", "bray"), color = "State", shape = 'Growth') + theme_few() +
  guides(colour = guide_legend(override.aes = list(size = 5))) +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  geom_point(size = 4) +
  scale_color_manual(values = paleta_alive)

# 2. PERMANOVA
# Calculate bray curtis distance matrix
ps.perc_bray <- phyloseq::distance(ps.perc, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.perc))

# Adonis test
PERMANOVA = adonis2(ps.perc_bray ~ Species*State*Unknow, data = sampledf)

#                     Df   SumOfSqs         R2         F Pr..F.
#Species               1  3.6367566 0.10783422 13.506969  0.001
#State                 2  3.1775533 0.09421829  5.900740  0.001
#Unknow                1  0.3810370 0.01129821  1.415177  0.064
#Species:State         2  2.5532745 0.07570767  4.741450  0.001
#Species:Unknow        1  0.3850165 0.01141620  1.429957  0.071
#State:Unknow          2  0.5833161 0.01729603  1.083222  0.270
#Species:State:Unknow  2  0.6607046 0.01959069  1.226933  0.110
#Residual             83 22.3477822 0.66263870        NA     NA
#Total                94 33.7254409 1.00000000        NA     NA


PERMANOVA2 = data.frame(PERMANOVA)
variable  = c('variable','variable','variable','variable','variable','variable','variable','variable','variable')
PERMANOVA2$Significance <- "No Significant"
pval_thres <- 0.05
PERMANOVA2$Significance[which(PERMANOVA2$Pr..F. < pval_thres)] <- "Significant"
PERMANOVA2$Significance <- PERMANOVA2$Significance %>% factor

LC = data.frame(cbind(variable, row.names(PERMANOVA2), PERMANOVA2[,c(3,6)]))
colnames(LC) = c('variable','Effect', 'R2','Significance')
LC$Effect[which(LC$Effect == 'Unknow')] <- "Conditions"
LC$Effect[which(LC$Effect == 'Species:Unknow')] <- "Species:Conditions"
LC$Effect[which(LC$Effect == 'Location:Unknow')] <- "Location:Conditions"
LC$Effect[which(LC$Effect == 'Species:Location:Unknow')] <- "Species:Location:Conditions"


#Plot Variance
dim(LC)
colnames(LC)
dim(LC)
data_melt <- LC[-(9),]

#Reorder x axis: the order of the factor will be the same as in the data.csv file
data_melt$Effect <- as.character(data_melt$Effect)#Turn your 'treatment' column into a character vector
data_melt$Effect <- factor(data_melt$Effect, levels=unique(rev(data_melt$Effect)))#Then turn it back into a factor with the levels in the correct order

mypal = c('white','#c0047b','#eaabcd','#a9d11c','#528501','#2d4800','#dcda65','#53c0cc','#1f4f58')

ggplot(data=data_melt, aes(x=variable, y=R2,fill=Effect)) +
  geom_bar(stat="identity",aes(color=Significance), size = 0.3,width = 1,linewidth=0.4) +# 
  scale_fill_manual(values = mypal) +
  theme_few() + #guides() + #color = FALSE, fill=FALSE
  labs(y="Variance explained (%)") +
  scale_color_manual(values = c('black',"red"),na.value =  "transparent",name = "Significance vs Water")+
  theme(legend.title=element_blank(), legend.margin=margin(c(0,0,0,0)),
        axis.text.x   = element_blank())

ggsave("Bacteria_variance.pdf", height=5, width=3.5, units='in')


# 3. Alfadiversity
# let's calcuolate Shannon diversity and Richness (Observed)
alpha.div <- estimate_richness(ps.1, measures=c("Shannon", "Observed",'Chao', 'Simpson'))
even <- evenness(ps.1, 'pielou')
write.table(alpha.div, "Alfa bact.txt")

plot_richness(ps.1, x = "Species", color = 'Location',
              title = 'Alphadiversity', scales = "free_y", nrow = 1,
              measures = c("Shannon"), sortby = NULL)

Metadata <- ps.1@sam_data

# it is easiest to append the alpha-diversity estimates to the metadata file for downstream analyses
Metadata$Shannon <- paste(alpha.div$Shannon)
Metadata$Observed <- paste(alpha.div$Observed)
Metadata$Simpson <- paste(alpha.div$Simpson)
Metadata$Chao1 <- paste(alpha.div$Chao1)
Metadata$Observed <- as.numeric(Metadata$Observed)
Metadata$Shannon <- as.numeric(Metadata$Shannon)
Metadata$Simpson <- as.numeric(Metadata$Simpson)
Metadata$Chao1 <- as.numeric(Metadata$Chao1)
Metadata$Location <- as.factor(Metadata$Location)
Metadata$Species <- as.factor(Metadata$Species)
Metadata$Replicates <- as.factor(Metadata$Replicates)
Metadata$Field <- as.factor(Metadata$Field)

#write.table(Metadata, "Alfadiversity Bact.txt")

Metadata2 = data.frame(Metadata)

aovShannon <- aov(Shannon~Species*State*Unknow, data=Metadata2)
summary(aovShannon)

lsdpH <- LSD.test(aovShannon, c('State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH
#,'Species','State''Unknow'

#                          Df Sum Sq Mean Sq F value   Pr(>F)    
#Species               1  3.254   3.254  57.219 4.76e-11 ***
#State                 2  0.464   0.232   4.076  0.02048 *  
#Unknow                1  0.280   0.280   4.916  0.02934 *  
#Species:State         2  0.287   0.144   2.524  0.08625 .  
#Species:Unknow        1  0.068   0.068   1.189  0.27868    
#State:Unknow          2  0.667   0.333   5.862  0.00416 ** 
#Species:State:Unknow  2  0.656   0.328   5.770  0.00451 ** 
#Residuals            83  4.720   0.057  


lsdpH <- LSD.test(aovpH_H2O, c('Plant_state','Species','State'), alpha = 0.05, p.adj="bonferroni", group=TRUE) 
lsdpH


Alpha_SP = merge(Table, Metadata2[,-c(2:8)], by = 'ID')
Alpha_SP2 = Alpha_SP[Alpha_SP$ID!='B212',]

summary(lm(Dbh~Observed.y, data=Alpha_Soil2[Alpha_Soil2$Species == "NRO",]))
summary(lm(Dbh~Observed.y, data=Alpha_Soil2[Alpha_Soil2$Species == "BW",]))


ggplot(Alpha_SP2, aes(x=Plant_state, y=Shannon.y)) + 
  geom_point() + 
  facet_grid(State~Species,space = "fixed",scales = "free") +
  geom_smooth(method="loess", se=T,aes(col=Species))+
  theme_few()

###Table ALphadiv
##Figure Alfa Diversidad Compartment
Metadata$Unknow = factor(Metadata$Unknow, c('G','B'))
safe_colorblind_palette <- c("#DDCC77","#CC6677")
ggplot(Metadata, aes(Species, Shannon,  fill = Unknow) ) +
  geom_boxplot() +
  theme_bw()+
  theme(legend.position="none")+
  ylim(c(4,6))+
  facet_grid(State~., scales = "free", space = "free")+
  geom_jitter(show.legend = FALSE, width = 0.25, shape=21, color='black')+
  scale_fill_manual(values = safe_colorblind_palette) +
  scale_color_manual(values = safe_colorblind_palette)

Metadata2 = data.frame(Metadata)

summary(shannon.model<-lmer(Shannon ~ State*Species*Unknow  + (1|Replicates), data = Metadata2))
anova(shannon.model)
emmeans(shannon.model, pairwise ~ State*Species, adjust = "none")

# 4. Plotting Relative Abundance Bar Charts####
# phylum-level
ps.compositional <- microbiome::transform(ps.1, "compositional")
ps.phyla.perc <-taxa_level(ps.compositional, "Rank6")

# identify the 10 most abundant phylum
phylum.10 <- names(sort(taxa_sums(ps.phyla.perc), TRUE)[1:10])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.phylum.10 <- prune_taxa(phylum.10, ps.phyla.perc)

melt.phylum <- psmelt(ps.phylum.10)

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

melt.phylum_mean <- melt.phylum%>%group_by(Species, State, Unknow, OTU)%>%
  summarise_all(mean)

melt.phylum_mean$Unknow = factor(melt.phylum_mean$Unknow, c('G','B'))

ggplot(melt.phylum_mean, aes(x = Unknow, y = Abundance, fill = OTU)) + 
  theme_bw() +
  facet_grid(.~Species*State, scales = "free",space = "fixed")+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = safe_colorblind_palette) + 
  labs(fill = "Phylum")

melt.phylum$Sample= factor(melt.phylum$Sample, levels=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,49,50,51,52,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72))
sample <- c('BW/Autumn/0-10',	'RO/Autumn/0-10',	'BS/Autumn/0-10','BW/Spring/0-10',	'RO/Spring/0-10',	'BS/Spring/0-10',	'BW/Autumn/10-20',	'RO/Autumn/10-20',	'BS/Autumn/10-20',	'BW/Spring/10-20',	'RO/Spring/10-20',	'BS/Spring/10-20')


library(ANCOMBC)
library(DT)

tse2 = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.1)
print(tse2)

out2 = ancombc(data = tse2, assay_name = "counts", 
               tax_level = "Rank6", phyloseq = NULL, 
               formula = "Unknow", 
               p_adj_method = "holm", prv_cut = 0.05, lib_cut = 1000, 
               group = "Unknow", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
               n_cl = 1, verbose = TRUE)

res = out2$res
res_global = out2$res_global

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "RO_BW")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

samp_frac = out2$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out2$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")

df_lfc2=merge(tab_lfc,tab_q[,c(1,3)], by='Taxon' )
colnames(df_lfc2) = c( 'taxon_id', '(Intercept)', ' PlantRO', 'pvalue')
head(df_lfc2)
df_lfc3 = df_lfc2[df_lfc2$pvalue<0.05, ] 
df_lfc4 = df_lfc3[df_lfc3$pvalue>0, ] 

nrow(df_lfc3)

df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(col_name[-1], "SE")
head(df_se)



df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, SpeciesNRO, RO_BWSE) %>%
  dplyr::filter(RO_BWSE != 0) %>% 
  dplyr::arrange(desc(SpeciesNRO)) %>%
  dplyr::mutate(direct = ifelse(SpeciesNRO > 0, "Q. Rubra", "J. Nigra"))

df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id)
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Q. Rubra", "J. Nigra"))

p_bac = ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = SpeciesNRO, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = SpeciesNRO - RO_BWSE, ymax = SpeciesNRO + RO_BWSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_bac

#5. Heatmap
ps.genus.perc.bac <-taxa_level(ps.perc, "Rank6")

# identify the 50 most abundant genus
genus.10.bac <- names(sort(taxa_sums(ps.genus.perc.bac), TRUE)[2:51])
# now we can use the list of the top phylum and subset those from the phyloseq object
ps.genus.10.bac <- prune_taxa(genus.10.bac, ps.genus.perc.bac)

Table.bac = cbind(ps.genus.10.bac@sam_data, log10(ps.genus.10.bac@otu_table+1)) 
Table_mean.bac <- Table.bac%>%group_by(Species, State, Unknow)%>%
  summarise_all(mean)

Mean1 <- data.frame(Table_mean.bac)
Mean0_15<- data.frame(t(Mean1[,-c(1:8)]))
names (Mean0_15) = c('BW-Michigan-B','BW-Michigan-G','BW-South_Indiana-B','BW-South_Indiana-G',
                     'BW-West_Indiana-B',"BW-West_Indiana-G",'NRO-Michigan-B','NRO-Michigan-G','NRO-South_Indiana-B','NRO-South_Indiana-G',
                     'NROW-West_Indiana-B',"NRO-West_Indiana-G")
annotation_col = data.frame(
  Genera =substr(colnames(Mean0_15),1,3))
rownames(annotation_col)=colnames(Mean0_15)
colnames(Mean0_15)
row.names(Mean1$ID)


pheatmap(Mean0_15, show_rownames=TRUE,show_colnames=TRUE,
         annotation_col=annotation_col, annotation_legend = TRUE,
         scale = "none",clustering_method="ward.D2",
         fontsize_row = 6, #breaks = 0:4,
         clustering_distance_cols="euclidean")


##6.  Mantel Tests ####
ps.T3 = subset_samples(ps.perc, ID !=  "B212") #
ps.RO = subset_samples(ps.perc, Species ==  "NRO") #
ps.BW = subset_samples(ps.T3, Species !=  "NRO") #

asv.table.BW <- data.frame(otu_table(ps.BW))
asv.table.RO <- data.frame(otu_table(ps.RO))

sam_dataBW <- data.frame(sam_data(ps.BW))
sam_dataRO <- data.frame(sam_data(ps.RO))

xdist.prod.BW <- vegdist(t(asv.table.BW), method = "bray")
xdist.prod.RO <- vegdist(t(asv.table.RO), method = "bray")
xdist.prod.BW_v <- data.frame(as.vector(xdist.prod.BW))
xdist.prod.RO_v <-as.vector(xdist.prod.RO)

BW = c(BW,1128)

nrow(xdist.prod.BW_v)


sam_data$Mg <- as.numeric(sam_data$Mg)
ydist.Mg <- vegdist(sam_data$Mg, method = "euclid")
#Mantel statistic r: 0.1155  Significance: 0.003 
#RO: Mantel statistic r: 0.2652 Significance: 0.001
#BW: Mantel statistic r: 0.0623 Significance: 0.091 

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Mg, method = "spearman")

sam_data$Ca <- as.numeric(sam_data$Ca)
ydist.Ca <- vegdist(sam_data$Ca, method = "euclid")
#Mantel statistic r: 0.1998  Significance: 0.001 
#RO: Mantel statistic r: 0.3142 Significance: 0.001 
#BW: Mantel statistic r: 0.0398 Significance: 0.192 

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Ca, method = "spearman")

sam_data$K <- as.numeric(sam_data$K)
ydist.K <- vegdist(sam_data$K, method = "euclid")
#Mantel statistic r: 0.01915  Significance: 0.29
#RO: Mantel statistic r: 0.1188 Significance: 0.0024 
#BW: Mantel statistic r: 0.2698 Significance: 0.001 

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.K, method = "spearman")

sam_data$SumNH4 <- as.numeric(sam_data$SumNH4)
ydist.SumNH4 <- vegdist(sam_data$SumNH4, method = "euclid")
#Mantel statistic r: 0.1827  Significance: 0.001
#RO: Mantel statistic r: 0.3169 Significance: 0.001 
#BW: Mantel statistic r: 0.0394 Significance: 0.19
# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.SumNH4, method = "spearman")

sam_data$CEC <- as.numeric(sam_data$CEC)
ydist.CEC <- vegdist(sam_data$CEC, method = "euclid")
#Mantel statistic r: 0.1436  Significance: 0.003
#RO: Mantel statistic r: 0.282 Significance: 0.001 
#BW: Mantel statistic r: 0.04865 Significance: 0.147

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.CEC, method = "spearman")

sam_data$Base_Sat <- as.numeric(sam_data$Base_Sat)
ydist.Base_Sat <- vegdist(sam_data$Base_Sat, method = "euclid")
#Mantel statistic r: 0.274  Significance: 0.001
#RO: Mantel statistic r: 0.2481 Significance: 0.002 
#BW: Mantel statistic r: 0.1068 Significance: 0.014

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Base_Sat, method = "spearman")

sam_data$Corg <- as.numeric(sam_data$Corg)
ydist.Corg <- vegdist(sam_data$Corg, method = "euclid")
#Mantel statistic r: 0.01913  Significance: 0.289
#RO: Mantel statistic r: 0.03801 Significance: 0.283 
#BW: Mantel statistic r: 0.01214 Significance: 0.393

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Corg, method = "spearman")

sam_data$pH_CaCl2 <- as.numeric(sam_data$pH_CaCl2)
ydist.pH_CaCl2 <- vegdist(sam_data$pH_CaCl2, method = "euclid")
#Mantel statistic r: 0.2363  Significance: 0.001
#RO: Mantel statistic r: 0.1975 Significance: 0.001 
#BW: Mantel statistic r: 0.05826 Significance: 0.116

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.pH_CaCl2, method = "spearman")

sam_data$pH_H2O <- as.numeric(sam_data$pH_H2O)
ydist.pH_H2O <- vegdist(sam_data$pH_H2O, method = "euclid")
#Mantel statistic r: 0.2467  Significance: 0.001
#RO: Mantel statistic r: 0.1759 Significance: 0.001 
#BW: Mantel statistic r: 0.06487 Significance: 0.084

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.pH_H2O, method = "spearman")


sam_data$PBray <- as.numeric(sam_data$PBray)
ydist.PBray <- vegdist(sam_data$PBray, method = "euclid")
#Mantel statistic r: -0.001031  Significance: 0.469
#RO: Mantel statistic r: 0.06074 Significance: 0.153 
#BW: Mantel statistic r: 0.07685 Significance: 0.046

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.PBray, method = "spearman")

sam_data$Clay <- as.numeric(sam_data$Clay)
ydist.Clay <- vegdist(sam_data$Clay, method = "euclid")
#Mantel statistic r: 0.2219  Significance: 0.001
#RO: Mantel statistic r: 0.1993 Significance: 0.003 
#BW: Mantel statistic r: 0.122 Significance: 0.01

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Clay, method = "spearman")


sam_data$Silt <- as.numeric(sam_data$Silt)
ydist.Silt <- vegdist(sam_data$Silt, method = "euclid")
#Mantel statistic r: 0.2756  Significance: 0.001
#RO: Mantel statistic r: 0.3376 Significance: 0.001 
#BW: Mantel statistic r: 0.5803 Significance: 0.001

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Silt, method = "spearman")

sam_data$Sand <- as.numeric(sam_data$Sand)
ydist.Sand <- vegdist(sam_data$Sand, method = "euclid")
#Mantel statistic r: 0.2917  Significance: 0.001
#RO: Mantel statistic r: 0.3632 Significance: 0.001 
#BW: Mantel statistic r: 0.5267 Significance: 0.001

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod, ydist.Sand, method = "spearman")

sam_dataRO$Dbh <- as.numeric(sam_dataBW$Dbh)
ydist.Dbh.BW <- vegdist(sam_dataBW$Dbh, method = "euclid")
#Mantel statistic r: 0.0635  Significance: 0.011
#BW: Mantel statistic r: 0.1411 Significance: 0.002
#RO: Mantel statistic r: 0.0779 Significance: 0.028

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod.BW, ydist.Dbh.BW, method = "spearman")

sam_dataRO$Height <- as.numeric(sam_dataRO$Height)
ydist.HeightBW <- vegdist(sam_dataBW$Height, method = "euclid")
ydist.Height_v = as.vector(ydist.HeightBW)

#Mantel statistic r: 0.08379  Significance: 0.011
#BW: Mantel statistic r: 0.2038 Significance: 0.001
#RO: Mantel statistic r: 0.05363 Significance: 0.179

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod.BW, ydist.HeightBW, method = "spearman")


Trait_Root = data.frame(cbind(xdist.prod.BW_v, ydist.Height_v))
Trait_Root$ydist.Height_v <- as.numeric(Trait_Root$ydist.Height_v)
Trait_Root$ as.vector.xdist.prod.BW.  <- as.numeric(Trait_Root$ as.vector.xdist.prod.BW. )

ggplot(Trait_Root, aes(x= as.vector.xdist.prod.BW., y=ydist.Height_v))+
  geom_point() + 
  geom_smooth(method="lm", se=F, colour='black')+
  labs(y='Height', x='Microbiome')+
  theme_few()
ncol(sam_data)

sam_dataBW$Height <- as.numeric(sam_dataBW$Height)
ydist.HeightBW <- vegdist(sam_dataBW$Height, method = "euclid")
ydist.Height_v = as.vector(ydist.HeightBW)

#Mantel statistic r: 0.08379  Significance: 0.011
#BW: Mantel statistic r: 0.2038 Significance: 0.001
#RO: Mantel statistic r: 0.05363 Significance: 0.179

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod.BW, ydist.HeightBW, method = "spearman")


Trait_Root = data.frame(cbind(xdist.prod.RO_v, ydist.ScompoRO))
Trait_Root$ydist.ScompoRO <- as.numeric(Trait_Root$ydist.ScompoRO)
Trait_Root$ xdist.prod.RO_v  <- as.numeric(Trait_Root$xdist.prod.RO_v)

ggplot(Trait_Root, aes(xdist.prod.RO_v, ydist.ScompoRO))+
  geom_point() + 
  geom_smooth(method="lm", se=F, colour='black')+
  labs(x='Soil characteristic', y='Microbiome')+
  theme_few()
ncol(sam_data)

TablaRO2 = data.frame(lapply(TablaRO[,19:32], as.numeric))
TablaBW2 = data.frame(lapply(TablaBW[-c(11),19:32], as.numeric))


ydist.ScompoRO <- vegdist(TablaRO2, method = "euclid")
ydist.ScompoBW <- vegdist(TablaBW2, method = "euclid")
#RO: Mantel statistic r: 0.3746 Significance: 0.001 
#BW: Mantel statistic r: 0.5525 Significance: 0.001

# mantel test with bray-curtis dissimilarity
vegan::mantel(xdist.prod.BW, ydist.ScompoBW, method = "spearman")
vegan::mantel(xdist.prod.RO, ydist.ScompoRO, method = "spearman")



#### log2foldchange
library(reshape2)
ps.perc
tableSyncom = cbind(ps.perc@sam_data[,c(3,6,7)], t(ps.perc@otu_table))
tableSyncom2=tableSyncom

library(reshape2)
melted_tableSyncom <- tableSyncom2 %>% melt

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}


log2foldchange <- c()

for(plant in melted_tableSyncom$Species  %>% unique){
  melted_sub <- melted_tableSyncom %>% subset(Species ==  plant) %>% droplevels
  for(ASV in melted_sub$variable  %>% unique){
    melted_sub2 <- melted_sub %>% subset(variable ==  ASV) %>% droplevels
    
    Good = data.frame(t(melted_sub2[melted_sub2$Growth=='G',][,5]))
    Bad = data.frame(t(melted_sub2[melted_sub2$Growth=='B',][,5]))
    
    Good_mean = apply(Good, 1, mean) 
    Bad_mean = apply(Bad, 1, mean) 
    
    Good_Bad <- log2(Good_mean+1) - log2(Bad_mean+1) 
    
    Good_Bad_statistic = t.test(Good,Bad)
    pvalue = Good_Bad_statistic$p.value#Water_mean+1
    
    result = cbind(Good_Bad,pvalue)
    result2 = cbind(ASV,result)
    result3 = cbind(plant,result2)
    log2foldchange <- rbind(log2foldchange,result3)
  }}

colnames(log2foldchange)=c('Species','ASV','diff','pvalue')
log2foldchange2 = data.frame(log2foldchange)

####Adjust the p-value
log2foldchange2$Significance <- "No Significant"
pval_thres <- 0.05
log2foldchange2$Significance[which(log2foldchange2$pvalue < pval_thres)] <- "q < 0.05"
log2foldchange2$Significance <- log2foldchange2$Significance %>% factor

#log2foldchange2_selected =log2foldchange2[log2foldchange2$Significance == "q < 0.05",] 
nrow(log2foldchange2_selected)

write.table(log2foldchange2_selected_taxa, "Comparasion good bad.txt")

tax_table = cbind(rownames(ps.perc@tax_table), ps.perc@tax_table)
head(tax_table)
colnames(tax_table) = c("ASV","Kingdom", "Phylum", "Classes", "Order", "Family","Genus")

log2foldchange2_selected_taxa = merge(log2foldchange2_selected, tax_table, by='ASV')
#log2foldchange2_selected_taxa[log2foldchange2_selected_taxa$Species=='NRO',]$diff>0,
#88 taxas in BW y 122


Good_Bad = read.table("Comparasion good bad.txt")

library(ggdendro)
###Clusterization
display <- Good_Bad  %>%
  acast(formula = Species~ASV,fill = 0,
        value.var = "diff") %>%
  scale
display2=data.frame(display)

dend_nutr <- as.dendrogram(hclust(dist(t(display2))))
dend_nutr_data <- dendro_data(dend_nutr)
nutr_order = dend_nutr_data$labels[,3]

# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data <- with(
  segment(dend_nutr_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

# Use the dendrogram label data to position the gene labels
gene_pos_table <- with(
  dend_nutr_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

# Table to position the samples
sample_pos_table <- data.frame(sample = nutr_order) %>%
  mutate(x_center = (1:n()), width = 1)

# Neglecting the gap parameters
heatmap_data <- display2 %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  cross_join(gene_pos_table) %>%
  cross_join(sample_pos_table)

# Limits for the vertical axes
gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 1 * height))) + 0.1 * c(-1, 1) # extra spacing: 0.1

species_order = dend_nutr_data$labels[,3]
dend_nutr_data$labels
Good_Bad$ASV=factor(Good_Bad$ASV, nutr_order)
Good_Bad$Species=factor(Good_Bad$Species, c('BW','NRO')) #species_order
nrow(list_log2foldchange2)

####Heatmap
plt_hmap = ggplot(data = Good_Bad, aes(Species,ASV)) +
  geom_raster(aes(fill = diff))+
  theme_few() +
  #geom_text(aes(label = Phylum),color = "black",size = 5)+
  scale_fill_paletteer_c("pals::kovesi.diverging_bwr_40_95_c42", #kovesi.diverging_bwr_55_98_c37
                         limits = c(-0.16,0.05),na.value = "#D9D9D9")+
  geom_tile(aes(color = Significance),fill = '#00000000', linewidth = 0.05,width = 0.9,height = 0.95) + #
  scale_color_manual(values = c('transparent',"black"),na.value =  "transparent",name = "Significance vs Water") + #Significance Genotype vs Col-0
  theme(axis.text.x = element_text(angle = -45, hjust=-0.05),#axis.text.y=element_blank(),
        axis.title = element_blank())


# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_reverse() + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_tufte() + 
  theme(panel.grid.minor = element_blank(),axis.text.y=element_blank(),
        axis.title = element_blank()) #

library(cowplot)
plot_grid(plt_dendr, plt_hmap, align = 'h', rel_widths = c(0.2, 1))












##############
library(reshape2)
overall_p <- function(my_model) {
  f <- summary(my_model)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes (p) <- NULL
  return(p)
}

tableSyncom = cbind(ps.phyla.perc@sam_data[,c(1,3,6,7)], ps.phyla.perc@otu_table)
tableSyncom2=tableSyncom

melted_tableSyncom <- tableSyncom2 %>% melt
growth = data.frame(ps.phyla.perc@sam_data[,c(1,9,10)])
melted_tableSyncom2 = merge(melted_tableSyncom,growth,by='ID' )

list_corr_Dbh <- c()

for(id in melted_tableSyncom2$variable  %>% unique){
  melted_sub <- melted_tableSyncom2 %>% subset(variable == id) %>% droplevels
  
  MetaCorr = lm(Dbh~value, data=melted_sub)
  summary(MetaCorr)$adj.r.squared
  
  ASV = data.frame(cbind(id,summary(MetaCorr)$adj.r.squared))
  Corr = cbind(ASV,overall_p(MetaCorr))
  
  list_corr_Dbh <- rbind(list_corr_Dbh, Corr)
}

list_corr_Dbh2=data.frame(list_corr_Dbh)
colnames(list_corr_Dbh2)=c('Metabolite','R-ajust' ,'p.value')
head(list_corr_Dbh2)

####Adjust the p-value
list_corr_Dbh2$SignPvalue <- "NoSignificant"
pval_thres <- 0.05
list_corr_Dbh2$SignPvalue[which(list_corr_Dbh2$p.value < pval_thres)] <- "p < 0.05"
list_corr_Dbh2$SignPvalue <- list_corr_Dbh2$SignPvalue %>% factor

###Subsample
Corr_species_Dbh <- list_corr_Dbh2[list_corr_Dbh2$SignPval ==  'p < 0.05',]
nrow(Corr_species_Dbh)
#Bacillus, Gp3, GP1, Chitinophagaceae_unclassified, Aciditerrimonas, Arthrobacter

Bacillus = ggplot(melted_tableSyncom2[melted_tableSyncom2$variable=='Bacillus',], #
                   aes(x=value , y=Dbh))+
  theme_bw() +
  #facet_grid(.~Species*State, scales = "free")+
  geom_point(size=3, aes(pch =Species, col=State)) + 
  geom_smooth(method=lm, se=TRUE)+ #aes(col=State),
  theme(legend.position='none')+
  labs(y='Dhb (cm)', x='Abundance', title='Bacillus')+
  scale_color_manual(values = paleta_alive)


Gp3 = ggplot(melted_tableSyncom2[melted_tableSyncom2$variable=='Gp3',], #
                   aes(x=value , y=Dbh))+
  theme_bw() +
  #facet_grid(.~Species*State, scales = "free")+
  geom_point(size=3, aes(pch =Species, col=State)) + 
  geom_smooth(method=lm, se=TRUE)+#aes(col=State),
  theme(legend.position='none')+
  labs(y='Dhb (cm)', x='Abundance', title='Gp3')+
  scale_color_manual(values = paleta_alive)


Gp1 = ggplot(melted_tableSyncom2[melted_tableSyncom2$variable=='Gp1',], #
             aes(x=value , y=Dbh))+
  theme_bw() +
  #facet_grid(.~Species*State, scales = "free")+
  geom_point(size=3, aes(pch =Species, col=State)) + 
  geom_smooth(method=lm, se=TRUE)+ #aes(col=State),
  theme(legend.position='none')+
  labs(y='Dhb (cm)', x='Abundance', title='Gp1')+
  scale_color_manual(values = paleta_alive)


Aciditerrimonas = ggplot(melted_tableSyncom2[melted_tableSyncom2$variable=='Aciditerrimonas',], #
                         aes(x=value , y=Dbh))+
  theme_bw() +
  #facet_grid(.~Species*State, scales = "free")+
  geom_point(size=3, aes(pch =Species, col=State)) + 
  geom_smooth(method=lm, se=TRUE)+#aes(col=State),
  theme(legend.position='none')+
  labs(y='Dhb (cm)', x='Abundance', title='Aciditerrimonas')+
  scale_color_manual(values = paleta_alive)

Arthrobacter = ggplot(melted_tableSyncom2[melted_tableSyncom2$variable=='Gp3',], #
                         aes(x=value , y=Dbh))+
  theme_bw() +
  #facet_grid(.~Species, scales = "free")+
  geom_point(size=3, aes(pch =Species, col=State)) + 
  geom_smooth(method=lm, se=TRUE)+#aes(col=State),
  theme(legend.position='none')+
  labs(y='Dhb (cm)', x='Abundance', title='Arthrobacter')+
  scale_color_manual(values = paleta_alive)

Chitinophagaceae_unclassified = ggplot(melted_tableSyncom2[melted_tableSyncom2$variable=='Chitinophagaceae_unclassified',], #
                                     aes(x=value , y=Dbh))+
  theme_bw() +
  #facet_grid(.~Species, scales = "free")+
  geom_point(size=3, aes(pch =Species, col=State)) + 
  geom_smooth(method=lm, se=TRUE)+ #aes(col=State),
  theme(legend.position='none')+
  labs(y='Dhb (cm)', x='Abundance', title='Chitinophagaceae unclassified')+
  scale_color_manual(values = paleta_alive)

 
ggarrange(Bacillus, Gp3, Gp1, 
          Aciditerrimonas, Arthrobacter, Chitinophagaceae_unclassified,
          labels = c("A", "B","C", "D","E", "F"),  common.legend = TRUE)

############Randon forest 2
#
library(lavaan)
library(semPlot)
library(lavaanPlot)
library(randomForest)
library(car)
library(stargazer)
library(ggfortify)
library(ggcorrplot)
library(ggbiplot)
library(devtools)
library(rfPermute)


#Datos
Tabla <- read.table("RF version 2.txt", header = TRUE)

TableBW = Table[Table$Species == 'BW',]
# Split into Train and Validation sets
# Training Set : Validation Set = 70 : 30 (random)
set.seed(100)
train <- sample(nrow(TableBW), 0.7*nrow(TableBW), replace = FALSE)
TrainSet1 <- TableBW[train,]
ValidSet1 <- TableBW[-train,]
summary(TrainSet1)
summary(ValidSet1)

# Create a Random Forest model with default parameters
model <- randomForest(y = TableBW$Dbh, x = TableBW[,c('CEC','Ca','SumNH4','Mg','Base_Sat','pH_CaCl2','pH_H2O','Corg','Clay','K','PBray','Na','Silt','Sand')], #, 'Gp1','Arthrobacter','Aciditerrimonas','Bacillus','Gp3','Chitinophagaceae_unclassified' 
                      data = TrainSet1, importance = TRUE, ntree = 500, mtry = 8, proximity = TRUE)
model
plot(model)

# Predicting on train set
predTrain1 <- predict(model, TrainSet1, type = "class")
# Checking classification accuracy
table(predTrain1, TrainSet1$Dbh) 

# Predicting on Validation set
predValid1 <- predict(model, ValidSet1, type = "class")
# Checking classification accuracy
mean(predValid1 == ValidSet1$Dbh)
table(predValid1,ValidSet1$Dbh)

# To check important variables
BW_importance  = importance(model) #         
varImpPlot(model)

BW_importance = data.frame(cbind(rownames(BW_importance),BW_importance)) 
colnames(BW_importance) = c('V1','BW_IncMSE', 'BW_IncNodePurity')
BW_importance$BW_IncMSE = as.numeric(BW_importance$BW_IncMSE)
BW_importance$V1 = factor(BW_importance$V1, rev(c("Sand","Silt","CEC","Ca","pH_CaCl2","K","Corg","pH_H2O","Clay","Base_Sat","Mg","PBray","Na","SumNH4")))

ggplot(data=BW_importance, aes(x=V1, y=BW_IncMSE, fill=BW_IncMSE)) +
  geom_bar(stat="identity")+
  theme_few()+ labs(y= '%IncMSE',x= 'Soil Properties')+
  coord_flip()+ theme(legend.position = 'none')


