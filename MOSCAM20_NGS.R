necessary_packages <- c('tidyverse', 'data.table', 'scales', 'knitr', 'grid', 'magrittr', 'DT', 'here', 'gridtext',
                        'ggpubr', 'ggthemes', 'reshape2', 'viridis', 'pals', 'circlize',
                        'vegan', 'phyloseq', 'metagenomeSeq', 'ComplexHeatmap', 'decontam', 'ggplot2','DT')
#install.packages(necessary_packages) #Install necessary and do this just once 
lapply(necessary_packages, library, character.only = TRUE) #Load 

#' heatmap function
plot_abundance <- function (obj, n, norm = TRUE, log = TRUE, fun = sd, rowclust = "bray", colclust = "bray", ...) 
{
  mat <- returnAppropriateObj(obj, norm, log)
  otusToKeep <- which(rowSums(mat) > 0)
  otuStats <- apply(mat[otusToKeep, ], 1, fun)
  otuIndices <- otusToKeep[order(otuStats, decreasing = TRUE)[1:n]]
  mat2 <- mat[otuIndices, ]
  mat3 <- mat2[order(row.names(mat2)),]
  mat4 <- as.data.frame(t(mat3))
  Heatmap(mat3, clustering_distance_columns = vegdist(mat4, method = colclust), clustering_distance_rows = vegdist(mat3, method = rowclust), ...)
}


#' ***
#' # Eukaryotic virome analysis
#' ## Prepare the data
#' ### Load the OTU table, taxonomy file and metadata into R
OTU <- read.table ("moscam2020-abundance.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
#' Information on OTU table 
ncol(OTU) #Number of samples 
nrow(OTU) #Number of contigs 
summary(rowSums(OTU))
summary(colSums(OTU))
sum(rowSums(OTU))   ## Total number of reads (abundance)
sum(colSums(OTU))   ## Total number of reads (abundance)

tax <- read.table("MOS20tax_trimmed", header=T, row.names=1, sep="\t", dec=".")
#' Information on tax table 
ncol(tax) 
nrow(tax)

meta <- read.table("Metadata20.csv", header=TRUE, row.names = 1, sep=";", dec=".")

#Shorten sample names 
colnames(OTU) = gsub("[EBYaeud]*20_KC_","",colnames(OTU),fixed = F)
names(OTU) <- gsub(x = names(OTU), pattern = "\\.", replacement = "-")
rownames(meta) = gsub("[EBYaeud]*20_KC_","",rownames(meta),fixed = F)
meta <- cbind(Sample=rownames(meta), meta)
#+ echo=FALSE
datatable(meta)


#' ### Make a phyloseq object
OTU.UF <- otu_table(as.matrix(OTU), taxa_are_rows=T)
tax.UF <- tax_table(as.matrix(tax))
meta.UF <- sample_data(meta)

MOS20 <- phyloseq(OTU.UF, tax.UF, meta.UF)
MOS20

#' #### Subset only viruses and remove prokaryotic viruses and EVEs
EVE_phage <- c("Atrato Retro-like virus", "Gurupi chuvirus-like 1", "Aedes aegypti To virus 1",
               "Aedes aegypti To virus 2", "Guato virus", "Kaiowa virus", "Atrato Chu-like virus 1",
               "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp.", 
               "Aedes aegypti virga-like virus", "Verdadero virus", "Aedes rhabdo-like virus",
               "Kvarnon virus", "Aedes binegev-like virus 2", "unclassified 9", "unclassified 10")

MOS20.V1 <- subset_taxa(MOS20, !is.element(Species, EVE_phage))
MOS20.V1 <- subset_taxa(MOS20.V1, Order!="Caudovirales")
MOS20.V1 <- subset_taxa(MOS20.V1, Phylum!="Phage")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Verdadero virus ")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Aedes rhabdo-like virus ")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Shilevirus leptomonadis")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Usinis virus")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Drosophila-associated totivirus 3")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Grizzly bear anellovirus 10")
# # #MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Anelloviridae sp.")
# MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Hubei mosquito virus 4")
MOS20.V1 <- subset_taxa(MOS20.V1, Species!="Formosus virus")


#' #### Subset only viruses and remove prokaryotic viruses and EVEs
MOS20.V2 <- subset_taxa(MOS20.V1, Kingdom=="Viruses")
MOS20.V2

#' **Info on sample variables and level of taxonomic ranks:**
sample_variables(MOS20.V2)
rank_names(MOS20.V2)

#' #### Agglomerate taxa on Viral_species level
MOS20_species <- MOS20.V2 %>%
  tax_glom(taxrank = "Species")
MOS20_species


#' **Only keep taxa with more than 0 reads:**
MOS20_species <- prune_taxa(taxa_sums(MOS20_species) > 0, MOS20_species)
MOS20_species

#' **Only keep samples with more than 0 viral reads:**
MOS20_final <- prune_samples(sample_sums(MOS20_species) > 0, MOS20_species)
MOS20_final

#Sum info
rowSums(otu_table(MOS20_final))
colSums(otu_table(MOS20_final))
sum(otu_table(MOS20))

#' **Convert phyloseq object to metagenomeseq object to make heatmap:**
MOS20_metaseq <- phyloseq_to_metagenomeSeq(MOS20_final)

#' **Aggregate by species:**
MOS20_metaseq_species <- aggregateByTaxonomy(MOS20_metaseq, lvl = "Species", norm = F, 
                                             aggfun = colSums, out = "MRexperiment", alternate = T)

#' **Count number of unique Viral_species:**
n_species <- length(unique(featureData(MOS20_metaseq_species)$Species))
n_species


#Count how many mosquito species there are and assign colors
unique(sample_data(MOS20_final)$Species)
length(unique(meta$Species))
#myColors <- c(viridis::plasma(7,0.8, begin=0, end = 1, direction = 1))
myColors <- c("#338068", "#814C8A")
#myColors <- c("#814C8A", "#EBC71B", "#E181AC", "#009E73", "#824CDE", "#EF3924", "#338068",  "#4A67B1", "#5D1F0B" )
names(myColors) <- levels(as.factor(meta$Species))
myColors
pal.bands(myColors, main="Mosquito species colors")


#Assign colors to location
length(unique(meta$Location))
locColors <- c(viridis::viridis(8,1, begin=0, end = 1, direction = 1))
locColors <- c("#F68B68", "#78ADBE", "#CDABD1","#5CBA47")
names(locColors) <- levels(as.factor(meta$Location))
locColors
pal.bands(locColors, main="location colors")

#' **Set heatmap colors:**
heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)

#' **Assign mosquito species for each sample:**
mosquito_species <- pData(MOS20_metaseq_species)$Species

#' **Assign colors to location:**
location <- pData(MOS20_metaseq_species)$Location


#' **Calculate average BLASTx values:**
blastx <- read.table("mean_blastx.tsv", header=T, row.names=1, dec=".", sep="\t")
blastx.UF <- otu_table(as.matrix(blastx), taxa_are_rows=T)
blastx.ps <- phyloseq(blastx.UF, tax_table(prune_taxa(taxa_sums(MOS20_final) >= 1, MOS20_final)))
blastx_metaseq <- phyloseq_to_metagenomeSeq(blastx.ps)
blastx_metaseq

#Aggregate by species
blastx_mean = aggregateByTaxonomy(blastx_metaseq, lvl = "Species", norm = F, aggfun = mean, out = "MRexperiment", alternate = T)
blastx_mean

#' **Store average BLASTx identities in dataframe:**
rowanno <- as.data.frame(returnAppropriateObj(blastx_mean, log=F, norm=F))
colnames(rowanno)[1] <- "blastx"

#' **Create color function for BLASTx values:**
col_fun=colorRamp2(c(0,100), c("white","deepskyblue4"))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(col_fun(0:100), main="Blastx % identity colors")

#' ## Eukaryotic virome heatmap
left_ra=rowAnnotation("Blastx"= rowanno$blastx,
                      col=list("Blastx"=col_fun),
                      show_annotation_name=T,
                      annotation_name_side = "top",
                      annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
                      annotation_name_rot = 270,
                      annotation_legend_param=list("Blastx" = list(title ="Blastx % Identity", 
                                                                   direction="horizontal",
                                                                   at=c(0,25,50,75,100))))



column_ha = HeatmapAnnotation(Location=location,
                              'Mosquito species' = mosquito_species,
                              show_annotation_name = T,
                              annotation_label = gt_render(c("Location", "Species")),
                              annotation_name_side = "right",
                              annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
                              annotation_name_rot = 0,
                              annotation_legend_param=list('Mosquito species'= list(title ='Mosquito species', 
                                                                                    labels_gp = gpar(fontface="italic")),
                                                           'Location' = list(title='Location')),
                              col = list('Mosquito species'=c(myColors), Location=c(locColors)))
#n depends on the amount of taxa in 'MOS20_metaseq_species'
hm <- plot_abundance(MOS20_metaseq_species, n = n_species, log = T, norm = F, colclust = "bray",
                     col = heatmapCols,
                     name = "Log2 Read Counts", 
                     top_annotation=column_ha,
                     left_annotation=left_ra,
                     row_names_gp = gpar(fontsize = 7),
                     show_column_names = T,
                     column_names_rot = 45,
                     column_names_gp=gpar(fontface=1, fontsize=7),
                     heatmap_legend_param = list(direction = "horizontal"),
                     cluster_rows = T,
                     row_title_gp = gpar(fontsize=7),
                     row_title_rot=0,
                     border=F)
#+ echo=TRUE, fig.width=15, fig.height=9
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", 
     merge_legend = T, legend_grouping="original")



#'*Alpha diversity*
#'*Alpha diversity*
#'*Alpha diversity*
#Wilcox
#Based on mosquito species 
#Based on mosquito species 

mosColors <- c("#338068", "#814C8A")
names(mosColors) <- levels(as.factor(meta$Species))
mosColors

alpha_mos <- plot_richness(MOS20_final, measures=c("Observed", "Simpson"), x="Species", color = "Species")+
  geom_boxplot()+
  geom_jitter(width=0.05, size=3)+
  stat_compare_means(aes(label= ..p.signif..), show.legend = F, size = 10) +  # Set size of significance asterisks
  theme_bw() +
  scale_fill_manual(values = mosColors) +  # Assigning colors manually
  scale_color_manual(values = mosColors) +  # Assigning colors manually
  scale_color_manual(values = myColors, name = "Mosquito Species") +
  theme(axis.text.y=element_text(angle=90,hjust=0.5, size = 22))+
  theme(axis.title.y = element_text(size = 22, face = "bold")) + # Increase size of y-axis title
  theme(strip.text.x = element_text(size = 22),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 22, face = "italic"),
        legend.title = element_text(face = "bold",size = 22),  # Change legend title
        plot.title = element_text(size = 30, face = "bold")) +  # Increase title size and bold
  labs(title = "a") 

alpha_mos 


#Wilcox
#Based on location species 
#Based on location species 

locColors <- c("#F68B68", "#78ADBE", "#CDABD1","#5CBA47")
names(locColors) <- levels(as.factor(meta$Location))
locColors

alpha_loc <- plot_richness(MOS20_final, measures = c("Observed", "Simpson"), x = "Location", color = "Location") +
  geom_boxplot() +
  geom_jitter(width = 0.05, size = 3) +
  theme_bw() +
  scale_fill_manual(values = locColors) +
  scale_color_manual(values = locColors) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size = 22)) +
  theme(axis.title.y = element_text(size = 22, face = "bold")) +
  theme(strip.text.x = element_text(size = 22),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22, face = "bold"),
        plot.title = element_text(size = 30, face = "bold")) +
  labs(title = element_text("b", face = "bold")) +
  scale_y_continuous(limits = c(0, NA)) +
  stat_compare_means(method = "wilcox.test", aes_string(x = "Climatic_zone", y = "value", group = "Climatic_zone"),
                     label = "p.signif", hide.ns = TRUE, geom = "text", size = 10, show.legend = FALSE, tip.length = 0.0,
                     comparisons = list(c(1, 2), c(1, 3)))

alpha_loc



library(gridExtra)
library(cowplot)
Alpha <- grid.arrange(alpha_mos, alpha_loc, nrow = 1)
Alpha

## 2. PCoA
PCA <- ordinate(MOS20_final, "PCoA")
p <- plot_ordination(MOS20_final, PCA, type = "samples", color = "Species", shape = "Location") +
 geom_point(size = 3, aes(group = "Species")) +
  geom_text(aes(label = Sample), vjust = 1.5) + 
  theme_bw()+
  scale_fill_manual(values = myColors) +  # Assigning colors manually
  scale_color_manual(values = myColors, name = "Mosquito Species", 
                     labels = expression(italic("Ae. africanus"), italic("Ae. albopictus"))) +  # Assigning custom italicized labels for legend
  theme(legend.text.align = 0,
        legend.key.size = unit(1.5, "lines"),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(size = 14),
        axis.text = element_text(size = 14),   # Increase size of axis text
        axis.title = element_text(size = 16, face = "bold")) + # Increase size of axis title
  ggtitle("") #Add title to plot b

p <- p +   stat_ellipse(type = "norm", linetype = 3, aes_string(group = "Species"), show.legend = FALSE)  # Adds circles 
p


# Load taxonomy info from clean phyloseq object (MOS20.V2)
MOS20_NGS_mastertable <- tax_table(MOS20.V2)
MOS20_NGS_mastertable <- as.data.frame(MOS20_NGS_mastertable)

# Extract row names, "Family," and "Species" columns into a new table
Master_table <- data.frame(
  Contigs = rownames(MOS20_NGS_mastertable),
  Family = MOS20_NGS_mastertable$Family,
  Species = MOS20_NGS_mastertable$Species
)


# Extract the sample names from the Contig column
Master_table$Sample_name <- sub(".*_P([^_]+)$", "P\\1", Master_table$Contigs)
#sample_names <- sapply(strsplit(Master_table$Contig, "_"), function(x) paste(x[7:length(x)], collapse = "_"))


# Extract all the length from the contig names
Contig_lengths <- str_match(Master_table$Contigs, "length_(\\d+)_")[, 2]
Master_table$Contig_length <- as.numeric(Contig_lengths)

# Calculate the total abundance per sample
Total_abundance <- rowSums(otu_table(MOS20.V2))
# Add the total abundance as a new column in the Master_table
Master_table$Total_abundance <- Total_abundance

#Extract coverage from contig name
library(stringr)
cov <- str_match(Master_table$Contigs, "NODE_\\d+_length_\\d+_cov_([^.]+)")[, 2]
Master_table$Coverage <- cov


# Assuming that the column in the `sample_data` containing mosquito species is named "Species"
sample_data <- sample_data(MOS20.V2)
Mosquito_species <- sample_data$Species
Mosquito_species <- sample_data(MOS20.V2)$Species
# Match the sample names and merge the species information into the Master_table
Master_table$Mosquito_species <- Mosquito_species[match(Master_table$Sample_name, rownames(sample_data))]


# Match the sample names and merge the Location information into the Master_table
Location <- sample_data(MOS20.V2)$Location
Master_table$Location <- Location[match(Master_table$Sample_name, rownames(sample_data))]

# Calculate the total abundance in percentage
Relative_abundance <- (Master_table$Total_abundance / sum(Master_table$Total_abundance)) * 100
# Add the total abundance percentage as a new column in the Master_table
Master_table$Relative_abundance <- Relative_abundance


#'*MOVE TO RELATIVE ABUNDANCE SCRIPT*
#'*MOVE TO RELATIVE ABUNDANCE SCRIPT*


#'*1. Treemap*
#'*1. Treemap*
#'
# Aggregate the Relative_Abundance values by Family
Aggregated_Family <- Master_table %>%
  group_by(Family) %>%
  summarise(Relative_abundance = sum(Relative_abundance))

# Sort the aggregated data by Relative_Abundance in descending order
Aggregated_Family <- Aggregated_Family %>%
  arrange(desc(Relative_abundance))
Aggregated_Family


#Assign colors to  
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE) #Color blind friendly 
display.brewer.pal(n = 11, name = 'Paired') # Change 'n' and 'name' depending on the color palette 
# Hexadecimal color specification 
brewer.pal(n = 11, name = "Paired")

# count the number of unique families
n_families <- Master_table %>%
  distinct(Family) %>%
  nrow()
n_families

n <- 17 #Number of families 

mypalette<- c("#C43B80", "#1B9E77", "#7570B3","#E31A1C", "#824CDE", "#A6761D", "#666666", "#82906F",
              "#FFFF99","#A6CEE3","#E05fAD","#74EA56","#DBA875","#E6DC50","#27408B", "#FF5733", "#5FBBEA")

names(mypalette) <- levels(as.factor(Master_table$Family))
mypalette
pal.bands(mypalette, main="Virus Family colors")


#Create a treemap 
### Make treemap (TM)
library(treemapify)
library(randomcoloR)
library(ggplot2)
library(treemap)
library(forcats)

TM <- ggplot(data = Aggregated_Family, aes(area = Relative_abundance, fill = Family, label = paste0(Family, "\n"))) +
  geom_treemap(color = "black") +
  geom_treemap_text(fontface = "italic", colour = "black", place = "centre", size = 14) +
  labs(title = element_text("a", face = "bold")) +
  theme(legend.position = "right")+
  scale_fill_manual(values = mypalette) +  # Custom color palette
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "Virus families",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    title.theme = element_text(face = "bold", size = 14),
    label.theme = element_text(face = "italic", size = 12)))

TM


### Make treenap for only unclassified Family 
###
# Filter data to include only "Unclassified" rows
Unclassified_fam <- Master_table %>%
  filter(Family == "unclassified")
# Recalculate relative abundance of filtered data (unclassified) in percentage and add to new column
Relative_abundance_unclass <- (Unclassified_fam$Relative_abundance / sum(Unclassified_fam$Relative_abundance)) * 100
# Add the total abundance percentage as a new column in the Unclassified_fam data frame
Unclassified_fam$Relative_abundance_unclass <- Relative_abundance_unclass
# Aggregate Viral_species values
Aggregated_Unclassified_fam <- Unclassified_fam %>%
  group_by(Species) %>%
  summarise(Relative_abundance_unclass = sum(Relative_abundance_unclass))

# # Reverse the order of the 'Viral_species' factor levels
# Aggregated_Unclassified_fam$Viral_species <- fct_reorder(Aggregated_Unclassified_fam$Viral_species, Aggregated_Unclassified_fam$Relative_abundance_unclass)
# Aggregated_Unclassified_fam$Species <- fct_reorder(Aggregated_Unclassified_fam$Species, -Aggregated_Unclassified_fam$Relative_abundance_unclass)

# Sort the aggregated data by Relative_Abundance in descending order
Aggregated_Unclassified_fam<- Aggregated_Unclassified_fam %>%
  arrange(desc(Relative_abundance_unclass))
Aggregated_Unclassified_fam


n_families_unclass <- Aggregated_Unclassified_fam %>%
  distinct(Species) %>%
  nrow()
n_families_unclass
n_fam <- 11
palette <- distinctColorPalette(n_fam)

# To choose your own colors 
palette <- c("#D86ED9", "#E04E8E", "#75E3A0", "#88954D", "#DDACA5", "#715CD1", "#D96255","#D4E691", "#AD39E5", "#80E7E1", "#976790", "#E6A764", "#DAA7DF", "#D0E1BD","#69AE9D", "#828ADD", "#DBE450", "#D3D9E2", "#7CB2D9", "#76E760")
names(palette) <- levels(as.factor(Master_table$Family$unclassified))
palette
pal.bands(mypalette, main="Virus Family colors")

TMU <- ggplot(data = Aggregated_Unclassified_fam, aes(area = Relative_abundance_unclass, fill = Species, label = paste0(Species, "\n"))) +
  geom_treemap(colour = "black") +
  geom_treemap_text(colour = "black", place = "centre") +
  labs(title = "b") +
  theme(legend.position = "right")+
  scale_fill_manual(values = palette) +  # Custom color palette
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "unclassified viruses",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    title.theme = element_text(size = 14, face = "bold"),
    label.theme = element_text(size = 12)))

TMU

library(gridExtra)
library(cowplot)
# Assuming TM and TMU are already defined ggplot objects
Treemap <- grid.arrange(TM, TMU, ncol = 1)
Treemap

combined_Treemap <- ggdraw() +
  draw_plot(Treemap, width = 1, height = 1, x = 0, y = 0) +
  draw_label("Combined Plot Title", x = 0.5, y = 1, 
             size = 16, vjust = 1, hjust = 0.5, fontface = "bold")
combined_Treemap


#'*2. Piechart*
#'*2. Piechart*
#'
# Calculate the total relative abundance for each mosquito species
Reads_mos <- aggregate(Relative_abundance ~ Mosquito_species, data = Master_table, sum)
Reads_mos

Pie_Mos <- ggplot(Reads_mos, aes(x = "", y = Relative_abundance, fill = Mosquito_species, label = Mosquito_species)) +
  geom_bar(colour = "black", stat = "identity", width = 1) +
  geom_text(aes(label = Mosquito_species), position = position_stack(vjust = 0.7), fontface = "italic") +  # Add labels
  coord_polar("y", start = 0) + 
  labs(title = "a") +
  theme_void() +
  scale_fill_manual(values = mosColors) +  # Assigning colors manually
  scale_color_manual(values = mosColors) + # Assigning colors manually
  #scale_fill_manual(values = palette)
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "Virus families",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    title.theme = element_text(face = "bold", size = 14),
    label.theme = element_text(face = "italic", size = 12)))
Pie_Mos

Pie_Mos <- ggplot(Reads_mos, aes(x = "", y = Relative_abundance, fill = Mosquito_species, label = paste0(Mosquito_species, "\n (", round(Relative_abundance, 2), "%)"))) +
  geom_bar(colour = "black", stat = "identity", width = 1) +
  geom_text(size = 4, position = position_stack(vjust = 0.65), fontface = "italic") +  # Add labels and adjust position
  coord_polar("y", start = 0) + 
  labs(title = "a") +
  theme_void() +
  scale_fill_manual(values = mosColors) +  # Assigning colors manually
  scale_color_manual(values = mosColors) + # Assigning colors manually
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "Virus families",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    title.theme = element_text(face = "bold", size = 14),
    label.theme = element_text(face = "italic", size = 12)))
Pie_Mos



#'*Details on Pie*
#'*Details on Pie*
# Filter data for Ae. africanus
# Filter data for Ae. africanus
africanus_data <- subset(Master_table, Mosquito_species == "Ae. africanus")
# Get unique families and their abundances
unique_families_afri <- unique(africanus_data$Family)
abundances_afri <- sapply(unique_families_afri, function(x) sum(africanus_data$Relative_abundance[africanus_data$Family == x]))
abundances_afri
#Number of virus spp. in Ae. africanus 
num_species_africanus <- n_distinct(filter(Master_table, Mosquito_species == "Ae. africanus")$Species)
num_species_africanus 

# Filter data for Ae. albopictus
# Filter data for Ae. albopictus
albopictus_data <- subset(Master_table, Mosquito_species == "Ae. albopictus ")
# Get unique families and their abundances
unique_families_albo <- unique(albopictus_data$Family)
abundances_albo <- sapply(unique_families_albo, function(x) sum(albopictus_data$Relative_abundance[albopictus_data$Family == x]))
abundances_albo
#Number of virus spp. in Ae. albopictus
num_species_albopictus <- n_distinct(filter(Master_table, Mosquito_species == "Ae. albopictus ")$Species)
num_species_albopictus


#Number of virus spp. in Ae. albopictus and Ae. africanus
# Extract the unique species for Ae. albopictus and Ae. africanus
species_albopictus <- unique(filter(Master_table, Mosquito_species == "Ae. albopictus ")$Species)
species_africanus <- unique(filter(Master_table, Mosquito_species == "Ae. africanus")$Species)
# Find the common species between the two
Aedes <- intersect(species_albopictus, species_africanus)
Aedes
# Count the number of species found in both
num_species_both_Aedes <- length(Aedes)
num_species_both_Aedes



# Calculate the total relative abundance for each mosquito species
Reads_Loc <- aggregate(Relative_abundance ~ Location, data = Master_table, sum)
Reads_Loc

Pie_Loc <- ggplot(Reads_Loc, aes(x = "", y = Relative_abundance, fill = Location, label = paste0(Location, "\n (", round(Relative_abundance, 0), "%)"))) +
  geom_bar(colour = "black", stat = "identity", width = 1) +
  geom_text(size = 2.8, colour = "black",position = position_stack(vjust = 0.55)) +  # Add labels
  coord_polar("y", start = 0) + 
  labs(title = "b") +
  theme_void() +
  scale_fill_manual(values = locColors) +  # Assigning colors manually
  scale_color_manual(values = locColors) + # Assigning colors manually
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "unclassified viruses",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    title.theme = element_text(size = 14, face = "bold"),
    label.theme = element_text(size = 12)))
Pie_Loc

# Assuming TM and TMU are already defined ggplot objects
Pie<- grid.arrange(Pie_Mos, Pie_Loc, nrow = 1)
Pie

combined_Pie <- ggdraw() +
  draw_plot(Pie, width = 1, height = 1, x = 0, y = 0) +
  draw_label("Combined Plot Title", x = 0.5, y = 1, 
             size = 16, vjust = 1, hjust = 0.5, fontface = "bold")
combined_Pie


#'*Details on Pie*
#'*Details on Pie*
# Filter data for Ae. africanus
# Filter data for Ae. africanus
buea_data <- subset(Master_table, Location == "Buea")
# Get unique families and their abundances
unique_families_buea <- unique(africanus_data$Family)
abundances_buea  <- sapply(unique_families_buea , function(x) sum(buea_data$Relative_abundance[buea_data$Family == x]))
abundances_buea 

# Filter data for Ae. albopictus
# Filter data for Ae. albopictus
yaounde_data <- subset(Master_table, Location == "Yaounde ")
# Get unique families and their abundances
unique_families_yaounde <- unique(yaounde_data$Family)
abundances_yaounde <- sapply(unique_families_yaounde, function(x) sum(yaounde_data$Relative_abundance[yaounde_data$Family == x]))
abundances_yaounde

