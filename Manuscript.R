#Install and load packages
if (!require("BiocManager", quietly = TRUE))  #Do this just the first time 
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")


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

#setwd("~/Documents/MDKC Files/MOSCAM20/Manuscript MOS20/Manuscript2023/NGS/R/Inputs")
here::i_am("Manuscript.R")

#' ***
#' # Eukaryotic virome analysis
#' ## Prepare the data
#' ### Load the OTU table, taxonomy file and metadata into R
OTU <- read.table ("Data/moscam2020-abundance.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
#' Information on OTU table 
ncol(OTU) #Number of samples 
nrow(OTU) #Number of contigs 
summary(rowSums(OTU))
summary(colSums(OTU))
sum(rowSums(OTU))   ## Total number of reads (abundance)
sum(colSums(OTU))   ## Total number of reads (abundance)

tax <- read.table("Data/MOS20tax_final.tsv", header=T, row.names=1, sep="\t", dec=".")
#' Information on tax table 
ncol(tax) 
nrow(tax)

#tax <- read.table("moscam2020-oldtaxfile.tsv", header=TRUE, row.names=1, sep="\t", dec=".")
meta <- read.table("Data/Metadata.csv", header=TRUE, row.names = 1, sep=";", dec=".")

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
               "Chuvirus Mos8Chu0", "Chibugado virus", "Prokaryotic dsDNA virus sp.", "Aedes aegypti virga-like virus")

MOS20.V2 <- subset_taxa(MOS20, !is.element(Species, EVE_phage))
MOS20.V2 <- subset_taxa(MOS20.V2, Order!="Caudovirales")
MOS20.V2 <- subset_taxa(MOS20.V2, Phylum!="Phage")


#MOS20_wo_unclass <- subset_taxa(MOS20.V2, Family!="unclassified")

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
length(unique(meta$Species))
myColors <- c(viridis::plasma(2,0.8, begin=0, end = 0.9, direction = 1))
#myColors <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00")
names(myColors) <- levels(as.factor(meta$Species))
myColors

#Assign colors to location
length(unique(meta$Location))
locColors <- c(viridis::viridis(4,1, begin=0.4, end = 1, direction = 1))
#locColors <- brewer.pal(8, 'Set1')
#display.brewer.all(type = 'qual')
names(locColors) <- levels(as.factor(meta$Location))
locColors

#' **Set heatmap colors:**
heatmapCols <- colorRampPalette(brewer.pal(9, "YlOrRd"))(20)

#' **Assign mosquito species for each sample:**
mosquito_species <- pData(MOS20_metaseq_species)$Species

#' **Assign colors to location:**
location <- pData(MOS20_metaseq_species)$Location

#' **Calculate average BLASTx values:**
blastx <- read.table("blastx-pid.tsv", header=F, row.names=1, dec=".", sep="\t")
colnames(blastx)<-"Average_id"
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

# Get the taxonomy table from the phyloseq object
taxonomy_table <- tax_table(MOS20_final)
taxonomy_table <- as.data.frame(taxonomy_table)
rownames(taxonomy_table)

# Extract row names and the "Family" column into a new table
family.df <- data.frame(Contig = rownames(taxonomy_table), Family = taxonomy_table$Family)


#' **Create colorpalette for viral families:**
MOS20_smelt <- psmelt(MOS20_final)
FamLevel <- levels(as.factor(MOS20_smelt[(MOS20_smelt$Family!="unclassified"),]$Family))
FamLevel <- c(FamLevel, "unclassified")
MOS20_smelt$Family <- factor(MOS20_smelt$Family, levels = FamLevel)

myFamCol <- c(stepped(n=17))
names(myFamCol) <- levels(as.factor(MOS20_smelt$Family))
#+ echo=FALSE, fig.width=7, fig.height=4
pal.bands(myFamCol, main = "Viral families")

#OR.  OR   OR 
# Assuming 'family.df' is your data frame
unique_families <- unique(family.df$Family)
unique_families
#choose colors
#mypalette_fam <- c('#1B9E77','#666666','#7570B3','#FFFF99','#66A61E','#E31A1C','#A6761D','#D95F02','#C43B80','#82906F','#DBA875' ,'#824CDE' , '#A6CEE3', '#74EA56', '#E05fAD', '#E6DC50','#27408B')
mypalette_fam <- c("#27408B","#E6DC50","#A6761D", "#74EA56", "#7570B3","#E05fAD","#A6CEE3", "#DBA875","#824CDE","#E31A1C", "#FFFF99", "#82906F",  "#C43B80", "#1B9E77", "#666666") 

# Assuming 'family.df' is your data frame
unique_families <- unique(family.df$Family)
# Assign your custom color palette to the unique families
myFamCol <- mypalette_fam
names(myFamCol) <- unique_families
# Use the pal.bands function to display the color palette
pal.bands(myFamCol, main = "Virus family")

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


right_ra=rowAnnotation("Virus family" = family.df$Family,  # Remove one value from the "Family" annotation
                       col=list('Virus family'=myFamCol),
                       show_annotation_name=T,
                       annotation_name_side = "top",
                       annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
                       annotation_name_rot = 360,
                       annotation_legend_param=list("Virus family"=list(title="Virus family", 
                                                                  labels_gp = gpar(fontface="italic"))))


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
                     right_annotation=right_ra,
                     row_names_gp = gpar(fontsize = 7),
                     show_column_names = T,
                     column_names_rot = 45,
                     column_names_gp=gpar(fontface=1, fontsize=6),
                     heatmap_legend_param = list(direction = "horizontal"),
                     cluster_rows = T,
                     row_title_gp = gpar(fontsize=10),
                     row_title_rot=0,
                     border=F)
#+ echo=TRUE, fig.width=15, fig.height=9
draw(hm, heatmap_legend_side = "left", annotation_legend_side = "left", 
     merge_legend = T, legend_grouping="original")


##' *Make treemap*
##' *Make treemap*
##' 
#Create a treemap 
### Make treemap (TM)
library(treemapify)
library(randomcoloR)
library(ggplot2)
library(treemap)
#'*1. Master table*
#'
# Get the taxonomy table from the phyloseq object
taxonomy_table <- tax_table(MOS20_final)
taxonomy_table <- as.data.frame(taxonomy_table)
# Extract row names and the "Family" column into a new table
Master_table <- data.frame(Contig = rownames(taxonomy_table), Family = taxonomy_table$Family)

# Extract the sample names from the contig names
#Master_table$Sample_Name <- sub(".*_(\\d+)$", "\\1", Master_table$Contig)
# Extract the contig names from the contig names
#Master_table$Contig_Name <- sub(".*\\s", "", Master_table$Contig)

# Extract the sample names from the contig names
# Extract the sample names from the Contig column
Master_table$Sample_Name <- sub(".*_P([^_]+)$", "P\\1", Master_table$Contig)
#sample_names <- sapply(strsplit(Master_table$Contig, "_"), function(x) paste(x[7:length(x)], collapse = "_"))
#Master_table$Sample_Name <- sample_names
#Master_table$Sample_Name <- sub(".*_([^_]+)_([^_]+)$", "\\1_\\2", Master_table$Contig) # will exctract sample name with one underscore
#Master_table$Sample_Name <- str_extract(Master_table$Contig, "(?<=_)[^_]+$") # will extract sample name with no underscore

#Extract coverage from contig name
library(stringr)
cov <- str_match(Master_table$Contig, "NODE_\\d+_length_\\d+_cov_([^.]+)")[, 2]
Master_table$Coverage <- cov


# Extract all the length from the contig names
contig_lengths <- str_match(Master_table$Contig, "length_(\\d+)_")[, 2]
Master_table$Contig_length <- as.numeric(contig_lengths)

# Calculate the total abundance per sample
total_abundance <- rowSums(otu_table(MOS20_final))
# Add the total abundance as a new column in the Master_table
Master_table$Total_Abundance <- total_abundance

# Assuming that the column in the `sample_data` containing mosquito species is named "Species"
sample_data <- sample_data(MOS20_final)
mosquito_species <- sample_data$Species
Mosquito_species <- sample_data(MOS20_final)$Species
# Match the sample names and merge the species information into the Master_table
Master_table$Mosquito_species <- Mosquito_species[match(Master_table$Sample_Name, rownames(sample_data))]

# Now Master_table contains the mosquito species information
# Get the mosquito species information from the metadata in the phyloseq object and add to master table
#sample_data <- sample_data(MOS20_final)
# Create a new column with mosquito species for each sample
#Master_table$Mosquito_Species <- sample_data$Species[match(Master_table$Sample_Name, sample_data$'Sample ID')]

# Assuming that the column in the `sample_data` containing mosquito Location

# Match the sample names and merge the Location information into the Master_table
#Location <- sample_data$Location
Location <- sample_data(MOS20_final)$Location
Master_table$Location <- Location[match(Master_table$Sample_Name, rownames(sample_data))]
# Create a new column with Location for each sample
#Master_table$Location <- sample_data$Location[match(Master_table$Sample_Name, sample_data$Sample_ID)]

# Assuming that the column in the `sample_data` containing virus name
Viral_species <- taxonomy_table$Species
# Match the sample names and merge the species information into the Master_table
Master_table$Viral_species <- Viral_species[match(Master_table$Contig, taxonomy_table$Contig)]
# Create a new column with contig virus species name
#contig_species <- taxonomy_table$Contig_Species[match(Master_table$Contig, taxonomy_table$Contig)]
# Match Species from taxonomy_table to row names of Master_table
Master_table$Viral_species <- Viral_species[match(Master_table$Contig, row.names(taxonomy_table))]

# Calculate the total abundance in percentage
Relative_abundance <- (Master_table$Total_Abundance / sum(Master_table$Total_Abundance)) * 100
# Add the total abundance percentage as a new column in the Master_table
Master_table$Relative_Abundance <- Relative_abundance

#'*2. Treemap*
#'
# Aggregate the Relative_Abundance values by Family
Aggregated_Family <- Master_table %>%
  group_by(Family) %>%
  summarise(Relative_Abundance = sum(Relative_Abundance))

# Sort the aggregated data by Relative_Abundance in descending order
Aggregated_Family <- Aggregated_Family %>%
  arrange(desc(Relative_Abundance))
Aggregated_Family



#Assign colors to  
display.brewer.all()
display.brewer.all(colorblindFriendly = TRUE) #Color blind friendly 
display.brewer.pal(n = 11, name = 'Paired') # Change 'n' and 'name' depending on the color palette 
# Hexadecimal color specification 
brewer.pal(n = 11, name = "Paired")

mypalette <- c('#1B9E77','#666666','#7570B3','#FFFF99','#66A61E','#E31A1C','#A6761D','#D95F02','#C43B80','#82906F','#DBA875' ,'#824CDE' , '#A6CEE3', '#74EA56', '#E05fAD', '#E6DC50','#27408B')
mypalette<- c("#27408B","#E6DC50","#A6761D", "#74EA56", "#7570B3","#E05fAD","#A6CEE3", "#DBA875","#824CDE","#E31A1C", "#FFFF99", "#82906F",  "#C43B80", "#1B9E77", "#666666") 

mypalette<- c("#C43B80", "#1B9E77", "#7570B3","#E31A1C", "#824CDE", "#A6761D", "#666666", "#82906F","#FFFF99","#A6CEE3","#E05fAD","#74EA56","#DBA875","#E6DC50","#27408B")

#Create a treemap 
### Make treemap (TM)
library(treemapify)
library(randomcoloR)
library(ggplot2)
library(treemap)
library(forcats)

# count the number of unique families
n_families <- Master_table %>%
  distinct(Family) %>%
  nrow()
n_families

n <- 15 #Number of families 

palette <- distinctColorPalette(n) #this will choose colors for you but the next line permits you choose your colors

TM <- ggplot(data = Aggregated_Family, aes(area = Relative_Abundance, fill = Family, label = paste0(Family, "\n"))) +
  geom_treemap() +
  geom_treemap_text(colour = "black", place = "centre") +
  labs(title = "Treemap of Relative Abundance") +
  theme(legend.position = "right")+
  scale_fill_manual(values = mypalette) +  # Custom color palette
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "Viral families",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic")))

TM


### Make treenap for only unclassified Family 
###
# Filter data to include only "Unclassified" rows
Unclassified_fam <- Master_table %>%
  filter(Family == "unclassified")
# Recalculate relative abundance of filtered data (unclassified) in percentage and add to new column
Relative_abundance_unclass <- (Unclassified_fam$Relative_Abundance / sum(Unclassified_fam$Relative_Abundance)) * 100
# Add the total abundance percentage as a new column in the Unclassified_fam data frame
Unclassified_fam$Relative_Abundance_unclass <- Relative_abundance_unclass
# Aggregate Viral_species values
Aggregated_Unclassified_fam <- Unclassified_fam %>%
  group_by(Viral_species) %>%
  summarise(Relative_Abundance_unclass = sum(Relative_Abundance_unclass))

# Reverse the order of the 'Viral_species' factor levels
#Aggregated_Unclassified_fam$Viral_species <- fct_reorder(Aggregated_Unclassified_fam$Viral_species, Aggregated_Unclassified_fam$Relative_Abundance_unclass)
Aggregated_Unclassified_fam$Viral_species <- fct_reorder(Aggregated_Unclassified_fam$Viral_species, -Aggregated_Unclassified_fam$Relative_Abundance_unclass)


n_families_unclass <- Aggregated_Unclassified_fam %>%
  distinct(Viral_species) %>%
  nrow()
n_families_unclass
n_fam <- 20
palette <- distinctColorPalette(n_fam)
# To choose your own colors 
palette <- c("#D86ED9", "#E04E8E", "#75E3A0", "#88954D", "#DDACA5", "#715CD1", "#D96255","#D4E691", "#AD39E5", "#80E7E1", "#976790", "#E6A764", "#DAA7DF", "#D0E1BD","#69AE9D", "#828ADD", "#DBE450", "#D3D9E2", "#7CB2D9", "#76E760")



TMU <- ggplot(data = Aggregated_Unclassified_fam, aes(area = Relative_Abundance_unclass, fill = Viral_species, label = paste0(Viral_species, "\n"))) +
  geom_treemap() +
  geom_treemap_text(colour = "black", place = "centre") +
  labs(title = "Treemap of Relative Abundance") +
  theme(legend.position = "right")+
  scale_fill_manual(values = palette) +  # Custom color palette
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title = "Viral families",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic")))

TMU







TMUnclass <- ggplot(Aggregated_Unclassified_fam, aes(area = Relative_Abundance_unclass, fill = Viral_species, label = paste0(Viral_species, "\n"))) +
  geom_treemap() +
  geom_treemap_text(colour = "black", place = "centre") +
  labs(title = "Treemap of Relative Abundance") +
  theme(legend.position = "right") +
  scale_fill_manual(values = palette) +  # Custom color palette
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic")))
TMUnclass

TMUnclass <- ggplot(data = Aggregated_Unclassified_fam, aes(area = Relative_Abundance_unclass, fill = Viral_species, label = paste0(Viral_species, "\n"))) +
  geom_treemap() +
  geom_treemap_text(colour = "black", place = "centre") +
  labs(title = "Treemap of Relative Abundance") +
  theme(legend.position = "right",
        plot.title = element_text(size = 16),  # Adjust title size
        plot.subtitle = element_text(size = 14),  # Adjust subtitle size
        legend.text = element_text(size = 10),  # Adjust legend text size
        legend.title = element_text(size = 12),  # Adjust legend title size
        plot.caption = element_text(hjust = 0)) +
  scale_fill_manual(values = palette) +  # Custom color palette
  guides(fill = guide_legend(
    label.position = "right",
    title = "Family unclassified viral species",  # Change the legend title here
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic",
                               size = 10),
    direction = "vertical",  # Set legend labels to be vertical
    nrow = 21))  # Adjust the number of rows in the legend

TMUnclass  # Print the modified plot

ggsave("treemap.png", plot = TMUnclass)


#' ## **Barplot Mos_spp*
#' ## **Barplot Mos_spp*
#' ## **Barplot Mos_spp*

# Create the bar plot based on mosquito spp
Mos_Fam_bp <-ggplot(Master_table, aes(x = Mosquito_species, fill = Family, y = Relative_Abundance)) +
  geom_bar(stat = "identity") +
  labs(x = "Mosquito Species", y = "Relative abundance(%)", fill = "Family") +
  theme_minimal()+
  scale_fill_manual(values = mypalette)+
  guides(fill = guide_legend(
    label.position = "right",
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic", size = 10)))+
  #geom_text(aes(label = paste0(round(Relative_abundance, 2), "%")), position = position_stack(vjust = 0.5), size = 3)+
  theme(axis.text.x = element_text(face = "italic")) +
  ylim(0, 100)

Mos_Fam_bp 











#' ## **Venn diagram*
#' ## **Venn diagram*
#' ## **Venn diagram*

#Mkae venn diagram
library(VennDiagram)
library(gridExtra)

# Create a list of families for each mosquito species
family_lists <- lapply(unique(Master_table$Mosquito_species), function(species) {
  Master_table$Family[Master_table$Mosquito_species == species]
})

# Create the Venn diagram with legends
Mos_venn <- venn.diagram(
  x = family_lists,
  category.names = unique(Master_table$Mosquito_species),
  filename = NULL,
  fill = c("#440099", "#fcc200"),
  alpha = 0.5,
  cex = 0.8,
  cat.col = c("#440099", "#fcc200"),
  cat.cex = 0.8,
  cat.fontface = "bold.italic",
  cat.dist = 0.08,
  cat.pos = 0,
  margin = 0.1,
  main = "Mosquito Species Venn Diagram",
  main.cex = 1.2,
  main.fontface = "bold",
  legend = list(
    title = "Mosquito Species",
    labels = unique(Master_table$Mosquito_Species),
    col = c("#440099", "#fcc200"),
    cex = 0.8,
    fontface = "bold",
    pt.cex = 1.2,
    pos = "bottom"))

# Display the Venn diagram
grid.newpage()
grid.draw(Mos_venn)

# Extract unique family names from Master_table
family_names <- unique(Master_table$Family)

# Create a title for the family column names
family_label_title <- textGrob(
  label = "Family Labels",
  gp = gpar(fontface = "bold", fontsize = 12))

# Add labels in a column on the right side
right_labels <- tableGrob(
  family_names,
  theme = ttheme_default(
    core = list(bg_params = list(fill = NA)),
    colhead = list(bg_params = list(fill = NA)),
    rowhead = list(bg_params = list(fill = NA))
  ),
  rows = NULL
)

# Display the Venn diagram with labels in a column on the right side
grid.newpage()
grid.draw(
  arrangeGrob(
    Mos_venn,
    right = right_labels,
    widths = c(2, 2)))

# Create a vector of colors for each family based on Mosquito_Species
family_colors <- sapply(family_names, function(family) {
  species <- unique(Master_table$Mosquito_species[Master_table$Family == family])
  if ("Ae. africanus" %in% species && "Ae. albopictus " %in% species) {
    return("black")    # Set color to black if family is present in both species
  } else if ("Ae. africanus" %in% species) {
    return("#330099")  # Set color for Ae. aegypti
  } else if ("Ae. albopictus " %in% species) {
    return("#fcc200")  # Set color for Ae. albopictus
  } else {
    return("gray")     # Set default color for other families
  }})

# Create a tableGrob with colored family names
colored_labels <- tableGrob(
  family_names,
  theme = ttheme_default(
    core = list(fg_params = list(col = family_colors), bg_params = list(fill = NA)),
    colhead = list(bg_params = list(fill = NA)),
    rowhead = list(bg_params = list(fill = NA))
  ),
  rows = NULL)


# Combine the title and the colored family labels using gridExtra
colored_labels_with_title <- arrangeGrob(
  family_label_title,
  colored_labels,
  nrow = 2,
  heights = unit(c(1, 10), c("lines", "null")))

# Display the Venn diagram with colored family names and the title in a column on the right side
grid.newpage()
grid.draw(
  arrangeGrob(
    Mos_venn,
    center = colored_labels_with_title,
    widths = c(1, 1)))  # Adjust the widths here (e.g., decrease the width of the right column)




#' ## **Barplot Mos_spp*
#' ## **Barplot Mos_spp*
#' ## **Barplot Mos_spp*

# Create the bar plot based on mosquito spp
Mos_Fam_bp <-ggplot(Master_table, aes(x = Mosquito_species, fill = Family, y = Relative_Abundance)) +
  geom_bar(stat = "identity") +
  labs(x = "Mosquito Species", y = "Relative abundance(%)", fill = "Family") +
  theme_minimal()+
  scale_fill_manual(values = mypalette)+
  guides(fill = guide_legend(
    label.position = "right",
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic", size = 10)))+
  #geom_text(aes(label = paste0(round(Relative_abundance, 2), "%")), position = position_stack(vjust = 0.5), size = 3)+
  theme(axis.text.x = element_text(face = "italic")) +
  ylim(0, 100)

Mos_Fam_bp 



#' ## **piechart Mos_spp*
#' ## **piechart Mos_spp*
#' ## **piechart Mos_spp*
# Calculate the total relative abundance for each mosquito species
Mos_Species_Total <- aggregate(Relative_Abundance ~ Mosquito_species, data = Master_table, sum)

# Create a pie chart for tatal abundance 
Mos_Species_pie <- ggplot(Mos_Species_Total, aes(x = "", y = Relative_Abundance, fill = Mosquito_species)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) + 
  labs(title = "Relative Abundance of Mosquito Species") +
  theme_void() +
  scale_fill_manual(values = mypalette)+
  geom_text(aes(label = paste0(round(Relative_Abundance, 2), "%")), position = position_stack(vjust = 0.5), size = 4)

Mos_Species_pie

# Create a pie chart for tatal abundance for every mosquito spp

Mos_Species_pie <- ggplot(Master_table, aes(x = "", y = Relative_Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) + 
  labs(title = "Relative Abundance of Mosquito Species") +
  theme_void() +
  scale_fill_manual(values = mypalette)

Mos_Species_pie


# Calculate the total relative abundance for each family within each species
Family_Total <- aggregate(Relative_Abundance ~ Family + Mosquito_species, data = Master_table, sum)

# Create separate pie charts for each species
p1 <- ggplot(subset(Family_Total, Mosquito_species == "Ae. africanus"), aes(x = "", y = Relative_Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) + 
  labs(title = "Ae. africanus: Relative Abundance of Families") +
  theme_void() +
  scale_fill_manual(values = palette) 
#geom_text(aes(label = paste0(round(Relative_Abundance, 2), "%")), position = position_stack(vjust = 0.5), size = 4)
p1

p2 <- ggplot(subset(Family_Total, Mosquito_species == "Ae. albopictus "), aes(x = "", y = Relative_Abundance, fill = Family)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) + 
  labs(title = "Ae. albopictus: Relative Abundance of Families") +
  theme_void() +
  scale_fill_manual(values = palette) 
  #geom_text(aes(label = paste0(round(Relative_Abundance, 2), "%")), position = position_stack(vjust = 0.5), size = 4)
p2

#TM
# Calculate the total relative abundance for each family within each species
Family_Total <- aggregate(Relative_Abundance ~ Family + Mosquito_species, data = Master_table, sum)


p2 <- ggplot(subset(Family_Total, Mosquito_species == "Ae. albopictus "), aes(area = Relative_Abundance, fill = Family), label = paste0(Family, "\n")) +
  geom_treemap()+
  #geom_treemap_text(colour = "black", place = "centre") +
  labs(title = "Treemap of Relative Abundance") +
  theme(legend.position = "right") +
  scale_fill_manual(values = palette) +  # Custom color palette
  theme(plot.caption = element_text(hjust = 0)) +
  guides(fill = guide_legend(
    label.position = "right",
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic")))
p2


#' ## **Alpha diversity**
#' ## **Alpha diversity**
#' ## **Alpha diversity**
#' 
##This method uses just asterisks 
A <- plot_richness(MOS20_final,x="Species",color="Species",measures=c("Observed", "Shannon","Simpson", "Chao1")) + 
  geom_boxplot()+
  geom_jitter(width=0.1, size=1)+
  stat_compare_means(aes(label= ..p.signif..), show.legend = F)+
  scale_color_discrete(labels=c(expression(paste(italic("Aedes africanus"))),
                                expression(paste(italic("Aedes albopictus")))))+
  theme_bw()+
  theme(axis.text.y=element_text(angle=90,hjust=0.5))+
  theme(strip.text.x = element_text(size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text.align = 0)+
  labs(color = "Mosquito species")
A

A_Final <- A + labs(caption = "Figure 3. Proportion of viral sequences relative to the number of mosquito host sequences. Wilcox Test: p < 0.05 (*), p < 0.01 (**), p < 0.001 (***), p < 0.0001 (****), ns (not significant)", size =0.1 ) +
  theme(plot.caption = element_text(size=9, hjust = 0))
A_Final
ggsave("~/Documents/MDKC Files/MOSCAM20/R/NGS/Heatmap_+_moscam21_species_identification (2)/R_sessions/R_MOS20/MOS20_images/Alpha_diversity.pdf", 
       plot = A_Final, width = 12, height = 8, units = "in", dpi = 300)

#'  *Beta diversity*
## 1. NMDS 
#' Just OTUs
M.ord <- ordinate(MOS20_final, "NMDS","bray")
pcoa1 <- plot_ordination(MOS20_final, M.ord, type = "taxa", color = "Genus", title = "taxa")
print(pcoa1)
pcoa1 + facet_wrap(~Family,5)

#'Just samples
M.ord <- ordinate(MOS20_final, "NMDS","bray")
M.ord
NMDS <- plot_ordination(MOS20_final, M.ord, type = "samples", color = "Location", shape = "Species") +
  scale_shape(labels=c(expression(paste(italic("Aedes africanus"))),
                       expression(paste(italic("Aedes albopictus")))))+
  theme_bw()+
  theme(legend.text.align = 0,
        legend.key.size = unit(1.5, "lines"))+
  #legend.text = element_text(face = "italic")) +
  geom_point(size = 2) +
  ggtitle("") + #Add title to plot 
  stat_ellipse(type = "norm", linetype = 3, aes_string(group ="Species"), show.legend = F ) #Adds circles 
NMDS

#Color Mosquito species
NMDS <- plot_ordination(MOS20_final, M.ord, type = "samples", color = "Species", shape = "Location") +
  geom_point(size = 2) +
  theme_bw()+
  #stat_ellipse(type = "norm", linetype = 2, aes_string(group="Species"), show.legend = F)+
  theme(legend.text.align = 0,
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(face = "italic")) +
  scale_shape_manual(values = 1:8)+
  geom_point(size = 2) +
  ggtitle("") +#Add title to plot a 
  stat_ellipse(type = "norm", linetype = 3, aes_string(group ="Species"), show.legend = F ) #Adds circles 

NMDS


# Calculate stress and R-square values
# Perform PERMANOVA
metadata_vegan <- as(sample_data(MOS20_final), "data.frame")
perm <- adonis2(distance(MOS20_final, method="bray") ~ Species*Location,
                data = metadata_vegan)
perm 

NMDS_caption <- NMDS + labs(caption = "Figure 2. Treemap plot of the proportion of reads recovered among non-host contigs assembled from the 54 pools mosquitoes") +
  theme(plot.caption = element_text(hjust = 0)). ## adds caption only to plot a
NMDS_caption

ggsave("~/Documents/MDKC Files/MOSCAM20/R/NGS/Heatmap_+_moscam21_species_identification (2)/R_sessions/R_MOS20/MOS20_images/MY-NMDS.pdf",
       plot = NMDS,  width = 12, height = 8, units = "in", dpi = 300)

## 2. PCoA
PCA <- ordinate(MOS20_final, "PCoA")
p <- plot_ordination(MOS20_final, PCA, type = "samples", color = "Species", shape = "Location") +
  scale_shape(labels=c(expression(paste(italic("Aedes aegypti"))),
                       expression(paste(italic("Aedes albopictus")))))+
  theme_bw()+
  theme(legend.text.align = 0) + xlab("PCoA_33.9 %") +
  ylab("PCoA_17.8 %") +
  ggtitle("b") #Add title to plot b

p




