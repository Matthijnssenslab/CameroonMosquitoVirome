library(ggforce)

#Load qPCR data
qMOS20 <- read.table("qMOS20_final.tsv", header=T, row.names=1, sep="\t", dec=".")

# Remove rows where Sample is "MOS20qEd20KC_087"
#qMOS20_table <- qMOS20[qMOS20$Sample != "AlEd20KC087", ]
qMOS20_table <- qMOS20[qMOS20$Sample_number != "87", ]

# Colors for mosquito species
species_colors <- c("Ae. africanus" = "#338068", "Ae. albopictus" = "#804C89", "Ae. simpsoni" = "#E181AC", "Ae. aegypti" = "#EF3924")



###'*Bafoussam mosquito solemovirus*
###'*Bafoussam mosquito solemovirus*

# Plot code with color assignment
Plot_GS <- ggplot(qMOS20_table, aes(x = Location, y = GS_copies, fill= Mosquito_species)) +
  geom_violin(aes(color = Mosquito_species, fill = Mosquito_species), scale = "width", alpha = 0.3, trim = FALSE , adjust = 1) +  
  labs(title = "Bafoussam mosquito solemovirus",
       y = "GS copies per mosquito") +  
  geom_sina(aes(color = Mosquito_species, fill = Mosquito_species), scale = TRUE, size = 2) +
  stat_summary(fun = "median", geom = "point", color = "red", size = 5, position = position_dodge(0.9), show.legend = FALSE) +
  geom_hline(yintercept = 1E+01, linetype = "dashed", color = "black") +  
  theme(
    axis.title = element_text(face = "bold", size = 16),  
    axis.text = element_text(size = 14),  
    panel.background = element_rect(fill = 'transparent'),
    plot.background = element_rect(fill = 'transparent', color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_rect(fill = 'transparent'),
    axis.line = element_line(color = "black"),  
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) +
  guides(fill = guide_legend(
    label.position = "right",
    title.position = "top",
    title.hjust = 0.5,
    label.hjust = 0,
    label.vjust = 0.5,
    label.theme = element_text(face = "italic", size = 10))) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0), add = c(0, 1)),
    trans = scales::pseudo_log_trans(base = 10), 
    breaks = c(0, 10^seq(1, 10, 1)),
    labels = scales::scientific_format(),
    limits = c(0, 1E+10)
  )

Plot_GS

# Calculate the number and percentage of positive samples for each Location and Mosquito species
positive_percent_GS <- qMOS20 %>%
  group_by(Location, Mosquito_species) %>%
  summarize(
    Positive_Samples = sum(GS_copies > 0),
    Percentage_Positive = (sum(GS_copies > 0) / n()) * 100
  )
positive_percent_GS

positive_percent_GS <- qMOS20 %>%
  group_by(Location, Mosquito_species) %>%
  summarize(
    Positive_Samples = sum(GS_copies > 0),
    Total_Tested = n(),
    Percentage_Positive = (Positive_Samples / Total_Tested) * 100
  )

positive_percent_GS

# Calculate mean and median of GS_copies column per Location and Mosquito species
mean_median_GS_copies <- qMOS20 %>%
  group_by(Location, Mosquito_species) %>%
  summarize(
    Mean_GS_copies = mean(GS_copies, na.rm = T),
    Median_GS_copies = median(GS_copies, na.rm = T)
  )

# Merge mean_median_GS_copies and positive_percent_GSV
qGS <- merge(mean_median_GS_copies, positive_percent_GS, by = c("Location", "Mosquito_species"), all = TRUE)
# Round all numerical columns to one decimal place
qGS[, -c(1, 2)] <- round(qGS[, -c(1, 2)], 1)
# Round all numerical columns to the nearest whole number
# qGS[, -c(1, 2)] <- ceiling(qGS[, -c(1, 2)])

# Convert Mean_GS_copies and Median_GS_copies back to numeric
qGS$Mean_GS_copies <- as.numeric(qGS$Mean_GS_copies)
qGS$Median_GS_copies <- as.numeric(qGS$Median_GS_copies)

# Convert Mean_GS_copies and Median_GS_copies to exponential form without any values after the decimal point
qGS$Mean_GS_copies <- sprintf("%.1e", qGS$Mean_GS_copies)
qGS$Median_GS_copies <- sprintf("%.1e", qGS$Median_GS_copies)

# Print the modified dataset
print(qGS)


###Circle heatmap
GS_HM <- positive_percent_GS %>%
  ggplot(aes(x = Location, y = Mosquito_species, fill = Percentage_Positive)) +
  geom_point(aes(size = Total_Tested), color = "black", shape = 21, stroke = 0.5) +
  scale_fill_gradient(low = "white", high = "#B21B20", name = "Prevalence (%)", limits = c(0, 100)) +
  scale_size_continuous(range = c(3, 35)) +  # Adjust the range as needed
  ggtitle("Bafoussam mosquito solemovirus") +
  labs(title = "Bafoussam mosquito solemovirus",
       x = "Location",
       y = "Mosquito Species") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "#E6E7E8", color = NA),
    axis.text = element_text(size = 16),
    axis.text.y = element_text(size = 16, face = "italic"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()  # Remove grid lines
  )

GS_HM

GS_HM <- positive_percent_GS %>%
  ggplot(aes(x = Location, y = Mosquito_species, fill = Percentage_Positive, label = Positive_Samples)) +
  geom_point(aes(size = Total_Tested), color = "black", shape = 21, stroke = 0.5) +
  geom_text(size = 3, color = "black", show.legend = FALSE) +  # Add text labels for Positive Samples
  scale_fill_gradient(low = "white", high = "#B21B20", name = "Prevalence (%)", limits = c(0, 100)) +
  scale_size_continuous(breaks = seq(0, 40, by = 10), range = c(3, 35)) +  # Adjust the breaks and range
  ggtitle("Bafoussam mosquito solemovirus") +
  labs(title = "Bafoussam mosquito solemovirus",
       x = "Location",
       y = "Mosquito Species") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = "#E6E7E8", color = NA),
    axis.text = element_text(size = 16),
    axis.text.y = element_text(size = 16, face = "italic"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid = element_blank()  # Remove grid lines
  )

GS_HM




###'*Bafoussam mosquito Solemovirus*
###'Replace GS with BMSV

###'*Bafoussam mosquito Rhabdovirus *
###'Replace GS with BMRV

###'*Bafoussam mosquito Bunyavirus 1*
###'Replace GS with BMBV1

###'*Bafoussam mosquito Bunyavirus 2*
###'Replace GS with BMBV2

###'*Bafoussam mosquito orthomyxovirus 1 *
###'Replace GS with BMOV1

###'*Bafoussam mosquito orthomyxovirus 2 *
###'Replace GS with BMOV2




