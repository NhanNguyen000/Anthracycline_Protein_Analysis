

# Load R functions --------------------------------------------------------

get.path <- function(list_categories) {
  require("here")
  library(here)
  path <- here::here(list_categories)
  return(path)
}
source(get.path("R_functions.R"))


# In vitro data -----------------------------------------------------------

## Clean in vitro data:
Invitro_data_log <- Loaded_data(folder = get.path("data"), file_name = "Hecatos_Cardio_Px_ANTvsFluctDMSO_log2_median_normalized.txt")
Invitro_data_log <- get.rowname(Invitro_data_log)
Invitro_data_log <- Invitro_data_log[, !colnames(Invitro_data_log) %in% c("Fluct.DMSO_T0_0959", "Fluct.DMSO_T0_0960", "Fluct.DMSO_T0_0958",
                                                                          "DOX_Therapeutic_T8_278.1", "DOX_Therapeutic_T8_279.1", "DOX_Therapeutic_T8_277.1")]
Invitro_data_log <- Invitro_data_log[, -grep("^DAU", colnames(Invitro_data_log))]

colnames(Invitro_data_log) <- Change_a_Part_in_Names(list_input = colnames(Invitro_data_log), select_chatacter = "[.]", wanted_character = "_")
Invitro_data_log           <- get.uniques_Uniprot(Invitro_data_log)
rownames(Invitro_data_log) <- get.Uniprot_IDs(rownames(Invitro_data_log))
colnames(Invitro_data_log) <- make.sample_names(colnames(Invitro_data_log))
Invitro_data_log           <- Round_number_in_data(data = Invitro_data_log, singificant_number = 6) 

# Run WGCNA for the in vitro data: ----------------------------------------
library(WGCNA)									# load WGCNA library
options(stringsAsFactors = FALSE) 						# Set the global option for all functions
WGCNA_Invitro_data   <- Filter_data(data = t(Invitro_data_log))
#Made_sample_tree(data = WGCNA_Invitro_data)
categorial_samples <- c(rep(1, 21), rep(2, 42), rep(3, 42), rep(4, 34))
op <- par(mar=c(4,1,1,12))
Made_sample_tree2(WGCNA_Invitro_data, categorial_samples)
rm(op)
dev.off()

Pick_WGCNA_power(data = WGCNA_Invitro_data)  # power = 2 --> R square = 0.64, mean connectivity = 59.6
Protein_Module_color <- Build_WGCNA_geneTree(data = WGCNA_Invitro_data, network_type = "signed", picked_power = 2)

MEList       <- moduleEigengenes(WGCNA_Invitro_data, Protein_Module_color$dynamicColors)	# Calculate eigengenes
MEs          <- MEList$eigengenes
get.MEtree_plot(MEs, 0.25) # no modules lower than threshold

get.module_protein_number(MEList)
Protein_MEs      <- Module_Eigengen_Tree(data = WGCNA_Invitro_data, Module_color = MEList$validColors)

Invitro_module_color_matrix           <- as.matrix(MEList$validColors)
rownames(Invitro_module_color_matrix) <- colnames(WGCNA_Invitro_data)
Invitro_module_color                  <- get.Module_proteins(Invitro_module_color_matrix)
write.table(Invitro_module_color$turquoise, "Invitro_module_turquoise_NN_2020Mar03.txt", col.names = FALSE, row.names = FALSE)
write.table(Invitro_module_color$blue, "Invitro_module_blue_NN_2020Mar03.txt", col.names = FALSE, row.names = FALSE)
write.table(Invitro_module_color$green, "Invitro_module_green_NN_2020Mar03.txt", col.names = FALSE, row.names = FALSE)

Invitro_High_MM_proteins <- get.High_ModuleMembership_proteins(WGCNA_Invitro_data, MEs, MEList, cutoff = 0.8)
write.table(rownames(Invitro_High_MM_proteins$turquoise), "Invitro_module_turquoise_HighMM_NN_2020Mar30.txt", col.names = FALSE, row.names = FALSE)
write.table(rownames(Invitro_High_MM_proteins$blue), "Invitro_module_blue_HighMM_NN_2020Mar30.txt", col.names = FALSE, row.names = FALSE)
write.table(rownames(Invitro_High_MM_proteins$green), "Invitro_module_green_HighMM_NN_2020Mar30.txt", col.names = FALSE, row.names = FALSE)



# Making figure for the in vitro data -------------------------------------
Protein_metadata <- as.data.frame(get.metadata(rownames(WGCNA_Invitro_data)))
rownames(Protein_metadata) <- rownames(WGCNA_Invitro_data)
Protein_metadata$Dose      <- paste0(Protein_metadata$Drug, "_", Protein_metadata$Dose)

# Figure S1:
Print_ME_in_pdf(MEs = Protein_MEs, metadata = Protein_metadata, 
                save_folder = get.path("outcome"), file_name = "FigureS1_2020Jul01.pdf")
# figure S2:
pdf("FigureS2_2020Jul01.pdf", width = 7, height = 12)

library(factoextra)
library("gridExtra")
setwd(get.path("outcome"))
#pca_total <- get.pca_with_selected_time(WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")
get.pca_with_selected_time(WGCNA_Invitro_data,
                           c(substring(rownames(WGCNA_Invitro_data)[1:21], 1,8), substring(rownames(WGCNA_Invitro_data)[22:length(rownames(WGCNA_Invitro_data))], 1,7)),
                           name = "Total Invitro")

Invitro_module_color_matrix           <- as.matrix(MEList$validColors)
rownames(Invitro_module_color_matrix) <- colnames(WGCNA_Invitro_data)
Invitro_module_color                  <- get.Module_proteins(Invitro_module_color_matrix)

plist <- list()
for (module in names(Invitro_module_color)) {
  data <- Select_if_in_specific_list(t(WGCNA_Invitro_data), 
                                     selected_list = Invitro_module_color[[module]])
  data <- t(data)
  plist[[module]] <- get.pca_with_selected_time(data, 
                                                c(substring(rownames(data)[1:21], 1,8), substring(rownames(data)[22:length(rownames(data))], 1,7)), 
                                                name = module)
}

ml <- marrangeGrob(plist, nrow=3, ncol=1)
print(ml)
dev.off()



# made figure 3 in paper:
# ME values
library("plyr")	
library("ggplot2")
library("gridExtra")
ME_tem    <- as.data.frame(cbind(Protein_metadata, setNames(Protein_MEs["MEturquoise"], "MEs")))
Mean_tem  <- ddply(ME_tem, ~Dose+Time, summarise, Mean = mean(MEs, na.rm = TRUE), sd = sd(MEs, na.rm = TRUE))
b1<-ggplot(Mean_tem, aes(x = Time, y = Mean, ymin = -0.25, ymax = 0.25, colour = Dose, group = Dose)) + 
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
  xlab("Time (hours)") + ylab("Eigengene values") + ggtitle("ME of turquoise") +  theme_bw() 
 
g <- ggplot_build(b1); unique(g$data[[1]]["colour"]) # exrtact colors list for dendrogram

ME_tem    <- as.data.frame(cbind(Protein_metadata, setNames(Protein_MEs["MEblue"], "MEs")))
Mean_tem  <- ddply(ME_tem, ~Dose+Time, summarise, Mean = mean(MEs, na.rm = TRUE), sd = sd(MEs, na.rm = TRUE))
b2<-ggplot(Mean_tem, aes(x = Time, y = Mean, ymin = -0.25, ymax = 0.25, colour = Dose, group = Dose)) + 
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
  xlab("Time (hours)") + ylab("Eigengene values") + ggtitle("ME of blue") +  theme_bw() 

ME_tem    <- as.data.frame(cbind(Protein_metadata, setNames(Protein_MEs["MEgreen"], "MEs")))
Mean_tem  <- ddply(ME_tem, ~Dose+Time, summarise, Mean = mean(MEs, na.rm = TRUE), sd = sd(MEs, na.rm = TRUE))
b3<-ggplot(Mean_tem, aes(x = Time, y = Mean, ymin = -0.25, ymax = 0.25, colour = Dose, group = Dose)) + 
  geom_line() + geom_point() + 
  geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
  xlab("Time (hours)") + ylab("Eigengene values") + ggtitle("ME of green") +  theme_bw() 

pdf("Figure3B_2020Jul01.pdf", width = 8, height = 2.5)
ml <- marrangeGrob(list(b1,b2), nrow=1, ncol=2)
print(ml)
dev.off()


#PCA
library(factoextra)
data <- t(Select_if_in_specific_list(t(WGCNA_Invitro_data), 
                                   selected_list = Invitro_module_color[["turquoise"]]))
c1<-plot(get.pca_with_selected_time(data, 
                                c(substring(rownames(data)[1:21], 1,8), substring(rownames(data)[22:length(rownames(data))], 1,7)), 
                                name = "turquoise"))

data <- t(Select_if_in_specific_list(t(WGCNA_Invitro_data), 
                                     selected_list = Invitro_module_color[["blue"]]))
c2<-plot(get.pca_with_selected_time(data, 
                                c(substring(rownames(data)[1:21], 1,8), substring(rownames(data)[22:length(rownames(data))], 1,7)), 
                                name = "blue"))
data <- t(Select_if_in_specific_list(t(WGCNA_Invitro_data), 
                                     selected_list = Invitro_module_color[["green"]]))
c3<-plot(get.pca_with_selected_time(data, 
                                c(substring(rownames(data)[1:21], 1,8), substring(rownames(data)[22:length(rownames(data))], 1,7)), 
                                name = "green"))
pdf("Figure3C_2020Jul01.pdf", width = 7, height = 12)
ml <- marrangeGrob(list(c1,c2,c3), nrow=3, ncol=1)
print(ml)
dev.off()

#sampling tree
dend <- as.dendrogram(hclust(dist(WGCNA_Invitro_data), method = "average"))
conditions <- c("Con", "DOX_The", "DOX_Tox", "EPI_The", "EPI_Tox", "IDA_The", "IDA_Tox")
plot_colors <- as.vector(unlist(unique(g$data[[1]]["colour"])))
plot_colors_order <-rep(NA, length(labels(dend)))
for(i in 1:length(conditions)) {
  plot_colors_order[grep(conditions[i], labels(dend))] <- plot_colors[i]
}

labels_colors(dend) <- plot_colors_order

pdf("Figure3A_2020Jul01.pdf", width = 7, height = 11)
op <- par(mar=c(4,1,1,12))
plot(dend, xlab = "Height", xlim = c(40, 0), horiz = TRUE)
rm(op)
dev.off()

# Biopsies data -----------------------------------------------------------

## Clean Biopsise data:
Biospies_data_log <- Loaded_data(folder = get.path("data"), file_name = "Hecatos_Px_Cardio_Biopsies_21samples_log2_median_norm.txt")
Biospies_data_log <- get.rowname(Biospies_data_log)
Biospies_data_log <- get.uniques_Uniprot(Biospies_data_log )

rownames(Biospies_data_log) <- get.Uniprot_IDs(rownames(Biospies_data_log))
Biospies_data_log           <- Round_number_in_data(data = Biospies_data_log , singificant_number = 6) 

Biospies_info <- read.table(paste0(get.path("data"), "/Hecatos_Px_Cardio_Biopsies_SampleAnnot.txt"), header = TRUE, sep = "\t", fill= TRUE)

Biospies_info$Biopsies_type <- get.biopsies_type(Biospies_info)
Biospies_info$Biopsies_type <- paste0(Biospies_info$Biopsies_type,"_", Biospies_info$SID)
Sample_info                 <- Biospies_info[, c("SampleID", "Biopsies_type")]
Sample_info$SampleID        <- paste0("X", Sample_info$SampleID)
Biospies_data_log_withSampleInfo   <- merge(Sample_info, t(Biospies_data_log), by.x = "SampleID", by.y = "row.names")
rownames(Biospies_data_log_withSampleInfo) <- make.names(Biospies_data_log_withSampleInfo[, 2], unique = TRUE)

## Making table S2:
write.csv(Biospies_info[c("SID", "SampleID", "Control...Cardiotoxicity",
                          "Sex", "Cancer.type", "Chemotherapeutic.agents")],
          file = "TableS2.csv", row.names = FALSE)

# Run WGCNA for the biopsies data -----------------------------------------

library(WGCNA)									# load WGCNA library
options(stringsAsFactors = FALSE) 						# Set the global option for all functions

WGCNA_Biospies_data    <- Filter_data(data = Biospies_data_log_withSampleInfo[, -c(1,2)])
dev.off()
#Made_sample_tree(data = WGCNA_Biospies_data)
categorial_samples <- c(1, 2, 3, rep(4, 6), 2, 3, 2, 4, rep(2, 3), rep(4, 2), 2, 4, 2)
op <- par(mar=c(4,4,4,16))
Made_sample_tree2(WGCNA_Biospies_data,categorial_samples)
rm(op)

Removed_samples <- c("Control_patient_10033",
                     "Patient_ANTtreatment_10027",
                     "Control_patient_10059",
                     "Patient_ANTtreatment_614",
                     "Control_patient_10021")
Biospies_data_log_withSampleInfo_v2 <- Biospies_data_log_withSampleInfo[-which(rownames(Biospies_data_log_withSampleInfo) %in%Removed_samples),]

WGCNA_Biospies_data_v2 <- Filter_data(data = Biospies_data_log_withSampleInfo_v2[, -c(1,2)])
dev.off()
#Made_sample_tree(data = WGCNA_Biospies_data_v2)
categorial_samples <- c(rep(4, 5), 2, 3, 4, 2, 3, rep(2, 2), 4, rep(2,2), 4)
op <- par(mar=c(4,4,4,16))
Made_sample_tree2(WGCNA_Biospies_data_v2,categorial_samples)
rm(op)
dev.off()

Pick_WGCNA_power(data = WGCNA_Biospies_data_v2)  # power = 5 --> R square = 0.73, mean connectivity = 56.9
Protein_Module_color   <- Build_WGCNA_geneTree(data = WGCNA_Biospies_data_v2, network_type = "signed", picked_power = 5)

MEList <- moduleEigengenes(WGCNA_Biospies_data_v2, Protein_Module_color$dynamicColors)	# Calculate eigengenes
MEs    <- MEList$eigengenes
get.MEtree_plot(MEs, 0.25) # 2 modules lower than threshold
MEs_v2  <- get.MEs_merge(WGCNA_Biospies_data_v2, Protein_Module_color$dynamicColors, MEDissThres = 0.25, Protein_Module_color$geneTree)
get.MEtree_plot(MEs_v2$MEs, 0.25)

MEList_v2                 <- moduleEigengenes(WGCNA_Biospies_data_v2, MEs_v2$moduleColors)
get.module_protein_number(MEList_v2) 

Biopsies_selected_module <- c("turquoise", "blue", "brown", "yellow", "green")
# get protein from selected module
Biopsies_module_color_matrix           <- as.matrix(MEs_v2$moduleColors)
rownames(Biopsies_module_color_matrix) <- colnames(WGCNA_Biospies_data_v2)
Biopsies_module_color                  <- get.Module_proteins(Biopsies_module_color_matrix)
for(module in Biopsies_selected_module) {
  write.table(Biopsies_module_color[[module]], paste0("Biopsies_module_", module, "_NN_2020Mar30.txt"),
              col.names = FALSE, row.names = FALSE)
}

Biospies_High_MM_proteins <- get.High_ModuleMembership_proteins(WGCNA_Biospies_data_v2, MEs_v2$MEs, MEList_v2, cutoff = 0.8)
for(module in Biopsies_selected_module) {
  write.table(rownames(Biospies_High_MM_proteins[[module]]), paste0("Biopsies_module_", module, "_HighMM_NN_2020Mar30.txt"),
              col.names = FALSE, row.names = FALSE)
}

# Making figure for the biopsies data -------------------------------------

# plot figure S3:
pdf("FigureS3.pdf", width = 12, height = 6)
dend <- as.dendrogram(hclust(dist(WGCNA_Biospies_data), method = "average"))
conditions <- c("Control", "Patient_ANT", "Patient_nonANT")
plot_colors <- c("#F8766D", "#C49A00", "#53B400")
plot_colors_order <-rep(NA, length(labels(dend)))
for(i in 1:length(conditions)) {
  plot_colors_order[grep(conditions[i], labels(dend))] <- plot_colors[i]
}

labels_colors(dend) <- plot_colors_order
par(mar=c(4,4,4,16))
plot(dend, xlab = "Height", xlim=c(140,0), horiz = TRUE)
dev.off()

# figure S4:
pdf("FigureS4_2020Jul01.pdf", width = 21, height = 12)

library(factoextra)
library("gridExtra")
setwd(get.path("outcome"))
plist <- list()
condition_biospies <- paste0(unlist(lapply(strsplit(rownames(MEList_v2$eigengenes), "[_]"), `[[`, 1)),
                             "_", unlist(lapply(strsplit(rownames(MEList_v2$eigengenes), "[_]"), `[[`, 2)))
plist[["pca_total"]] <- get.pca(WGCNA_Biospies_data_v2, condition_biospies, name = "total biopsies")

Biospies_module_color_matrix           <- as.matrix(MEs_v2$moduleColors)
rownames(Biospies_module_color_matrix) <- colnames(WGCNA_Biospies_data_v2)
Biospies_module_color                  <- get.Module_proteins(Biospies_module_color_matrix)

for (module in names(Biospies_module_color)) {
  data <- Select_if_in_specific_list(t(WGCNA_Biospies_data_v2), 
                                     selected_list = Biospies_module_color[[module]])
  plist[[module]] <- get.pca(t(data), condition_biospies, name = module)
}
ml <- marrangeGrob(plist, nrow=3, ncol=3)
print(ml)
dev.off()


# plot figure 4:

pdf("Figure4.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
#sampling 
dend <- as.dendrogram(hclust(dist(WGCNA_Biospies_data_v2), method = "average"))
conditions <- c("Control", "Patient_ANT", "Patient_nonANT")
plot_colors <- c("#F8766D", "#C49A00", "#53B400")
plot_colors_order <-rep(NA, length(labels(dend)))
for(i in 1:length(conditions)) {
  plot_colors_order[grep(conditions[i], labels(dend))] <- plot_colors[i]
}

labels_colors(dend) <- plot_colors_order
par(mar=c(4,4,4,16))
plot(dend, xlab = "Height",  horiz = TRUE)

# ME values:
condition_biospies <- paste0(unlist(lapply(strsplit(rownames(MEList_v2$eigengenes), "[_]"), `[[`, 1)),
                             "_", unlist(lapply(strsplit(rownames(MEList_v2$eigengenes), "[_]"), `[[`, 2)))
MEs_biopsies <- as.data.frame(cbind(condition_biospies, MEList_v2$eigengenes))

library("plyr")	
Mean_tem <- MEs_biopsies %>% dplyr::group_by(condition_biospies) %>% dplyr::summarise_all(.funs = c(mean = "mean", Sd="sd"), na.rm = TRUE)
Mean_tem <- get.rowname(as.matrix(Mean_tem))
class(Mean_tem) <- "numeric"
Mean_tem_v1 <- Mean_tem[, grep("mean", colnames(Mean_tem))]
biopsies_modules <-unlist(lapply(strsplit(colnames(Mean_tem_v1), "[_]"), `[[`, 1))
biopsies_modules<- substring(biopsies_modules, 3)
colnames(Mean_tem_v1) <- paste0("ME of ", biopsies_modules)
par(mar=c(9,4,4,4))
barplot(Mean_tem_v1, col= plot_colors, 
        ylab = "Eigengene values", ylim = c(-0.15, 0.20), las = 2, beside = TRUE)
legend("topright", legend=rownames(Mean_tem_v1), fill= plot_colors)
dev.off()


# Make figure 6 + fingure S6: combine both in vitro and biopsies data ------------------
# in vitro
selected_Proteins <- c("Q16698", "P06576", "P38117", "P29692") # DECR1, ATP5F1B, ETFB, EEF1D
outcome1 <- Print_protein_expression_log2FC_in_pdf(WGCNA_Invitro_data, selected_Proteins, Protein_metadata, 
                                       get.path("outcome"), "Figure6A_invitro_2020Jul01.pdf")

selected_Proteins <- c("P11021", "P04792") # HSPA5, HSPB1
outcome2 <- Print_protein_expression_log2FC_in_pdf(WGCNA_Invitro_data, selected_Proteins, Protein_metadata, 
                                       get.path("outcome"), "Figure6A_biopsies_2020Jul01.pdf")

# biopsies
pdf("Figure6B_2020Jul01.pdf", width = 8, height = 6)
selected_Proteins <- c("Q16698", "P06576", "P38117", "P29692", "P11021", "P04792") # DECR1, ATP5F1B, ETFB, EEF1D, HSPA5, HSPB1
outcome3 <- Print_protein_expression_log2FC_biopsies(WGCNA_Biospies_data_v2, selected_Proteins) 
dev.off()


# Using DisGeNET database -------------------------------------------------

#Heart_Failure2           <- Loaded_data(folder = get.path("data"), file_name = "C0018801_disease_gda_summary.tsv") # Disgenet data 2019
#Heart_Failure_Evidences <- Loaded_data(folder = get.path("data"), file_name = "C0018801_disease_gda_evidences.tsv") #  Disgenet data 2019

Heart_Failure           <- Loaded_data(folder = get.path("data"), file_name = "C0018801_disease_gda_summary_2020Aug26.tsv") # Disgenet data 2020
length(unique(Heart_Failure$UniProt))
Heart_Failure_Evidences <- Loaded_data(folder = get.path("data"), file_name = "C0018801_disease_gda_evidences_2020Aug26.tsv") # Disgenet data 2020
Heart_Failure_Evidences <- Heart_Failure_Evidences[Heart_Failure_Evidences$Gene %in% Heart_Failure$Gene, ] 

Heart_Failure_Type_of_Evidences <- list()
Heart_Failure_Type_of_Evidences[["All_proteins"]] <- unique(Heart_Failure$UniProt)
Heart_Failure_Type_of_Evidences[["All_genes"]]    <- unique(Heart_Failure$Gene)
Heart_Failure_Type_of_Evidences[["All_genes_have_evidence"]] <- intersect(Heart_Failure_Evidences$Gene, Heart_Failure$Gene)
for (evidence_type in unique(Heart_Failure_Evidences$Association_Type)) {
  Heart_Failure_Type_of_Evidences[[evidence_type]] <- unique(Heart_Failure_Evidences$Gene[Heart_Failure_Evidences$Association_Type == evidence_type])
}

Invitro_HeartFailure_Proteins  <- Select_if_in_specific_list(t(WGCNA_Invitro_data), 
                                                             selected_list = unique(Heart_Failure$UniProt))
Invitro_HeartFailure_Proteins_noMisingValue <- get.protein_without_missing_value(WGCNA_Invitro_data, 
                                                                                 rownames(Invitro_HeartFailure_Proteins))
Module_Invitro_ovelap_Heart_Failure <- get.Module_Proteins_overlap_DisGeNET(Invitro_module_color, 
                                                                            Heart_Failure, 
                                                                            WGCNA_Invitro_data)

Biospies_HeartFailure_Proteins <- Select_if_in_specific_list(t(WGCNA_Biospies_data_v2), 
                                                             selected_list = unique(Heart_Failure$UniProt))
Biospies_HeartFailure_Proteins_noMisingValue <- get.protein_without_missing_value(WGCNA_Biospies_data_v2, 
                                                                                  rownames(Biospies_HeartFailure_Proteins))
Module_Biospies_ovelap_Heart_Failure <- get.Module_Proteins_overlap_DisGeNET(Biospies_module_color, 
                                                                             Heart_Failure, 
                                                                             WGCNA_Biospies_data_v2)

# module high MM vs. Heart failure
Invitro_HeartFailure_HighMM <-list()
for (module_inVitro in names(Invitro_High_MM_proteins)) {
  Invitro_HeartFailure_HighMM[[module_inVitro]]<- intersect(rownames(Invitro_High_MM_proteins[[module_inVitro]]), 
                                                              Module_Invitro_ovelap_Heart_Failure$Protein_overlap_DisGeNet[[module_inVitro]])
}

Biopsies_HeartFailure_HighMM <-list()
for (module_biospies in names(Biospies_High_MM_proteins)) {
  Biopsies_HeartFailure_HighMM[[module_biospies]]<- intersect(rownames(Biospies_High_MM_proteins[[module_biospies]]), 
                                                                Module_Biospies_ovelap_Heart_Failure$Protein_overlap_DisGeNet[[module_biospies]])
}

Protein_HighMM_overlap<-list()
for (module_inVitro in names(Invitro_High_MM_proteins)) {
  for (module_biospies in names(Biospies_High_MM_proteins)) {
    names_module <- paste0(module_inVitro, "_", module_biospies)
    Protein_HighMM_overlap[[names_module]] <- intersect(rownames(Invitro_High_MM_proteins[[module_inVitro]]), 
                                                          rownames(Biospies_High_MM_proteins[[module_biospies]]))
  }
}
unlist(Protein_HighMM_overlap)
intersect(unlist(Protein_HighMM_overlap), unique(Heart_Failure$UniProt))


# Making data input for cytoscape -----------------------------------------

Proteins <- merge(Invitro_module_color_matrix, Biospies_module_color_matrix, by = "row.names", all = TRUE)
Proteins <- get.rowname(Proteins)
colnames(Proteins) <- c("Invitro_protein", "Biopsies_protein")

# 1 - Protein_detected_invitro, 2 - Protein_detected_biopsies, 3 - Protein_in_both_invitro_biopsies
Protein_detected_invitro    <- rownames(Proteins) %in% colnames(WGCNA_Invitro_data)
Protein_detected_biopsies   <- rownames(Proteins) %in% colnames(WGCNA_Biospies_data_v2)
Protein_detection <- c()
for(i in 1:length(rownames(Proteins))) {
  if (Protein_detected_invitro[i] == TRUE & Protein_detected_biopsies[i] == FALSE) Protein_detection[i]     <- 1
  else if (Protein_detected_invitro[i] == FALSE & Protein_detected_biopsies[i] ==TRUE) Protein_detection[i] <- 2
  else if (Protein_detected_invitro[i] == TRUE & Protein_detected_biopsies[i] ==TRUE) Protein_detection[i] <- 3
}

Heart_Failure_info_protein  <- rownames(Proteins) %in% unique(Heart_Failure$UniProt)
High_MM_Protein <- rownames(Proteins) %in% c(get.High_MM_proteins_list(Biospies_High_MM_proteins), 
                                             get.High_MM_proteins_list(Invitro_High_MM_proteins))
High_MM_Protein_overlap <- rownames(Proteins) %in% unlist(Protein_HighMM_overlap)
output <- cbind(Proteins, Protein_detection, Heart_Failure_info_protein, High_MM_Protein, High_MM_Protein_overlap)

output[which(rownames(Proteins) %in% unlist(Protein_HighMM_overlap)),] 
write.csv(output, "Protein_inputCytoscape_NN_2020Aug26.csv")

dim(output) # total 1708 protein used for WGCNA analysis in invitro and biopsies
length(which(output$Protein_detection==3)) # 704 proteins detected in both invitro and biopsies

a<-output[which(output$Protein_detection==3),] # 704 proteins detected in both invitro and biopsies
dim(a)
b<-a[which(a$High_MM_Protein ==TRUE),] # 242 proteins are high module membership
dim(b)
length(which(b$Heart_Failure_info_protein == TRUE))
