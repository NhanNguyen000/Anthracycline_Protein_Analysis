# R scripts
get.path <- function(list_categories) {
  require("here")
  library(here)
  path <- here::here(list_categories)
  return(path)
}
source(get.path("R_functions.R"))

## Proteomic data in vitro:
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

## WGCNA for proteomics data in invitro:
library(WGCNA)									# load WGCNA library
options(stringsAsFactors = FALSE) 						# Set the global option for all functions
WGCNA_Invitro_data   <- Filter_data(data = t(Invitro_data_log))
dev.off()
#Made_sample_tree(data = WGCNA_Invitro_data)
categorial_samples <- c(rep(1, 21), rep(2, 42), rep(3, 42), rep(4, 34))
op <- par(mar=c(14,4,4,2))
Made_sample_tree2(WGCNA_Invitro_data, categorial_samples)
rm(op)
dev.off()

Pick_WGCNA_power(data = WGCNA_Invitro_data)  # power = 2 --> R square = 0.64, mean connectivity = 59.6
Protein_Module_color <- Build_WGCNA_geneTree(data = WGCNA_Invitro_data, network_type = "signed", picked_power = 2)

MEList       <- moduleEigengenes(WGCNA_Invitro_data, Protein_Module_color$dynamicColors)	# Calculate eigengenes
MEs          <- MEList$eigengenes
get.MEtree_plot(MEs, 0.25) # no modules lower than threshold

get.module_protein_number(MEList)

Invitro_module_color_matrix           <- as.matrix(MEList$validColors)
rownames(Invitro_module_color_matrix) <- colnames(WGCNA_Invitro_data)
Invitro_module_color                  <- get.Module_proteins(Invitro_module_color_matrix)
write.table(Invitro_module_color$turquoise, "Invitro_module_turquoise_NN_2020Mar03.txt", col.names = FALSE, row.names = FALSE)
write.table(Invitro_module_color$blue, "Invitro_module_blue_NN_2020Mar03.txt", col.names = FALSE, row.names = FALSE)
write.table(Invitro_module_color$green, "Invitro_module_green_NN_2020Mar03.txt", col.names = FALSE, row.names = FALSE)


#Invitro_Hub_proteins     <- get.Hub_proteins(WGCNA_Invitro_data, MEList)
Invitro_High_MM_proteins <- get.High_ModuleMembership_proteins(WGCNA_Invitro_data, MEs, MEList, cutoff = 0.8)
write.table(rownames(Invitro_High_MM_proteins$turquoise), "Invitro_module_turquoise_HighMM_NN_2020Mar30.txt", col.names = FALSE, row.names = FALSE)
write.table(rownames(Invitro_High_MM_proteins$blue), "Invitro_module_blue_HighMM_NN_2020Mar30.txt", col.names = FALSE, row.names = FALSE)
write.table(rownames(Invitro_High_MM_proteins$green), "Invitro_module_green_HighMM_NN_2020Mar30.txt", col.names = FALSE, row.names = FALSE)


#Invitro_High_MM_proteins_no_missing_data <- list()
#for(module_color in names(Invitro_High_MM_proteins)) {
#  Invitro_High_MM_proteins_no_missing_data[[module_color]]<- get.protein_without_missing_value(WGCNA_Invitro_data, rownames(Invitro_High_MM_proteins[[module_color]]))
#}



# plot the module eigengenes
Protein_MEs      <- Module_Eigengen_Tree(data = WGCNA_Invitro_data, Module_color = MEList$validColors)
Protein_metadata <- as.data.frame(get.metadata(rownames(WGCNA_Invitro_data)))
rownames(Protein_metadata) <- rownames(WGCNA_Invitro_data)
Protein_metadata$Dose      <- paste0(Protein_metadata$Drug, "_", Protein_metadata$Dose)

Print_ME_in_pdf(MEs = Protein_MEs, metadata = Protein_metadata, 
                save_folder = get.path("outcome"), file_name = "WGCNA_Protein_Anthacycline_new_code_NN_2020Mar26.pdf")

## run PCA to identify the sample
library(factoextra)
setwd(get.path("outcome"))

pdf(file = "PCA_Invitro_NN_2020Mar26.pdf")
get.pca(WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")
get.pca_with_time(WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")
get.pca_with_selected_time (WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")

Invitro_module_color_matrix           <- as.matrix(MEList$validColors)
rownames(Invitro_module_color_matrix) <- colnames(WGCNA_Invitro_data)
Invitro_module_color                  <- get.Module_proteins(Invitro_module_color_matrix)

for (module in names(Invitro_module_color)) {
  data <- Select_if_in_specific_list(t(WGCNA_Invitro_data), 
                                     selected_list = Invitro_module_color[[module]])
  data <- t(data)
  plot(get.pca(data, substring(rownames(data), 1,8), name = module))
  plot(get.pca_with_time(data, substring(rownames(data), 1,8), name = module))
  plot(get.pca_with_selected_time(data, substring(rownames(data), 1,8), name = module))
}
dev.off()

# Print protein expression & Log2FC
#Print_protein_expression_log2FC_in_pdf(WGCNA_Invitro_data, Invitro_Hub_proteins, Protein_metadata, 
#                                       get.path("outcome"), "WGCNA_Hub_Protein_invitro_NN_2030Jan13.pdf")

## Proteomic data in biospies:
Biospies_data_log <- Loaded_data(folder = get.path("data"), file_name = "Hecatos_Px_Cardio_Biopsies_21samples_log2_median_norm.txt")
Biospies_data_log <- get.rowname(Biospies_data_log)
Biospies_data_log <- get.uniques_Uniprot(Biospies_data_log )

rownames(Biospies_data_log) <- get.Uniprot_IDs(rownames(Biospies_data_log))
Biospies_data_log           <- Round_number_in_data(data = Biospies_data_log , singificant_number = 6) 

Biospies_info <- read.table(paste0(get.path("data"), "/Hecatos_Px_Cardio_Biopsies_SampleAnnot.txt"), header = TRUE, sep = "\t", fill= TRUE)

Biospies_info$Biopsies_type <- get.biopsies_type(Biospies_info)
Sample_info                 <- Biospies_info[, c("SampleID", "Biopsies_type")]
Sample_info$SampleID        <- paste0("X", Sample_info$SampleID)
Biospies_data_log_withSampleInfo   <- merge(Sample_info, t(Biospies_data_log), by.x = "SampleID", by.y = "row.names")
rownames(Biospies_data_log_withSampleInfo) <- make.names(Biospies_data_log_withSampleInfo[, 2], unique = TRUE)

## WGCNA for proteomics data in biospies:
library(WGCNA)									# load WGCNA library
options(stringsAsFactors = FALSE) 						# Set the global option for all functions

WGCNA_Biospies_data    <- Filter_data(data = Biospies_data_log_withSampleInfo[, -c(1,2)])
dev.off()
#Made_sample_tree(data = WGCNA_Biospies_data)
categorial_samples <- c(1, 2, 3, rep(4, 6), 2, 3, 2, 4, rep(2, 3), rep(4, 2), 2, 4, 2)
op <- par(mar=c(14,4,4,2))
Made_sample_tree2(WGCNA_Biospies_data,categorial_samples)
rm(op)

Removed_samples <- c("Control_no_data", "Control.8", "Cardiotoxicity_with_ANT.6", "Control.9", "Cardiotoxicity_with_ANT.7")
Biospies_data_log_withSampleInfo_v2 <- Biospies_data_log_withSampleInfo[-which(rownames(Biospies_data_log_withSampleInfo) %in%Removed_samples),]

WGCNA_Biospies_data_v2 <- Filter_data(data = Biospies_data_log_withSampleInfo_v2[, -c(1,2)])
dev.off()
#Made_sample_tree(data = WGCNA_Biospies_data_v2)
categorial_samples <- c(rep(4, 5), 2, 3, 4, 2, 3, rep(2, 2), 4, rep(2,2), 4)
op <- par(mar=c(14,4,4,2))
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

#Try use eigenegen to categorise the biospies data:
condition_biospies <- unlist(lapply(strsplit(rownames(MEList_v2$eigengenes), "[.]"), `[[`, 1))
MEs_biopsies <- as.data.frame(cbind(condition_biospies, MEList_v2$eigengenes))

library("plyr")	
Mean_tem <- MEs_biopsies %>% dplyr::group_by(condition_biospies) %>% dplyr::summarise_all(.funs = c(mean = "mean", Sd="sd"), na.rm = TRUE)
Mean_tem <- get.rowname(as.matrix(Mean_tem))
class(Mean_tem) <- "numeric"
Mean_tem_v1 <- Mean_tem[, grep("mean", colnames(Mean_tem))]
op <- par(mar=c(14,4,4,2))
barplot(Mean_tem_v1, col=c("red", "green", "blue"),  las = 2, beside = TRUE)
barplot(Mean_tem_v1, col=c("red", "green", "blue"), legend = rownames(Mean_tem_v1),  las = 2, beside = TRUE)
rm(op)

Biopsies_selected_module <- c("turquoise", "blue", "brown", "yellow", "green")
# get protein from selected module
Biopsies_module_color_matrix           <- as.matrix(MEs_v2$moduleColors)
rownames(Biopsies_module_color_matrix) <- colnames(WGCNA_Biospies_data_v2)
Biopsies_module_color                  <- get.Module_proteins(Biopsies_module_color_matrix)
for(module in Biopsies_selected_module) {
  write.table(Biopsies_module_color[[module]], paste0("Biopsies_module_", module, "_NN_2020Mar30.txt"),
              col.names = FALSE, row.names = FALSE)
}
# get high module membership protein from selected module
#Biospies_Hub_proteins     <- chooseTopHubInEachModule(WGCNA_Biospies_data_v2, MEList_v2$validColors)
Biospies_High_MM_proteins <- get.High_ModuleMembership_proteins(WGCNA_Biospies_data_v2, MEs_v2$MEs, MEList_v2, cutoff = 0.8)
for(module in Biopsies_selected_module) {
  write.table(rownames(Biospies_High_MM_proteins[[module]]), paste0("Biopsies_module_", module, "_HighMM_NN_2020Mar30.txt"),
              col.names = FALSE, row.names = FALSE)
}

#Biospies_High_MM_proteins_no_missing_data <- list()
#for(module_color in names(Biospies_High_MM_proteins)) {
#  Biospies_High_MM_proteins_no_missing_data[[module_color]]<- get.protein_without_missing_value(WGCNA_Biospies_data_v2, rownames(Biospies_High_MM_proteins[[module_color]]))
#}

## run PCA to identify the sample
library(factoextra)
setwd(get.path("outcome"))
pdf(file = "PCA_Biopsies_NN_2020Mar30.pdf")

condition_biospies <- unlist(lapply(strsplit(rownames(WGCNA_Biospies_data_v2), "[.]"), `[[`, 1))
get.pca(WGCNA_Biospies_data_v2, condition_biospies, name = "Total Biopsies")

Biospies_module_color_matrix           <- as.matrix(MEs_v2$moduleColors)
rownames(Biospies_module_color_matrix) <- colnames(WGCNA_Biospies_data_v2)
Biospies_module_color                  <- get.Module_proteins(Biospies_module_color_matrix)

for (module in names(Biospies_module_color)) {
  data <- Select_if_in_specific_list(t(WGCNA_Biospies_data_v2), 
                                     selected_list = Biospies_module_color[[module]])
  plot(get.pca(t(data), condition_biospies, name = module))
}
dev.off()

## print protein expression
Print_protein_expression_log2FC_biopsies(WGCNA_Biospies_data_v2, Biospies_Hub_proteins) 
Print_protein_expression_log2FC_biopsies(WGCNA_Biospies_data_v2, c("Q9Y5L4", "Q07021", "P35232", "P15311", "P30044", "P14324")) 


## DisGeNET database:
Heart_Failure           <- Loaded_data(folder = get.path("data"), file_name = "C0018801_disease_gda_summary.tsv")
Heart_Failure_Evidences <- Loaded_data(folder = get.path("data"), file_name = "C0018801_disease_gda_evidences.tsv")
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


#Print_protein_expression_log2FC_in_pdf(WGCNA_Invitro_data, unlist(Protein_High_MM_opverlap), Protein_metadata, 
#                                       save_folder = get.path("outcome"), 
#                                       file_name = "WGCNA_Overlap_HighMM_Protein_invitro_NN_2020Jan13.pdf")
#Print_protein_expression_log2FC_biopsies(WGCNA_Biospies_data_v2, unlist(Protein_High_MM_opverlap)) 

# save data for cytoscape
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
output <- cbind(Proteins, Protein_detection, Heart_Failure_info_protein, High_MM_Protein)
write.csv(output, "Protein_inputCytoscape_NN_2020Mar31.csv")

dim(output) # total 1708 protein used for WGCNA analysis in invitro and biopsies
length(which(output$Protein_detection==3)) # 704 proteins detected in both invitro and biopsies

a<-output[which(output$Protein_detection==3),] # 704 proteins detected in both invitro and biopsies
dim(a)
b<-a[which(a$High_MM_Protein ==TRUE),] # 242 proteins are high module membership
dim(b)
length(which(b$Heart_Failure_info_protein == TRUE))

###================================= code not used
## Venn diagram
library(VennDiagram)
venn.diagram(
  x = list(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2), unique(Heart_Failure$UniProt)),
  fill = c("red", "orange", "blue"), sub.cex = 2,
  category.names = c("Invitro" , "Biopsies", "DisGeNet_Heart_Failure"),
  filename = 'Venn_diagramm_total_protein_NN_20120Jan13.png', lwd = 2, apha = 0.5,  output=FALSE
)

selected_module_inVitro <- c("turquoise", "green", "blue")
selected_module_biospies <- c("turquoise", "blue", "brown", "yellow", "green", "pink")
for(module_inVitro in selected_module_inVitro) {
  for(module_biospies in selected_module_biospies) {
    plist[i] <- venn.diagram(x = list(Invitro_module_color[[module_inVitro]], 
                                      Biospies_module_color[[module_biospies]],
                                      unique(Heart_Failure$UniProt)),
                             fill = c("red", "orange", "blue"), sub.cex = 2,
                             category.names = c(paste0("Invitro_", module_inVitro),
                                                paste0("Biopsies_", module_biospies),
                                                "DisGeNet_Heart_Failure"),
                             filename = paste0(module_inVitro, "_", module_biospies, "_DisGeNet", "_NN_2020Jan13.png"),
                             lwd = 2, apha = 0.5,  output=FALSE
    )
  }
}

# check the overlap between module
selected_module_inVitro <- c("turquoise", "blue", "green")
selected_module_biospies <- c("turquoise", "blue", "brown", "yellow", "green", "pink")
for(module_inVitro in selected_module_inVitro) {
  for(module_biospies in selected_module_biospies) {
    value <- intersect(Invitro_module_color[[module_inVitro]], Biospies_module_color[[module_biospies]])
    value2 <- intersect(value, unique(Heart_Failure$UniProt))
    #print(paste("Invitro_", module_inVitro, paste0("Biopsies_", module_biospies), length(value)))
    print(paste("Invitro_", module_inVitro, paste0("Biopsies_", module_biospies), length(value2)))
  }
}

for(module_biospies in selected_module_biospies) {
  value <- intersect(colnames(WGCNA_Invitro_data), Biospies_module_color[[module_biospies]])
  value2 <- intersect(value, unique(Heart_Failure$UniProt))
  #print(paste("Invitro_protein", paste0("Biopsies_", module_biospies), length(value)))
  print(paste("Invitro_protein", paste0("Biopsies_", module_biospies), length(value2)))
}

for(module_inVitro in selected_module_inVitro) {
  value <- intersect(Invitro_module_color[[module_inVitro]], colnames(WGCNA_Biospies_data_v2))
  value2 <- intersect(value, unique(Heart_Failure$UniProt))
  #print(paste("Invitro_", module_inVitro, "Biopsies_protein", length(value)))
  print(paste("Invitro_", module_inVitro,"Biopsies_protein", length(value2)))
}

value <- intersect(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2))
value2 <- length(intersect(value, unique(Heart_Failure$UniProt)))


## combine database
Proteins <- merge(Invitro_module_color_matrix, Biospies_module_color_matrix, by = "row.names", all = TRUE)
Proteins <- get.rowname(Proteins)
colnames(Proteins) <- c("Invitro_protein", "Biopsies_protein")

Protein_detected_invitro    <- rownames(Proteins) %in% colnames(WGCNA_Invitro_data)
Protein_detected_biopsies   <- rownames(Proteins) %in% colnames(WGCNA_Biospies_data_v2)
Protein_only_in_invitro     <- rownames(Proteins) %in% setdiff(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2))
Protein_only_in_biopsies    <- rownames(Proteins) %in% setdiff(colnames(WGCNA_Biospies_data_v2), colnames(WGCNA_Invitro_data))
Protein_in_both_invitro_biopsies <- rownames(Proteins) %in% intersect(colnames(WGCNA_Biospies_data_v2), colnames(WGCNA_Invitro_data))
Heart_Failure_info_protein  <- rownames(Proteins) %in% unique(Heart_Failure$UniProt)
High_MM_Protein_invitro     <- rownames(Proteins) %in% get.High_MM_proteins_list(Invitro_High_MM_proteins)
High_MM_Protein_biopsies    <- rownames(Proteins) %in% get.High_MM_proteins_list(Biospies_High_MM_proteins)
High_MM_Protein <- rownames(Proteins) %in% c(get.High_MM_proteins_list(Biospies_High_MM_proteins), get.High_MM_proteins_list(Invitro_High_MM_proteins))

output <- cbind(Proteins,
                Protein_detected_invitro, Protein_detected_biopsies,
                Protein_only_in_invitro, Protein_only_in_biopsies,
                Protein_in_both_invitro_biopsies, Heart_Failure_info_protein, 
                High_MM_Protein_invitro, High_MM_Protein_biopsies)
write.csv(output, "protein_data_NN_2020Jan13.csv")


Proteins <- merge(Invitro_module_color_matrix, Biospies_module_color_matrix, by = "row.names", all = TRUE)
Proteins <- get.rowname(Proteins)
colnames(Proteins) <- c("Invitro_protein", "Biopsies_protein")
Protein_detected_invitro    <- rownames(Proteins) %in% colnames(WGCNA_Invitro_data)
Protein_detected_biopsies   <- rownames(Proteins) %in% colnames(WGCNA_Biospies_data_v2)

Protein_detection <- c()
for(i in 1:length(rownames(Proteins))) {
  if (Protein_detected_invitro[i] == TRUE & Protein_detected_biopsies[i] == FALSE) Protein_detection[i]     <- 1
  else if (Protein_detected_invitro[i] == FALSE & Protein_detected_biopsies[i] ==TRUE) Protein_detection[i] <- 2
  else if (Protein_detected_invitro[i] == TRUE & Protein_detected_biopsies[i] ==TRUE) Protein_detection[i] <- 3
}

# 1 - Protein_detected_invitro, 2 - Protein_detected_biopsies, 3 - Protein_in_both_invitro_biopsies
Heart_Failure_info_protein  <- rownames(Proteins) %in% unique(Heart_Failure$UniProt)
High_MM_Protein <- rownames(Proteins) %in% c(get.High_MM_proteins_list(Biospies_High_MM_proteins), get.High_MM_proteins_list(Invitro_High_MM_proteins))

output <- cbind(Proteins, Protein_detection, Heart_Failure_info_protein, High_MM_Protein)
write.csv(output, "protein_data_v2_NN_2020Jan13.csv")

## Protein in both in vitro and biopsies data
#a<-output[which(output$Protein_detection == 3),] # 704 proteins detected in both in vitro and biopsies data
#b<-a[which(a$High_MM_Protein == TRUE),] # 242 proteins are high module membership
#d<-b[which(b$Heart_Failure_info_protein ==TRUE),] # 36 protein are high module membership

## Protein in DisGeNET
#a<-output[which(output$Heart_Failure_info_protein == TRUE),] # 169 proteins are in DisGeNET
#b<-a[which(a$High_MM_Protein == TRUE),] # 61 proteins are high module membership

## Protein high module membership
#a<-output[which(output$Heart_Failure_info_protein == TRUE),] # 477 proteins are in DisGeNET
#b<-a[which(a$Heart_Failure_info_protein == TRUE),]  # 61 proteins are high module membership

## protein in vitro
a<-output[which(output$Protein_detection != 2),] # 810 protien in vitro
b<-a[which(a$High_MM_Protein == TRUE),] # 247 proteins are high module membership
d<-b[which(b$Heart_Failure_info_protein == TRUE),] # 36 proteins in DIsGeNET

## protein in vitro
a<-output[which(output$Protein_detection != 1),] # 1602 protien in vitro
b<-a[which(a$High_MM_Protein == TRUE),] # 472 proteins are high module membership
d<-b[which(b$Heart_Failure_info_protein == TRUE),] # 61 proteins in DIsGeNET


## Check the overlap in proteomics data between in vitro vs. biospies:
library(VennDiagram)
venn.diagram(
  x = list(rownames(Invitro_data_log), rownames(Biospies_data_log)),
  category.names = c("ANT_Proteins" , "Bio_Proteins"),
  filename = 'Invitro_vs_Biopsies_NN_2020Jan13.png', lwd = 2, output= FALSE
)


library(VennDiagram)
venn.diagram(
  x = list(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2)),
  category.names = c("ANT_Proteins" , "Bio_Proteins"),
  filename = 'Invitro_vs_Biopsies_v2_NN_2020Jan13.png', lwd = 2, output= FALSE
)



library(VennDiagram)
venn.diagram(
  x = list(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2), unique(Heart_Failure$UniProt)),
  category.names = c("ANT_Proteins" , "Bio_Proteins", "DisGeNet_Proteins"),
  filename = 'Invitro_vs_Biopsies_v2_vs_DisGeNet_NN_2020Jan13.png', lwd = 3, output= FALSE
)

library(VennDiagram)
venn.diagram(
  x = list(Invitro_module_color$turquoise, colnames(WGCNA_Biospies_data_v2), unique(Heart_Failure$UniProt)),
  category.names = c("ANT_Proteins_Turquoise" , "Bio_Proteins", "DisGeNet_Proteins"),
  filename = 'Invitro_Turquoise_vs_Biopsies_v2_vs_DisGeNet_NN_2020Jan13.png', lwd = 3, output= FALSE
)

## Test

a<-intersect(Invitro_module_color$turquoise, colnames(WGCNA_Biospies_data_v2))
data <- Select_if_in_specific_list(t(WGCNA_Biospies_data_v2), selected_list = a)
heatmap(data)
data <- Select_if_in_specific_list(t(WGCNA_Invitro_data), selected_list = a)
heatmap(data)
write.table(a, row.names = FALSE, "test.txt")

a2<-intersect(a, unique(Heart_Failure$UniProt))
write.table(a2, row.names = FALSE, "test2.txt")


library(VennDiagram)
venn.diagram(
  x = list(Invitro_module_color$turquoise, Biospies_module_color$turquoise,
           unique(Heart_Failure$UniProt)),
  category.names = c("ANT_Proteins_Turquoise" , "Bio_Proteins_turquoise",
                     "DisGeNet_Proteins"),
  filename = '#18_venn_diagramm.png',
  output=TRUE
)

# run Veen diagram

library(VennDiagram)
for(module_inVitro in names(Invitro_module_color)) {
  for(module_biospies in names(Biospies_module_color)) {
    venn.diagram(x = list(Invitro_module_color[[module_inVitro]], 
                          Biospies_module_color[[module_biospies]],
                          unique(Heart_Failure$UniProt)),
                 category.names = c(paste0("ANT_Proteins_", module_inVitro),
                                    paste0("Bio_Proteins_", module_biospies),
                                    "DisGeNet_Proteins"),
                 filename = paste0(module_inVitro, "_", module_biospies, "_DisGeNet", ".png"),
                 output=TRUE
    )
  }
}


venn.diagram(x = list(Invitro_module_color$turquoise, rownames(Biospies_data_log), 
                      unique(Heart_Failure$UniProt), unique(Heart_Failure$UniProt)),
             category.names = c("ANT_Proteins_Turquoise", "Bio_Proteins", 
                                "HeartFailure_Proteins", "DisGeNET_Protein"),
             filename = "ANT_Proteins_Turquoise_HeartFailure.png",  output=TRUE)





