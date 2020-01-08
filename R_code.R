## Proteomic data in vitro:
Invitro_data_log <- Loaded_data(folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", file_name = "Hecatos_Cardio_Px_ANTvsFluctDMSO_log2_median_normalized.txt")
Invitro_data_log <- get.rowname(Invitro_data_log)
Invitro_data_log <- Invitro_data_log[, !colnames(Invitro_data_log) %in% c("Fluct.DMSO_T0_0959", "Fluct.DMSO_T0_0960", "Fluct.DMSO_T0_0958",
                                                                          "DOX_Therapeutic_T8_278.1", "DOX_Therapeutic_T8_279.1", "DOX_Therapeutic_T8_277.1")]
Invitro_data_log <- Invitro_data_log[, -grep("^DAU", colnames(Invitro_data_log))]

colnames(Invitro_data_log)       <- Change_a_Part_in_Names(list_input = colnames(Invitro_data_log), select_chatacter = "[.]", wanted_character = "_")
#colnames(Invitro_data_log)       <- gsub('.{4}$', '', colnames(Invitro_data_log))
#colnames(Invitro_data_log)[1:21] <- gsub('.{1}$', '', colnames(Invitro_data_log)[1:21])

Removed_Proteins           <- Identify_if_have_speficic_symbol(data = rownames(Invitro_data_log), symbol = ":")
Invitro_data_log           <- Remove_if_have_specific_symbol(data = Invitro_data_log, remove_list = Removed_Proteins)
Invitro_data_log           <- Round_number_in_data(data = Invitro_data_log, singificant_number = 6) 
rownames(Invitro_data_log) <- Select_Part_Names(rownames(Invitro_data_log), separator = "|", selected_position = 2)
colnames(Invitro_data_log) <- make.sample_names(colnames(Invitro_data_log))

## WGCNA for proteomics data in invitro:
library(WGCNA)									# load WGCNA library
options(stringsAsFactors = FALSE) 						# Set the global option for all functions
WGCNA_Invitro_data   <- Filter_data(data = t(Invitro_data_log))
#Made_sample_tree(data = WGCNA_Invitro_data)
categorial_samples <- c(rep(1, 21), rep(2, 42), rep(3, 42), rep(4, 34))
op <- par(mar=c(14,4,4,2))
Made_sample_tree2(WGCNA_Invitro_data, categorial_samples)
rm(op)
Pick_WGCNA_power(data = WGCNA_Invitro_data)  # power = 2 --> R square = 0.64, mean connectivity = 59.6
Protein_Module_color <- Build_WGCNA_geneTree(data = WGCNA_Invitro_data, network_type = "signed", picked_power = 2)
# https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/


MEList       <- moduleEigengenes(WGCNA_Invitro_data, Protein_Module_color$dynamicColors)	# Calculate eigengenes
MEs          <- MEList$eigengenes
get.MEtree_plot(MEs, 0.25) # no modules lower than threshold
#MEs_v2       <- get.MEs_merge(WGCNA_Invitro_data, Protein_Module_color$dynamicColors, MEDissThres = 0.25, Protein_Module_color$geneTree)
#get.MEtree_plot(MEs_v2$MEs, 0.25)
#MEList_v2                <- moduleEigengenes(WGCNA_Invitro_data, MEs_v2$moduleColors)

Invitro_Hub_proteins     <- get.Hub_proteins(WGCNA_Invitro_data, MEList)
Invitro_High_MM_proteins <- get.High_ModuleMembership_proteins(WGCNA_Invitro_data, MEs, MEList, cutoff = 0.8)
Invitro_High_MM_proteins_no_missing_data <- list()
for(module_color in names(Invitro_High_MM_proteins)) {
  Invitro_High_MM_proteins_no_missing_data[[module_color]]<- get.protein_without_missing_value(WGCNA_Invitro_data, rownames(Invitro_High_MM_proteins[[module_color]]))
}

Invitro_module_color_matrix           <- as.matrix(MEList$validColors)
rownames(Invitro_module_color_matrix) <- colnames(WGCNA_Invitro_data)
Invitro_module_color                  <- get.Module_proteins(Invitro_module_color_matrix)


# plot the module eigengenes
Protein_MEs<-Module_Eigengen_Tree(data = WGCNA_Invitro_data, Module_color = MEList$validColors)
Protein_metadata <- as.data.frame(get.metadata(rownames(WGCNA_Invitro_data)))
rownames(Protein_metadata) <- rownames(WGCNA_Invitro_data)
Protein_metadata$Dose <- paste0(Protein_metadata$Drug, "_", Protein_metadata$Dose)

#Print_ME_in_pdf(MEs = Protein_MEs, metadata = Protein_metadata, 
#                save_folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", file_name = "WGCNA_Protein_Anthacycline_new_code_NN_20191216.pdf")

## run PCA to identify the sample
library(factoextra)
pdf(file = "PCA_Invitro_NN_20191216.pdf")
get.pca(WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")
get.pca_with_time(WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")
get.pca_with_selected_time (WGCNA_Invitro_data, substring(rownames(WGCNA_Invitro_data), 1,8), name = "Total Invitro")

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
Print_protein_expression_in_pdf(WGCNA_Invitro_data, Invitro_Hub_proteins, Protein_metadata, 
                                save_folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", 
                                file_name = "WGCNA_Hub_Protein_invitro_NN_20191217.pdf")

Print_protein_expression_log2FC_in_pdf(WGCNA_Invitro_data, Invitro_Hub_proteins, Protein_metadata, 
                                       save_folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", 
                                       file_name = "WGCNA_Hub_Protein_invitro_NN_20191218.pdf")

## Proteomic data in biospies:
Biospies_data_log <- Loaded_data("I:/HeCaToS/HeCaToS_info_Data/TGXdata/Human_Biopsies", file_name = "Hecatos_Px_Cardio_Biopsies_21samples_log2_median_norm.txt")
Biospies_data_log <- get.rowname(Biospies_data_log)

Removed_Proteins             <- Identify_if_have_speficic_symbol(data = rownames(Biospies_data_log), symbol = ":")
Biospies_data_log            <- Remove_if_have_specific_symbol(data = Biospies_data_log, remove_list = Removed_Proteins)
Biospies_data_log            <- Round_number_in_data(data = Biospies_data_log , singificant_number = 6) 
rownames(Biospies_data_log ) <- Select_Part_Names(rownames(Biospies_data_log ), separator = "|", selected_position = 2)

Biospies_info<- read.table("I:/HeCaToS/HeCaToS_info_Data/TGXdata/Human_Biopsies/Hecatos_Px_Cardio_Biopsies_SampleAnnot.txt",
                           header = TRUE, sep = "\t", fill= TRUE)
#Biospies_info<- read.csv("I:/HeCaToS/HeCaToS_info_Data/TGXdata/Human_Biopsies/Hecatos_Px_Cardio_Biopsies_SampleAnnot_NN_20191209.csv")

Biospies_info$Biopsies_type <- get.biopsies_type(Biospies_info)
Sample_info                 <- Biospies_info[, c("SampleID", "Biopsies_type")]
Sample_info$SampleID        <- paste0("X", Sample_info$SampleID)
Biospies_data_log_rotate    <- merge(Sample_info, t(Biospies_data_log), by.x = "SampleID", by.y = "row.names")
rownames(Biospies_data_log_rotate) <- make.names(Biospies_data_log_rotate[, 2], unique = TRUE)


## WGCNA for proteomics data in biospies:
library(WGCNA)									# load WGCNA library
options(stringsAsFactors = FALSE) 						# Set the global option for all functions

#WGCNA_Biospies_data    <- Filter_data(data = t(Biospies_data_log))
#Made_sample_tree(data = WGCNA_Biospies_data)
#WGCNA_Biospies_data    <- WGCNA_Biospies_data[-grep("^X2017", rownames(WGCNA_Biospies_data)),]

WGCNA_Biospies_data    <- Filter_data(data = Biospies_data_log_rotate[, -c(1,2)])
#Made_sample_tree(data = WGCNA_Biospies_data)
categorial_samples <- c(1, 2, 3, rep(4, 6), 2, 3, 2, 4, rep(2, 3), rep(4, 2), 2, 4, 2)
op <- par(mar=c(14,4,4,2))
Made_sample_tree2(WGCNA_Biospies_data,categorial_samples)
rm(op)

Removed_samples        <- c("Control_no_data", "Control.8", "Cardiotoxicity_with_ANT.6", 
                            "Control.9", "Cardiotoxicity_with_ANT.7")
Biospies_data_log_rotate_v2 <- Remove_if_have_specific_symbol(data = Biospies_data_log_rotate, remove_list = Removed_samples)


WGCNA_Biospies_data_v2 <- Filter_data(data = Biospies_data_log_rotate_v2[, -c(1,2)])
#Made_sample_tree(data = WGCNA_Biospies_data_v2)
categorial_samples <- c(rep(4, 5), 2, 3, 4, 2, 3, rep(2, 2), 4, rep(2,2), 4)
op <- par(mar=c(14,4,4,2))
Made_sample_tree2(WGCNA_Biospies_data_v2,categorial_samples)
rm(op)

Pick_WGCNA_power(data = WGCNA_Biospies_data_v2)  # power = 5 --> R square = 0.73, mean connectivity = 56.9
Protein_Module_color   <- Build_WGCNA_geneTree(data = WGCNA_Biospies_data_v2, network_type = "signed", picked_power = 5)

MEList <- moduleEigengenes(WGCNA_Biospies_data_v2, Protein_Module_color$dynamicColors)	# Calculate eigengenes
MEs    <- MEList$eigengenes
get.MEtree_plot(MEs, 0.25) # 2 modules lower than threshold
MEs_v2  <- get.MEs_merge(WGCNA_Biospies_data_v2, Protein_Module_color$dynamicColors, MEDissThres = 0.25, Protein_Module_color$geneTree)
get.MEtree_plot(MEs_v2$MEs, 0.25)


MEList_v2                 <- moduleEigengenes(WGCNA_Biospies_data_v2, MEs_v2$moduleColors)
#Biospies_Hub_proteins     <- get.Hub_proteins(WGCNA_Biospies_data_v2, MEList_v2)
Biospies_Hub_proteins     <- chooseTopHubInEachModule(WGCNA_Biospies_data_v2, MEList_v2$validColors)
Biospies_High_MM_proteins <- get.High_ModuleMembership_proteins(WGCNA_Biospies_data_v2, MEs_v2$MEs, MEList_v2, cutoff = 0.8)
Biospies_High_MM_proteins_no_missing_data <- list()
for(module_color in names(Biospies_High_MM_proteins)) {
  Biospies_High_MM_proteins_no_missing_data[[module_color]]<- get.protein_without_missing_value(WGCNA_Biospies_data_v2, rownames(Biospies_High_MM_proteins[[module_color]]))
}

Biospies_module_color_matrix           <- as.matrix(MEs_v2$moduleColors)
rownames(Biospies_module_color_matrix) <- colnames(WGCNA_Biospies_data_v2)
Biospies_module_color                  <- get.Module_proteins(Biospies_module_color_matrix)

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


## run PCA to identify the sample
library(factoextra)
condition_biospies <- unlist(lapply(strsplit(rownames(WGCNA_Biospies_data_v2), "[.]"), `[[`, 1))

pdf(file = "PCA_Biopsies_NN_20191217.pdf")
get.pca(WGCNA_Biospies_data_v2, condition_biospies, name = "Total Biopsies")

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
Heart_Failure         <- Loaded_data(folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", file_name = "C0018801_disease_gda_summary.tsv")
Heart_Failure_Evidences <- Loaded_data(folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", file_name = "C0018801_disease_gda_evidences.tsv")
Heart_Failure_Evidences <- Heart_Failure_Evidences[Heart_Failure_Evidences$Gene %in% Heart_Failure$Gene, ] 

Heart_Failure_Type_of_Evidences <- list()
Heart_Failure_Type_of_Evidences[["Unique_proteins"]] <- unique(Heart_Failure$UniProt)
Heart_Failure_Type_of_Evidences[["Unique_genes"]] <- unique(Heart_Failure$Gene)
Heart_Failure_Type_of_Evidences[["Unique_genes_have_evidence"]] <- intersect(Heart_Failure_Evidences$Gene, Heart_Failure$Gene)
for (evidence_type in unique(Heart_Failure_Evidences$Association_Type)) {
  Heart_Failure_Type_of_Evidences[[evidence_type]] <- unique(Heart_Failure_Evidences$Gene[Heart_Failure_Evidences$Association_Type == evidence_type])
}

Proteins_Invitro_ovelap_Heart_Failure  <- Select_if_in_specific_list(t(WGCNA_Invitro_data), selected_list = unique(Heart_Failure$UniProt))
Proteins_Invitro_ovelap_Heart_Failure_no_mising_value <- get.protein_without_missing_value(WGCNA_Invitro_data, rownames(Proteins_Invitro_ovelap_Heart_Failure))
Module_Invitro_ovelap_Heart_Failure <- get.Module_Proteins_overlap_DisGeNET(Invitro_module_color, Heart_Failure, WGCNA_Invitro_data)

Proteins_Biospies_ovelap_Heart_Failure <- Select_if_in_specific_list(t(WGCNA_Biospies_data_v2), selected_list = unique(Heart_Failure$UniProt))
Proteins_Biospies_ovelap_Heart_Failure_no_mising_value <- get.protein_without_missing_value(WGCNA_Biospies_data_v2, rownames(Proteins_Biospies_ovelap_Heart_Failure))
Module_Biospies_ovelap_Heart_Failure <- get.Module_Proteins_overlap_DisGeNET(Biospies_module_color, Heart_Failure, WGCNA_Biospies_data_v2)

# module high MM vs. Heart failure
Invitro_High_MM_Heart_Failure<-list()
for (module_inVitro in names(Invitro_High_MM_proteins)) {
  Invitro_High_MM_Heart_Failure[[module_inVitro]]<- intersect(rownames(Invitro_High_MM_proteins[[module_inVitro]]), 
                                                              Module_Invitro_ovelap_Heart_Failure$Protein_overlap_DisGeNet[[module_inVitro]])
}

Biopsies_High_MM_Heart_Failure<-list()
for (module_biospies in names(Biospies_High_MM_proteins)) {
  Biopsies_High_MM_Heart_Failure[[module_biospies]]<- intersect(rownames(Biospies_High_MM_proteins[[module_biospies]]), 
                                                                Module_Biospies_ovelap_Heart_Failure$Protein_overlap_DisGeNet[[module_biospies]])
}

Protein_High_MM_opverlap<-list()
for (module_inVitro in names(Invitro_High_MM_proteins)) {
  for (module_biospies in names(Biospies_High_MM_proteins)) {
    names_module <- paste0(module_inVitro, "_", module_biospies)
    Protein_High_MM_opverlap[[names_module]] <- intersect(rownames(Invitro_High_MM_proteins[[module_inVitro]]), 
                                                          rownames(Biospies_High_MM_proteins[[module_biospies]]))
  }
}
unlist(Protein_High_MM_opverlap)

intersect(unlist(Protein_High_MM_opverlap), unique(Heart_Failure$UniProt))

Print_protein_expression_log2FC_in_pdf(WGCNA_Invitro_data, unlist(Protein_High_MM_opverlap), Protein_metadata, 
                                       save_folder = "I:/HeCaToS/Paper/Paper1_Protein_Anthacycline/Data_and_Code", 
                                       file_name = "WGCNA_Overlap_HighMM_Protein_invitro_NN_20191219.pdf")
Print_protein_expression_log2FC_biopsies(WGCNA_Biospies_data_v2, unlist(Protein_High_MM_opverlap)) 

## Venn diagram
library(VennDiagram)
venn.diagram(
  x = list(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2), unique(Heart_Failure$UniProt)),
  fill = c("red", "orange", "blue"), sub.cex = 2,
  category.names = c("Invitro" , "Biopsies", "DisGeNet_Heart_Failure"),
  filename = 'Venn_diagramm_total_protein_NN_20191218.png', lwd = 2, apha = 0.5,  output=FALSE
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
                             filename = paste0(module_inVitro, "_", module_biospies, "_DisGeNet", "_NN_20191218.png"),
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
write.csv(output, "protein_data_NN_20191223.csv")


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
write.csv(output, "protein_data_v2_NN_20191223.csv")


## Check the overlap in proteomics data between in vitro vs. biospies:
library(VennDiagram)
venn.diagram(
  x = list(rownames(Invitro_data_log), rownames(Biospies_data_log)),
  category.names = c("ANT_Proteins" , "Bio_Proteins"),
  filename = '#14_venn_diagramm.png', lwd = 2,
  output=TRUE
)



library(VennDiagram)
venn.diagram(
  x = list(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2)),
  category.names = c("ANT_Proteins" , "Bio_Proteins"),
  filename = '#15_venn_diagramm.png', lwd = 2,
  output=TRUE
)



library(VennDiagram)
venn.diagram(
  x = list(colnames(WGCNA_Invitro_data), colnames(WGCNA_Biospies_data_v2), unique(Sum_heart_disease_from_Disgenet$UniProt)),
  category.names = c("ANT_Proteins" , "Bio_Proteins", "DisGeNet_Proteins"),
  filename = '#16_venn_diagramm.png', lwd = 3,
  output=TRUE
)

library(VennDiagram)
venn.diagram(
  x = list(Invitro_module_color$turquoise, colnames(WGCNA_Biospies_data_v2), unique(Sum_heart_disease_from_Disgenet$UniProt)),
  category.names = c("ANT_Proteins_Turquoise" , "Bio_Proteins", "DisGeNet_Proteins"),
  filename = '#17_venn_diagramm.png', lwd = 3,
  output=TRUE
)

a<-intersect(Invitro_module_color$turquoise, colnames(WGCNA_Biospies_data_v2))
data <- Select_if_in_specific_list(t(WGCNA_Biospies_data_v2), selected_list = a)
heatmap(data)
data <- Select_if_in_specific_list(t(WGCNA_Invitro_data), selected_list = a)
heatmap(data)
write.table(a, row.names = FALSE, "test.txt")

a2<-intersect(a, unique(Sum_heart_disease_from_Disgenet$UniProt))
write.table(a2, row.names = FALSE, "test2.txt")


library(VennDiagram)
venn.diagram(
  x = list(Invitro_module_color$turquoise, Biospies_module_color$turquoise,
           unique(Sum_heart_disease_from_Disgenet$UniProt)),
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
                          unique(Sum_heart_disease_from_Disgenet$UniProt)),
                 category.names = c(paste0("ANT_Proteins_", module_inVitro),
                                    paste0("Bio_Proteins_", module_biospies),
                                    "DisGeNet_Proteins"),
                 filename = paste0(module_inVitro, "_", module_biospies, "_DisGeNet", ".png"),
                 output=TRUE
    )
  }
}


venn.diagram(x = list(Invitro_module_color$turquoise, rownames(Biospies_data_log), 
                      unique(Heart_Failure$UniProt), unique(Sum_heart_disease_from_Disgenet$UniProt)),
             category.names = c("ANT_Proteins_Turquoise", "Bio_Proteins", 
                                "HeartFailure_Proteins", "DisGeNET_Protein"),
             filename = "ANT_Proteins_Turquoise_HeartFailure.png",  output=TRUE)


### Function: ================================================
Loaded_data <- function(folder, file_name) {
  setwd(folder)
  output<-read.table(file_name, header = TRUE, sep = "\t") 
  return(output)
}

get.rowname <- function(data) {
  rownames(data) <- data[, 1]
  data <- data[, -1]
  return(data)
}

Change_a_Part_in_Names <- function(list_input, select_chatacter, wanted_character) {
  output <- gsub(select_chatacter, wanted_character, list_input)
  return(output)
}

Identify_if_have_speficic_symbol <- function(data, symbol) {
  output <-data[which(grepl(symbol, data))]
  return(output)
}

Remove_if_have_specific_symbol <- function (data, remove_list) {
  output <- data[!rownames(data) %in% remove_list,]
  return(output)
}

Round_number_in_data <- function(data, singificant_number) {
  output <- signif(data, singificant_number)
  return(output)
}


Select_Part_Names <- function(list, separator, selected_position) {
  # eg. for separator = "|" --> "[|]"
  output<-c()
  for (i in 1: length(list)) output[i] <- strsplit(list[i], paste0("[", separator, "]"))[[1]][selected_position]
  return(output)
}


Filter_data <- function(data) {
  gsg <- goodSamplesGenes(data, verbose = 3);				# check data
  if (gsg$allOK) {
    print(paste0("The data pass the WGCNA requirment so we do not use the WGCNA fillter. Number of protein = ", ncol(data), "; Number of sample = ",  nrow(data)))
    
  }
  else {
    print(paste0("The data did not pass the WGCNA requirment so we used WGCNA fillter to remove low quality genes. Before filtering: Number of protein = ", ncol(data), "; Number of sample = ",  nrow(data)))
    data = data[gsg$goodSamples, gsg$goodGenes]
    print(paste0("The data did not pass the WGCNA requirment so we used WGCNA fillter to remove low quality genes. Afther filtering: Number of protein = ", ncol(data), "; Number of sample = ",  nrow(data)))
  }
  output <- data
  return(output)
}


Made_sample_tree <- function(data) {
  sampleTree = hclust(dist(data), method = "average") 
  plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
}

Made_sample_tree2 <- function(data, categorial_samples) {
  library(dendextend)
  
  dend <- as.dendrogram(hclust(dist(data), method = "average"))
  labels_colors(dend) <- categorial_samples
  plot(dend, main = "Sample clustering")
}


Pick_WGCNA_power <- function(data) {
  powers <- c(c(1:10), seq(from = 12, to=20, by=2))	
  sft <- pickSoftThreshold(data, powerVector = powers, verbose = 5) 
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  par(mfrow=c(2,2))
  plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"))
  text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], labels = powers, cex = 0.9, col = "red")
  abline(h=0.90,col="red") 							# this line corresponds to using an R^2 cut-off of h
  
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = paste("Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9 ,col = "red")
}


Build_WGCNA_geneTree <- function(data, network_type, picked_power) {
  adjacency <- adjacency(datExpr = data, type = network_type, power = picked_power)
  dissTOM = TOMdist(adjMat = adjacency, TOMType = network_type)
  geneTree = hclust(as.dist(dissTOM), method = "average")
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 3, pamRespectsDendro = FALSE, minClusterSize = 30)  # Module identification using dynamic tree cut
  dynamicColors = labels2colors(dynamicMods)			# Convert numeric lables into colors
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
  return(list(geneTree = geneTree, dynamicColors = dynamicColors))
}


Module_Eigengen_Tree<- function(data, Module_color) {
  MEList = moduleEigengenes(data, colors = Module_color)	# Calculate eigengenes
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs, use = "p") 	 # Calculate dissimilarity of module eigengenes,  use "pairwise complete observations", to omitted NAs
  MEDiss[is.na(MEDiss)]<-0           # convert the NaN value to 0
  METree = hclust(as.dist(MEDiss), method = "average")		# Cluster module eigengenes
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  return(MEs)
}


get.MEtree_plot <- function(MEs, MEDissThres) {
  MEDiss = 1-cor(MEs)
  METree = hclust(as.dist(MEDiss), method = "average");
  sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
}


get.MEs_merge <- function(datExpr, dynamicColors, MEDissThres, geneTree) {
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors = merge$colors;
  mergedMEs = merge$newMEs;
  
  sizeGrWindow(12, 9)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  moduleColors = mergedColors # Rename to moduleColors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs
  return(list("MEs"= mergedMEs, "moduleColors" = moduleColors))
}

get.Hub_proteins <- function(dataExpr, MElist) {
  TopHub_proteins<- chooseTopHubInEachModule(dataExpr, MEList$validColors)
  return(TopHub_proteins)
}

get.High_ModuleMembership_proteins <- function(dataExpr, MEs, MElist, cutoff) {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
  HighMM_protein <- list()
  for (i in substring(names(MEs), 3)) {
    MM_measurement <-geneModuleMembership[which(MEList$validColors == i),]
    HighMM_protein[[i]] <- MM_measurement[which(MM_measurement[, paste0("ME", i)] >= cutoff),]
  }
  return(HighMM_protein)
}
get.Module_proteins <- function(module_color) {
  Module_proteins <- list()
  for(color in unique(as.vector(module_color))) {
    Module_proteins[[color]] <- rownames(module_color)[module_color == color]
  }
  return(Module_proteins)
}

Select_if_in_specific_list <- function (data, selected_list) {
  output <- data[rownames(data) %in% selected_list,]
  return(output)
}

get.pca <- function(data, condition, name) {
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))
}

get.pca_with_time <- function(data, condition, name) {
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  
  library(ggrepel)
  metadata <- as.data.frame(get.metadata(rownames(data)))
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))  + geom_text_repel(label = metadata$Time, size = 3)
}

get.pca_with_selected_time <- function(data, condition, name) {
  data[is.na(data)] <- 0
  res.pca <- prcomp(data, scale = TRUE)
  
  library(ggrepel)
  metadata <- as.data.frame(get.metadata(rownames(data)))
  metadata$Time[which(metadata$Time == "002"| metadata$Time == "008"| metadata$Time == "024"| metadata$Time == "072")] <- ""
  fviz_pca_ind(res.pca, geom="point",  pointsize = 2, 
               habillage=condition, addEllipses=TRUE, ellipse.level=0.95,  
               title = paste0("PCA - ", name))  + geom_text_repel(label = metadata$Time, size = 3)
}

get.metadata <- function(name) {
  output <- matrix(unlist(strsplit(name, "_")), ncol=4, byrow=TRUE)
  colnames(output) <- c("Drug", "Dose", "Time", "Replicate")
  return(output)
}

get.protein_without_missing_value <- function(Protein_data, select_list) {
  Data_no_missing_value <- na.omit(t(Protein_data))
  output <- Data_no_missing_value[rownames(Data_no_missing_value) %in% select_list,]
  return(output)
}


make.sample_names <- function(sample_names) {
  convert_names <- matrix(unlist(strsplit(sample_names, "_")), ncol = 4, byrow = TRUE)
  colnames(convert_names) <- c("Drug", "Dose", "Time", "Replicate")
  
  convert_names[which(convert_names[, "Drug"] == "Fluct"), 1]       <- "Con"
  convert_names[which(convert_names[, "Dose"] == "Therapeutic"), 2] <- "The"
  convert_names[which(convert_names[, "Dose"] == "Toxic"), 2]       <- "Tox"
  convert_names[which(convert_names[, "Time"] == "T2"), 3]          <- "002"
  convert_names[which(convert_names[, "Time"] == "T8"), 3]          <- "004"
  convert_names[which(convert_names[, "Time"] == "T24"), 3]         <- "024"
  convert_names[which(convert_names[, "Time"] == "T72"), 3]         <- "072"
  convert_names[which(convert_names[, "Time"] == "T168"), 3]        <- "168"
  convert_names[which(convert_names[, "Time"] == "T240"), 3]        <- "240"
  convert_names[which(convert_names[, "Time"] == "T336"), 3]        <- "336"
  
  output <- apply(convert_names, 1, paste , collapse = "_" )
  return(output)
}


Print_ME_in_pdf <- function(MEs, metadata, save_folder, file_name) {
  if("plyr" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: plyr")
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: ggplot2")
  if("gridExtra" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: gridExtra")
  
  library("plyr")						
  library("ggplot2")		
  library("gridExtra")
  
  setwd(save_folder)
  pdf(file_name, onefile = TRUE, width = 15) # set upt the pdf file
  
  plist<- list()
  for (i in colnames(MEs))	{
    ME_tem    <- as.data.frame(cbind(metadata, setNames(MEs[i], "MEs")))
    Mean_tem  <- ddply(ME_tem, ~Dose+Time, summarise, Mean = mean(MEs, na.rm = TRUE), sd = sd(MEs, na.rm = TRUE))
    plist[[i]]<- ggplot(Mean_tem, aes(x = Time, y = Mean, ymin = -0.4, ymax = 0.4, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
      xlab("Time (hours)") + ylab("Eigengens value") + ggtitle(i) +  theme_bw() 
  }
  ml <- marrangeGrob(plist, nrow=3, ncol=2)
  print(ml)
  dev.off() # Close the pdf file
}

get.biopsies_type <- function(Biospies_info) {
  Samples <- c()
  Samples$Patient[Biospies_info$Control...Cardiotoxicity == "Late onset cardiotoxicity"] <- "Cardiotoxicity"
  Samples$Patient[Biospies_info$Control...Cardiotoxicity == "Control patients"]          <- "Control"
  
  Samples$Drug[Biospies_info$Control...Cardiotoxicity == "Control patients"]           <- "" # control samples
  Samples$Drug[Biospies_info$Chemotherapeutic.agents == ""]            <- "_no_data"
  Samples$Drug[grep("rubicin", Biospies_info$Chemotherapeutic.agents)] <- "_with_ANT"
  Samples$Drug[grep("cycline", Biospies_info$Chemotherapeutic.agents)] <- "_with_ANT"
  Samples$Drug[is.na(Samples$Drug)]                                    <- "_without_ANT"
  
  summary_sample           <- matrix(unlist(Samples), ncol =2, byrow = FALSE)
  output <- apply(summary_sample, 1, paste , collapse = "" )
  return(output)
}

Print_protein_expression_in_pdf <- function(expression_data, selected_list, metadata, save_folder, file_name) {
  if("plyr" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: plyr")
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: ggplot2")
  if("gridExtra" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: gridExtra")
  
  library("plyr")						
  library("ggplot2")		
  library("gridExtra")
  
  
  Protein_expression <- as.data.frame(t(Select_if_in_specific_list(t(expression_data), selected_list)))
  
  setwd(save_folder)
  pdf(file_name, onefile = TRUE, width = 15) # set upt the pdf file
  
  plist<- list()
  for (i in colnames(Protein_expression))	{
    Protein_expression_tem    <- as.data.frame(cbind(metadata, setNames(Protein_expression[i], "MEs")))
    Mean_tem  <- ddply(Protein_expression_tem, ~Dose+Time, summarise, Mean = mean(MEs, na.rm = TRUE), sd = sd(MEs, na.rm = TRUE))
    plist[[i]]<- ggplot(Mean_tem, aes(x = Time, y = Mean, ymin = -0.4, ymax = 0.4, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
      xlab("Time (hours)") + ylab("Eigengens value") + ggtitle(i) +  theme_bw() 
  }
  ml <- marrangeGrob(plist, nrow=3, ncol=2)
  print(ml)
  dev.off() # Close the pdf file
}

Print_protein_expression_log2FC_in_pdf <- function(expression_data, selected_list, metadata, save_folder, file_name) {
  if("plyr" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: plyr")
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: ggplot2")
  if("gridExtra" %in% rownames(installed.packages()) == FALSE) print("Error!!! Install R package: gridExtra")
  
  library("plyr")						
  library("ggplot2")		
  library("gridExtra")
  
  
  Protein_expression <- as.data.frame(t(Select_if_in_specific_list(t(expression_data), selected_list)))
  
  setwd(save_folder)
  pdf(file_name, onefile = TRUE, width = 15) # set upt the pdf file
  
  plist_mean   <- list()
  plist_Log2FC <- list()
  for (i in colnames(Protein_expression))	{
    Protein_expression_tem    <- as.data.frame(cbind(metadata, setNames(Protein_expression[i], "MEs")))
    Mean_tem  <- ddply(Protein_expression_tem, ~Dose+Time, summarise, Mean = mean(MEs, na.rm = TRUE), sd = sd(MEs, na.rm = TRUE))
    
    Log2FC_tem <- Mean_tem %>% dplyr::group_by(Dose) %>% mutate(Log2FC = Mean -  Mean_tem[Mean_tem$Dose == 'Con_DMSO', 'Mean'])
    Log2FC_tem <- Log2FC_tem[-which(Log2FC_tem$Dose == "Con_DMSO"),]
    
    plist_mean[[i]]<- ggplot(Mean_tem, aes(x = Time, y = Mean, ymin = -0.4, ymax = 0.4, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() + geom_errorbar(aes(ymin = Mean-sd, ymax= Mean+sd), width=.2,position=position_dodge(0.05)) +
      xlab("Time (hours)") + ylab("Expression value") + ggtitle(i) +  theme_bw() 
    
    plist_Log2FC[[i]]<- ggplot(Log2FC_tem, aes(x = Time, y = Log2FC, ymin = -0.4, ymax = 0.4, colour = Dose, group = Dose)) + 
      geom_line() + geom_point() +
      xlab("Time (hours)") + ylab("Log2FC value") + ggtitle(i) +  theme_bw()
  }
  ml <- marrangeGrob(cbind(plist_mean, plist_Log2FC), nrow=3, ncol=2)
  print(ml)
  dev.off() # Close the pdf file
}


get.Module_Proteins_overlap_DisGeNET <- function(module_color, DisGeNET, expression_data) {
  Module_Proteins_overlap_DisGeNET <- list()
  Module_Proteins_overlap_DisGeNET_no_mising_value <- list()
  for (module in names(module_color)) {
    Module_Proteins_overlap_DisGeNET[[module]] <- intersect(module_color[[module]], unique(DisGeNET$UniProt))
    Module_Proteins_overlap_DisGeNET_no_mising_value[[module]] <- get.protein_without_missing_value(expression_data, Module_Proteins_overlap_DisGeNET[[module]])
  }
  output <-list("Protein_overlap_DisGeNet" = Module_Proteins_overlap_DisGeNET, "Protein_overlap_DisGeNet_no_missing_value" = Module_Proteins_overlap_DisGeNET_no_mising_value)
  return(output)
}

Print_protein_expression_log2FC_biopsies <- function(expression_data, selected_list) {
  library("plyr")	
  
  Protein_expression <- as.data.frame(t(Select_if_in_specific_list(t(expression_data), selected_list)))
  Protein_expression$condition_biospies <- unlist(lapply(strsplit(rownames(Protein_expression), "[.]"), `[[`, 1))
  
  Mean_expression <- Protein_expression %>% dplyr::group_by(condition_biospies) %>% dplyr::summarise_all(.funs = c(mean = "mean", Sd="sd"), na.rm = TRUE)
  Mean_expression <- get.rowname(as.matrix(Mean_expression))
  class(Mean_expression) <- "numeric"
  Mean_expression_v1 <- Mean_expression[, grep("mean", colnames(Mean_expression))]
  
  op <- par(mar=c(14,4,4,2))
  barplot(Mean_expression_v1, col=c("red", "green", "blue"),  las = 2, beside = TRUE)
  barplot(Mean_expression_v1, col=c("red", "green", "blue"), legend = rownames(Mean_expression_v1),  las = 2, beside = TRUE)
  rm(op)
  
  
  Mean_Log2FC<- sweep(Mean_expression_v1[c("Cardiotoxicity_with_ANT", "Cardiotoxicity_without_ANT"),  ], 2, Mean_expression_v1["Control", ])
  op <- par(mar=c(14,4,4,2))
  barplot(Mean_Log2FC, col=c("red", "green"),  las = 2, beside = TRUE)
  barplot(Mean_Log2FC, col=c("red", "green"), legend = rownames(Mean_Log2FC),  las = 2, beside = TRUE)
  rm(op)
}

get.High_MM_proteins_list <- function(High_MM_proteins) {
  output <- c()
  for(module in names(High_MM_proteins)) {
    output <- c(output, rownames(High_MM_proteins[[module]]))
  }
  return(output)
}


