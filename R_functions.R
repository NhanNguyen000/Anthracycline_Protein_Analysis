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