library(maftools)
library(SAMBAR)
library(survival)
library(GSVAdata)
library(GSVA)
library(ComplexHeatmap)
library(readxl)
library(dplyr)
library(UpSetR)
library(ggsci)
library(table1)
library(survminer)
library(gtools)

source("function.R")

#
Sample2ID_109 <- read.table("data/109_Sample2ID.txt",header = T,sep = "\t",check.names = F)
rownames(Sample2ID_109) <- Sample2ID_109$id
clinical_109 <- read.table("data/109_clinical.xls",header = T,sep="\t",check.names = F)
clinical_109[,"Patient ID"] <- Sample2ID_109[clinical_109$`Patient ID`,"sample"]
clinical_109$`Patient ID` <- as.character(clinical_109$`Patient ID`)

#
exon_mutation_types <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
                         "Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Splice_Site")

promoter_mutation_types <- c("5'Flank","5'UTR")

LUSC_maf <- read.maf("data/LUSC.somatic.maf",vc_nonSyn = c(exon_mutation_types,promoter_mutation_types))

mutation_109 <- read.table("data/LUSC.somatic.maf",header = TRUE,sep="\t",quote = "")
mutation_109$Tumor_Sample_Barcode <- as.character(mutation_109$Tumor_Sample_Barcode)

mutation_109[,"Patient ID"] <- mutation_109$Tumor_Sample_Barcode
mutation_109$Cancer_Type <- "LUSC"

mutation_109 <- mutation_109[,c("Hugo_Symbol","Variant_Classification","Tumor_Sample_Barcode","Patient ID")]


#
list_109_exon <- add_mutation_data_to_clinical_exon(clinical_df = clinical_109,mutation_df = mutation_109)
list_109 <- add_mutation_data_to_clinical(clinical_df = clinical_109,mutation_df = mutation_109)

#

sigpw_genes <- as.data.frame(data.table::fread(system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")))



#
write_gmt(sigpw_genes, "oncogenic_sig_patwhays.gmt")

signatureset_sigpw <- "oncogenic_sig_patwhays.gmt"


library(gson)
signatureset_c2 <- "data/c2.all.v2023.1.Hs.symbols.gmt"
signatureset_c2_list <- read.gmt(signatureset_c2)

#Ten pathways for exon mutations

exon_gene_sigpw <- intersect(unique(sigpw_genes$Gene),unique(unique(list_109_exon[[2]][list_109_exon[[2]][,3]>6,1])))
promoter_gene_sigpw <- intersect(unique(sigpw_genes$Gene),unique(unique(list_109[[2]][list_109[[2]][,3]>6,1])))

exon_gene_mat_sigpw <- list_109_exon[[1]][,exon_gene_sigpw]

exon_pathway_mat_sigpw <- analyze_signature(list_109_exon[[1]],exon_gene_sigpw,signatureset_sigpw,0)

promoter_pathway_mat_sigpw <- analyze_signature(list_109[[1]],promoter_gene_sigpw,signatureset_sigpw,0)



#Ten pathways for promoter mutations

promoter_gene_mat_sigpw <- list_109[[1]][,promoter_gene_sigpw]


#Calculate the multivariate prognosis of exon mutant genes in cohort 109


multi_data_exon <- na.omit(list_109_exon[[1]])
multi_data_exon[multi_data_exon$Stage%in%c("I"),"Stage"] <- 1
multi_data_exon[multi_data_exon$Stage%in%c("II"),"Stage"] <- 2
multi_data_exon[multi_data_exon$Stage%in%c("III"),"Stage"] <- 3
multi_data_exon[multi_data_exon$Stage%in%c("IV"),"Stage"] <- 4
multi_data_exon$Stage <- as.numeric(multi_data_exon$Stage)

multi_data_exon[multi_data_exon$Gender=="Male","Gender"] <- 1
multi_data_exon[multi_data_exon$Gender=="Female","Gender"] <- 0
multi_data_exon$Gender <- as.numeric(multi_data_exon$Gender)

colnames(multi_data_exon) <- make.names(colnames(multi_data_exon))


multicox_OS_pathway_exon_genes <- run_multicox_OS(make.names(exon_gene_sigpw),multi_data_exon,cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))

multicox_RFS_pathway_exon_genes <- run_multicox_RFS(make.names(exon_gene_sigpw),multi_data_exon,cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))



#Calculate the multivariate prognosis of promoter mutant genes

multi_data_promoter <- na.omit(list_109[[1]])
multi_data_promoter[multi_data_promoter$Stage%in%c("I"),"Stage"] <- 1
multi_data_promoter[multi_data_promoter$Stage%in%c("II"),"Stage"] <- 2
multi_data_promoter[multi_data_promoter$Stage%in%c("III"),"Stage"] <- 3
multi_data_promoter[multi_data_promoter$Stage%in%c("IV"),"Stage"] <- 4
multi_data_promoter$Stage <- as.numeric(multi_data_promoter$Stage)

multi_data_promoter[multi_data_promoter$Gender=="Male","Gender"]<- 1
multi_data_promoter[multi_data_promoter$Gender=="Female","Gender"]<- 0
multi_data_promoter$Gender <- as.numeric(multi_data_promoter$Gender)

colnames(multi_data_promoter) <- make.names(colnames(multi_data_promoter))


multicox_OS_pathway_promoter_genes <- run_multicox_OS(make.names(promoter_gene_sigpw),multi_data_promoter,cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))
multicox_RFS_pathway_promoter_genes <- run_multicox_RFS(make.names(promoter_gene_sigpw),multi_data_promoter,cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))


#Summary of the multivariate prognosis of the 109 cohort mutant genes

multi_data_promoter_length <- colSums(multi_data_promoter[,make.names(promoter_gene_sigpw)])
multi_data_exon_length <- colSums(multi_data_exon[,make.names(exon_gene_sigpw)])

multicox_109_pathway_exon_genes <-merge(multicox_OS_pathway_exon_genes,multicox_RFS_pathway_exon_genes,by="Gene")

multicox_109_pathway_promoter_genes <- merge(multicox_OS_pathway_promoter_genes,multicox_RFS_pathway_promoter_genes,by="Gene")

multicox_109_pathway_exon_genes$cohort <- "109"
multicox_109_pathway_promoter_genes$cohort <- "109"

multicox_109_pathway_exon_genes$Mutation_Number <- as.vector(unlist(multi_data_exon_length[multicox_109_pathway_exon_genes$Gene]))
multicox_109_pathway_promoter_genes$Mutation_Number <- as.vector(unlist(multi_data_promoter_length[multicox_109_pathway_promoter_genes$Gene]))



#Prognosis of exon mutation pathways in cohort 109

exon_pathway_mat_sigpw_clinical <- as.data.frame(t(exon_pathway_mat_sigpw)) %>%  tibble::rownames_to_column(var = "Patient ID")
exon_pathway_mat_sigpw_clinical_merge <- merge(clinical_109,exon_pathway_mat_sigpw_clinical,by="Patient ID",all="TRUE")

exon_pathway_mat_sigpw_clinical_merge <- na.omit(exon_pathway_mat_sigpw_clinical_merge)
exon_pathway_mat_sigpw_clinical_merge[exon_pathway_mat_sigpw_clinical_merge$Stage%in%c("I"),"Stage"] <- 1
exon_pathway_mat_sigpw_clinical_merge[exon_pathway_mat_sigpw_clinical_merge$Stage%in%c("II"),"Stage"] <- 2
exon_pathway_mat_sigpw_clinical_merge[exon_pathway_mat_sigpw_clinical_merge$Stage%in%c("III"),"Stage"] <- 3
exon_pathway_mat_sigpw_clinical_merge[exon_pathway_mat_sigpw_clinical_merge$Stage%in%c("IV"),"Stage"] <- 4
exon_pathway_mat_sigpw_clinical_merge$Stage <- as.numeric(exon_pathway_mat_sigpw_clinical_merge$Stage)

exon_pathway_mat_sigpw_clinical_merge[exon_pathway_mat_sigpw_clinical_merge$Gender=="Male","Gender"]<- 1
exon_pathway_mat_sigpw_clinical_merge[exon_pathway_mat_sigpw_clinical_merge$Gender=="Female","Gender"]<- 0
exon_pathway_mat_sigpw_clinical_merge$Gender <- as.numeric(exon_pathway_mat_sigpw_clinical_merge$Gender)

colnames(exon_pathway_mat_sigpw_clinical_merge) <- make.names(colnames(exon_pathway_mat_sigpw_clinical_merge))
rownames(exon_pathway_mat_sigpw_clinical_merge) <- exon_pathway_mat_sigpw_clinical_merge$Patient.ID

exon_pathway_mat_sigpw_clinical_merge_length <- colSums(exon_pathway_mat_sigpw_clinical_merge[,12:ncol(exon_pathway_mat_sigpw_clinical_merge)])

exon_pathway_multicox_os_survial <- run_multicox_OS(make.names(rownames(exon_pathway_mat_sigpw)),exon_pathway_mat_sigpw_clinical_merge,
                                                    cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))

exon_pathway_multicox_rfs_survial <- run_multicox_RFS(make.names(rownames(exon_pathway_mat_sigpw)),exon_pathway_mat_sigpw_clinical_merge,
                                                      cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))

#The prognosis of the 109 cohort promoter mutation pathway

promoter_pathway_mat_sigpw_clinical <- as.data.frame(t(promoter_pathway_mat_sigpw)) %>%  tibble::rownames_to_column(var = "Patient ID")
promoter_pathway_mat_sigpw_clinical_merge <- merge(clinical_109,promoter_pathway_mat_sigpw_clinical,by="Patient ID",all="TRUE")

promoter_pathway_mat_sigpw_clinical_merge <- na.omit(promoter_pathway_mat_sigpw_clinical_merge)
promoter_pathway_mat_sigpw_clinical_merge[promoter_pathway_mat_sigpw_clinical_merge$Stage%in%c("I"),"Stage"] <- 1
promoter_pathway_mat_sigpw_clinical_merge[promoter_pathway_mat_sigpw_clinical_merge$Stage%in%c("II"),"Stage"] <- 2
promoter_pathway_mat_sigpw_clinical_merge[promoter_pathway_mat_sigpw_clinical_merge$Stage%in%c("III"),"Stage"] <- 3
promoter_pathway_mat_sigpw_clinical_merge[promoter_pathway_mat_sigpw_clinical_merge$Stage%in%c("IV"),"Stage"] <- 4
promoter_pathway_mat_sigpw_clinical_merge$Stage <- as.numeric(promoter_pathway_mat_sigpw_clinical_merge$Stage)

promoter_pathway_mat_sigpw_clinical_merge[promoter_pathway_mat_sigpw_clinical_merge$Gender=="Male","Gender"]<- 1
promoter_pathway_mat_sigpw_clinical_merge[promoter_pathway_mat_sigpw_clinical_merge$Gender=="Female","Gender"]<- 0
promoter_pathway_mat_sigpw_clinical_merge$Gender <- as.numeric(promoter_pathway_mat_sigpw_clinical_merge$Gender)

colnames(promoter_pathway_mat_sigpw_clinical_merge) <- make.names(colnames(promoter_pathway_mat_sigpw_clinical_merge))
rownames(promoter_pathway_mat_sigpw_clinical_merge) <- promoter_pathway_mat_sigpw_clinical_merge$Patient.ID
promoter_pathway_mat_sigpw_clinical_merge_length <- colSums(promoter_pathway_mat_sigpw_clinical_merge[,12:ncol(promoter_pathway_mat_sigpw_clinical_merge)])

promoter_pathway_multicox_os_survial <- run_multicox_OS(make.names(rownames(promoter_pathway_mat_sigpw)),promoter_pathway_mat_sigpw_clinical_merge,
                                                        cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))

promoter_pathway_multicox_rfs_survial <- run_multicox_RFS(make.names(rownames(promoter_pathway_mat_sigpw)),promoter_pathway_mat_sigpw_clinical_merge,
                                                          cli_factor = c("Age", "Stage", "Smoking.0.non.smoker.1.smoker."))


#Summary of 109 Cohort Pathway Analysis


exon_pathway_multicox_merge_survial <- merge(exon_pathway_multicox_os_survial,exon_pathway_multicox_rfs_survial,by="Gene")
exon_pathway_multicox_merge_survial$cohort <- "109"
exon_pathway_multicox_merge_survial$Mutation_Number <- as.vector(unlist(exon_pathway_mat_sigpw_clinical_merge_length[exon_pathway_multicox_merge_survial$Gene]))


promoter_pathway_multicox_merge_survial <- merge(promoter_pathway_multicox_os_survial,promoter_pathway_multicox_rfs_survial,by="Gene")
promoter_pathway_multicox_merge_survial$cohort <- "109"
promoter_pathway_multicox_merge_survial$Mutation_Number <- as.vector(unlist(promoter_pathway_mat_sigpw_clinical_merge_length[promoter_pathway_multicox_merge_survial$Gene]))




#TCGA data analysis


library(maftools)
#Exon mutations of the TCGA gene and prognosis

load("data/TCGA_tidy_Clinical.RData")

TCGA_Clinical.tidy.lusc <- as.data.frame(TCGA_Clinical.tidy[TCGA_Clinical.tidy$Project=="LUSC",])[,c("Tumor_Sample_Barcode","Age","Smoking_history","Tumor_stage","OS","OS.time","RFS","RFS.time")]
TCGA_Clinical.tidy.lusc <- TCGA_Clinical.tidy.lusc[substr(TCGA_Clinical.tidy.lusc$Tumor_Sample_Barcode,14,15) %in% "01",]
TCGA_Clinical.tidy.lusc$patients <- make.names(substr(TCGA_Clinical.tidy.lusc$Tumor_Sample_Barcode,1,12))
TCGA_Clinical.tidy.lusc <- TCGA_Clinical.tidy.lusc[,c("patients","Age","Smoking_history","Tumor_stage","OS","OS.time","RFS","RFS.time")]
TCGA_Clinical.tidy.lusc$OS.time <- TCGA_Clinical.tidy.lusc$OS.time/30
TCGA_Clinical.tidy.lusc$RFS.time <- TCGA_Clinical.tidy.lusc$RFS.time/30

colnames(TCGA_Clinical.tidy.lusc) <- c("patients","Age","Smoking_history","Stage","Survival.status","OS.months.","Recurrence.status","RFS.months.")
TCGA_Clinical.tidy.lusc <- unique(TCGA_Clinical.tidy.lusc)
rownames(TCGA_Clinical.tidy.lusc) <- TCGA_Clinical.tidy.lusc$patients
TCGA_Clinical.tidy.lusc$Stage <- as.vector(TCGA_Clinical.tidy.lusc$Stage)
TCGA_Clinical.tidy.lusc[TCGA_Clinical.tidy.lusc$Stage%in%"I","Stage"] <- 1
TCGA_Clinical.tidy.lusc[TCGA_Clinical.tidy.lusc$Stage%in%"II","Stage"] <- 2
TCGA_Clinical.tidy.lusc[TCGA_Clinical.tidy.lusc$Stage%in%"III","Stage"] <- 3
TCGA_Clinical.tidy.lusc[TCGA_Clinical.tidy.lusc$Stage%in%"IV","Stage"] <- 4
TCGA_Clinical.tidy.lusc[TCGA_Clinical.tidy.lusc$Stage%in%"X","Stage"] <- NA
TCGA_Clinical.tidy.lusc$Stage <- as.numeric(TCGA_Clinical.tidy.lusc$Stage)


#Read the exon mutation information
TCGA_LUSC_maf <- read.maf("data/TCGA_maf/TCGA.LUSC.mutect.95258183-63ea-4c97-ae29-1bae9ed06334.DR-10.0.somatic.maf.gz")

TCGA_LUSC_maf_data <- TCGA_LUSC_maf@data

TCGA_LUSC_maf_data_select <- TCGA_LUSC_maf_data[,c("Hugo_Symbol","Tumor_Sample_Barcode")]
TCGA_LUSC_maf_data_select$Tumor_Sample_Barcode <- as.vector(TCGA_LUSC_maf_data_select$Tumor_Sample_Barcode)
TCGA_LUSC_maf_data_select <- unique(TCGA_LUSC_maf_data_select)

#TCGAintersect(unique(sigpw_genes$Gene),names(which(table(TCGA_LUSC_maf_data_select$Hugo_Symbol) > 6)))

TCGA_exon_gene_sigpw <- intersect(unique(sigpw_genes$Gene),names(which(table(TCGA_LUSC_maf_data_select$Hugo_Symbol) > 6)))


TCGA_LUSC_mut_info <- genesToBarcodes(maf = TCGA_LUSC_maf,genes = TCGA_exon_gene_sigpw)

for (i in TCGA_exon_gene_sigpw) {
  colnames(TCGA_LUSC_mut_info[[i]]) <- c("patients","Frame_Shift_Del","Frame_Shift_Ins",
                                         "In_Frame_Del","In_Frame_Ins","Missense","Nonsense",
                                         "Nonstop","Splice_Site","Translation_Start_Site","total")
  TCGA_LUSC_mut_info[[i]][,"patients"] <- substr(as.vector(as.data.frame(TCGA_LUSC_mut_info[[i]])[,"patients"]),1,12)
}

TCGA_maf_all_patients <- make.names(unique(substr(TCGA_LUSC_maf@clinical.data$Tumor_Sample_Barcode,1,12)))

TCGA_Clinical.tidy.lusc.maf <- TCGA_Clinical.tidy.lusc[TCGA_Clinical.tidy.lusc$patients%in%TCGA_maf_all_patients,]

TCGA_Clinical.tidy.lusc.maf.2 <- TCGA_Clinical.tidy.lusc.maf[,1:8]

tcga_exon_gene_patients_list <- list()

for (i in TCGA_exon_gene_sigpw) {
  TCGA_Clinical.tidy.lusc.maf[make.names(unique(TCGA_LUSC_mut_info[[i]]$patients)),i] <- 1
  TCGA_Clinical.tidy.lusc.maf[is.na(TCGA_Clinical.tidy.lusc.maf[,i]),i] <- 0
  
  tcga_exon_gene_patients_list[[i]] <- make.names(unique(TCGA_LUSC_mut_info[[i]]$patients))
}

tcga_exon_gene_patients_list_length <- unlist(lapply(tcga_exon_gene_patients_list, length))


#Analysis of single-gene exon mutations in TCGA

tcga_exon_gene_multicox_os_survial <- run_multicox_OS(make.names(TCGA_exon_gene_sigpw),TCGA_Clinical.tidy.lusc.maf,cli_factor = c("Age", "Stage", "Smoking_history"))
tcga_exon_gene_multicox_rfs_survial <- run_multicox_RFS(make.names(TCGA_exon_gene_sigpw),TCGA_Clinical.tidy.lusc.maf,cli_factor = c("Age", "Stage", "Smoking_history"))
tcga_exon_gene_multicox_rfs_os_merge_survial <- merge(tcga_exon_gene_multicox_os_survial,tcga_exon_gene_multicox_rfs_survial,by="Gene")
tcga_exon_gene_multicox_rfs_os_merge_survial$cohort <- "TCGA (n=492)"
tcga_exon_gene_multicox_rfs_os_merge_survial$Mutation_Number <- as.vector(unlist(tcga_exon_gene_patients_list_length[tcga_exon_gene_multicox_rfs_os_merge_survial$Gene]))

#Prognostic analysis of the TCGA exon mutation pathway
tcga_exon_pathway_mat_sigpw <- analyze_signature(TCGA_Clinical.tidy.lusc.maf,TCGA_exon_gene_sigpw,signatureset_sigpw,0)


tcga_exon_pathway_mat_sigpw_clinical <- as.data.frame(t(tcga_exon_pathway_mat_sigpw)) %>%  tibble::rownames_to_column(var = "patients")
tcga_exon_pathway_mat_sigpw_clinical_merge <- merge(TCGA_Clinical.tidy.lusc.maf.2,tcga_exon_pathway_mat_sigpw_clinical,by="patients",all="TRUE")


colnames(tcga_exon_pathway_mat_sigpw_clinical_merge) <- make.names(colnames(tcga_exon_pathway_mat_sigpw_clinical_merge))
rownames(tcga_exon_pathway_mat_sigpw_clinical_merge) <- tcga_exon_pathway_mat_sigpw_clinical_merge$patients

tcga_exon_pathway_mat_sigpw_clinical_merge_length <- colSums(tcga_exon_pathway_mat_sigpw_clinical_merge[,9:ncol(tcga_exon_pathway_mat_sigpw_clinical_merge)])


exon_tcga_pathway_multicox_os_survial <- run_multicox_OS(make.names(rownames(tcga_exon_pathway_mat_sigpw)),tcga_exon_pathway_mat_sigpw_clinical_merge,
                                                         cli_factor = c("Age", "Stage", "Smoking_history"))

exon_tcga_pathway_multicox_rfs_survial <- run_multicox_RFS(make.names(rownames(tcga_exon_pathway_mat_sigpw)),tcga_exon_pathway_mat_sigpw_clinical_merge,
                                                           cli_factor = c("Age", "Stage", "Smoking_history"))

exon_tcga_pathway_multicox_os_rfs_merge_survial <- merge(exon_tcga_pathway_multicox_os_survial,exon_tcga_pathway_multicox_rfs_survial,by="Gene")
exon_tcga_pathway_multicox_os_rfs_merge_survial$cohort <- "TCGA (n=492)"
exon_tcga_pathway_multicox_os_rfs_merge_survial$Mutation_Number <- as.vector(unlist(tcga_exon_pathway_mat_sigpw_clinical_merge_length[exon_tcga_pathway_multicox_os_rfs_merge_survial$Gene]))




# cnv

lusc.109.gistic <- readGistic(gisticAllLesionsFile="CNV/109/result_109/521405/all_lesions.conf_90.txt",
                              gisticAmpGenesFile="CNV/109/result_109/521405/amp_genes.conf_90.txt",
                              gisticDelGenesFile="CNV/109/result_109/521405/del_genes.conf_90.txt",
                              gisticScoresFile="CNV/109/result_109/521405/scores.gistic",
                              isTCGA=FALSE)

cnv.summary <- as.data.frame(lusc.109.gistic@cnv.summary)
rownames(cnv.summary) <- cnv.summary$Tumor_Sample_Barcode
rownames(cnv.summary)[rownames(cnv.summary) %in% "133666"] <- "1336663"


cnMatrix <- as.data.frame(lusc.109.gistic@cnMatrix)

cnv.data <- as.data.frame(lusc.109.gistic@data)
cnv.data$Tumor_Sample_Barcode <- as.vector(cnv.data$Tumor_Sample_Barcode)

cnv.data[cnv.data$Tumor_Sample_Barcode %in% "133666","Tumor_Sample_Barcode"] <- "1336663"


gisticChromPlot(gistic = lusc.109.gistic, markBands = "all")
gisticChromPlot(gistic=lusc.109.gistic)

select_cytoband <- as.vector(unlist(lusc.109.gistic@cytoband.summary[order(lusc.109.gistic@cytoband.summary$qvalues),"Cytoband"][1:15]))

pdf("cnv.pdf",width = 7,height = 5)
gisticChromPlot(gistic = lusc.109.gistic,ref.build = "hg19",
                maf = LUSC_maf,
                fdrCutOff =  0.1,
                txtSize = 1.3,
                markBands = select_cytoband,
                mutGenes = c("ERBB4","SOX2","FGFR1","PIK3CA","CDKN2A"),
                mutGenesTxtSize = 1.3)
dev.off()



#out
clinical_109_2 <- clinical_109
colnames(clinical_109_2)[1] <- "Tumor_Sample_Barcode"
colnames(clinical_109_2)[4] <- "Smoking"
clinical_109_2[clinical_109_2$Smoking%in%1,"Smoking"] <- "Yes"
clinical_109_2[clinical_109_2$Smoking%in%0,"Smoking"] <- "No"



LUSC_maf_clinical <- read.maf("data/LUSC.somatic.maf",
                              vc_nonSyn = c(exon_mutation_types,promoter_mutation_types),
                              clinicalData = "data/109_clinical_2.xls"
                              
)


LUSC_maf_promoter <- LUSC_maf_clinical
LUSC_maf_exon <- LUSC_maf_clinical

exon_LUSC_maf <- subsetMaf(maf = LUSC_maf_clinical,
                           query="Variant_Classification %in% c('Missense_Mutation','Nonsense_Mutation',
                           'Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site')")

promoter_flag <- c("5'UTR","5'Flank")

promoter_LUSC_maf <- subsetMaf(maf = LUSC_maf_clinical, query = "Variant_Classification %in% promoter_flag",)


#Figure1


pdf("Figure_1A.pdf",width = 7,height = 8)

oncoplot(maf = exon_LUSC_maf,
         top = 20,
         clinicalFeatures = c("Gender","Age","Stage","Smoking"),
         fontSize = 1.25,
         titleFontSize = 2.5,
         legendFontSize = 1.8,
         annotationFontSize = 1.8,
         groupAnnotationBySize = TRUE,
         SampleNamefontSize = 1,
         anno_height = 2.5,
         legend_height = 5,
         showPct = TRUE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         logColBar =  FALSE,
         showTitle = TRUE,
         bgCol = "white"

)
dev.off()


pdf("Figure_1C.pdf",width = 7,height = 8)

oncoplot(maf = promoter_LUSC_maf,
         top = 20,
         clinicalFeatures = c("Gender","Age","Stage","Smoking"),
         fontSize = 1.25,
         titleFontSize = 2.5,
         legendFontSize = 1.8,
         annotationFontSize = 1.8,
         groupAnnotationBySize = TRUE,
         SampleNamefontSize = 1,
         anno_height = 2.5,
         legend_height = 5,
         showPct = TRUE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         logColBar =  FALSE,
         showTitle = TRUE,
         bgCol = "white"
)

dev.off()


library(eulerr)

calculate_custom_sets <- function(A, B, C, 
                                  nameA = "A", 
                                  nameB = "B", 
                                  nameC = "C", 
                                  nameA_and_B = "A&B", 
                                  nameA_and_C = "A&C", 
                                  nameB_and_C = "B&C", 
                                  nameA_and_B_and_C = "A&B&C") {
  
  onlyA <- setdiff(A, union(B, C))
  onlyB <- setdiff(B, union(A, C))
  onlyC <- setdiff(C, union(A, B))
  
  A_and_B <- intersect(A, B)
  A_and_C <- intersect(A, C)
  B_and_C <- intersect(B, C)
  
  A_and_B_only <- setdiff(A_and_B, C)
  A_and_C_only <- setdiff(A_and_C, B)
  B_and_C_only <- setdiff(B_and_C, A)
  
  A_and_B_and_C <- intersect(A_and_B, C)
  
  result <- list()
  result[[nameA]] <- onlyA
  result[[nameB]] <- onlyB
  result[[nameC]] <- onlyC
  result[[nameA_and_B]] <- A_and_B_only
  result[[nameA_and_C]] <- A_and_C_only
  result[[nameB_and_C]] <- B_and_C_only
  result[[nameA_and_B_and_C]] <- A_and_B_and_C
  
  lengths <- sapply(result, length)
  
  return(list(result = result, lengths = lengths))
}


venn_exon_mutation_genes <- unique(list_109_exon[[2]][list_109_exon[[2]][,3]>0,1])
venn_promoter_mutation_genes <- unique(list_109[[2]][list_109[[2]][,3]>0,1])
venn_cnv_mutation_genes <- unique(names(which(table(cnv.data$Hugo_Symbol)>0)))
venn_cnv_mutation_genes <- venn_cnv_mutation_genes[!grepl("hsa-",venn_cnv_mutation_genes)]
venn_cnv_mutation_genes <- gsub("\\[|\\]", "", venn_cnv_mutation_genes)


exon_promoter_cnv_venn <- calculate_custom_sets(venn_exon_mutation_genes, 
                                                venn_promoter_mutation_genes, 
                                                venn_cnv_mutation_genes, 
                                                nameA = "exon", 
                                                nameB = "promoter", 
                                                nameC = "CNV", 
                                                nameA_and_B = "exon&promoter", 
                                                nameA_and_C = "exon&CNV", 
                                                nameB_and_C = "promoter&CNV", 
                                                nameA_and_B_and_C = "exon&promoter&CNV")


pdf("Figure_1D.pdf",width = 5,height = 5)

plot(euler(exon_promoter_cnv_venn$lengths, shape = "circle"), 
     list(fill = c("#39a044", "#002a68","#bebcbc"), alpha = 0.9),
     quantities = list(type = "counts", cex=1.3, col = c("white")),  
     legend = list(side = "bottom"))


dev.off()




exon_TCGA_LUSC_maf <- subsetMaf(maf = TCGA_LUSC_maf,
                           query="Variant_Classification %in% c('Missense_Mutation','Nonsense_Mutation',
                           'Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins','Splice_Site')")

pdf("Figure_S1A.pdf",width = 14,height = 12)
oncoplot(exon_TCGA_LUSC_maf,
         top=20,
         fontSize = 2,
         titleFontSize = 2.5,
         legendFontSize = 2.7,
         annotationFontSize = 2.7,
         groupAnnotationBySize = TRUE,
         SampleNamefontSize = 1,
         anno_height = 2.5,
         legend_height = 5,
         showPct = TRUE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         logColBar =  FALSE,
         showTitle = TRUE,
         bgCol = "white"
         )
dev.off()




TCGA_venn_exon_mutation_genes <- unique(exon_TCGA_LUSC_maf@data$Hugo_Symbol)


pdf("Figure_S1B.pdf",width = 5,height = 5)

plot(euler(c("TCGA"=length(TCGA_venn_exon_mutation_genes),
             "109"=length(venn_exon_mutation_genes),
             "TCGA&109"=length(intersect(TCGA_venn_exon_mutation_genes,venn_exon_mutation_genes))
             ), shape = "circle"), 
     list(fill = c("#39a044", "#002a68"), alpha = 0.9),
     quantities = list(type = "counts", cex=1.3, col = c("white")),  
     legend = list(side = "bottom",cex=1.3))

dev.off()


#


pdf("Figure_2A.pdf",width = 7,height = 9)
oncoplot(exon_LUSC_maf,
         genes = exon_gene_sigpw,
         fontSize = 1.25,
         titleFontSize = 2.5,
         legendFontSize = 2,
         annotationFontSize = 2.7,
         groupAnnotationBySize = TRUE,
         SampleNamefontSize = 1,
         anno_height = 2.5,
         legend_height = 5,
         showPct = TRUE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         logColBar =  FALSE,
         showTitle = TRUE,
         bgCol = "white"
         )
dev.off()


Heatmap_color <- c("white","black")
names(Heatmap_color) <- c(0,1)

exon_oncoplot_pathway_109 <- ComplexHeatmap::Heatmap(exon_pathway_mat_sigpw,
                                     name = "109 exon pathway mutaion",
                                     column_title = "109 exon pathway mutaion",
                                     show_column_names = FALSE,
                                     show_heatmap_legend = FALSE,
                                     show_column_dend = FALSE,
                                     show_row_dend = FALSE,
                                     cluster_rows = TRUE,
                                     cluster_columns = TRUE,
                                     row_names_side = "left",
                                     row_names_gp = gpar(fontsize = 16),
                                     col = Heatmap_color,
                                     rect_gp = gpar(col = "white", lwd = 1),
                                     column_title_gp =  gpar(col = "black", fontsize = 18,fontface="bold")
)

exon_oncoplot_pathway_109_row_percent <- round(exon_pathway_mat_sigpw_clinical_merge_length[make.names(rownames(exon_pathway_mat_sigpw))]/109*100)

custom_legend <- Legend(title = "mutaion status",
                        at = c(0, 1), 
                        labels = c("WT","MUT"),
                        direction = "vertical",
                        legend_gp = gpar(fill = c("white","#34495E")),
                        title_gp = gpar(fontsize = 13, fontface = "bold"),
                        labels_gp = gpar(fontsize = 13),
                        title_gap = unit(2, "mm"),
                        row_gap = unit(1, "mm")
                        )

row_109_ha <- rowAnnotation(Percentage = anno_text(paste0(exon_oncoplot_pathway_109_row_percent, "%"),
                                                   gp = gpar(fontsize = 14),
                                                   just = "right",
                                                   location = 1,
                                                   width = unit(1, "cm"))
)



pdf("Figure_2B.pdf",width = 8,height = 8)

draw(exon_oncoplot_pathway_109 + row_109_ha,annotation_legend_list = list(custom_legend))

dev.off()


pdf("Figure_2C.pdf",width = 7,height = 9)

oncoplot(maf = promoter_LUSC_maf,
         genes = promoter_gene_sigpw,
         fontSize = 1.25,
         titleFontSize = 2.5,
         legendFontSize = 2,
         annotationFontSize = 2.7,
         groupAnnotationBySize = TRUE,
         SampleNamefontSize = 1,
         anno_height = 2.5,
         legend_height = 5,
         showPct = TRUE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         logColBar =  FALSE,
         showTitle = TRUE,
         bgCol = "white"
         )

dev.off()



promoter_oncoplot_pathway_109 <- Heatmap(promoter_pathway_mat_sigpw,
                                         name = "109 promoter pathway mutaion",
                                         column_title = "109 promoter pathway mutaion",
                                         show_column_names = FALSE,
                                         show_heatmap_legend = FALSE,
                                         show_column_dend = FALSE,
                                         show_row_dend = FALSE,
                                         cluster_rows = TRUE,
                                         cluster_columns = TRUE,
                                         row_names_side = "left",
                                         row_names_gp = gpar(fontsize = 16),
                                         col = Heatmap_color,
                                         rect_gp = gpar(col = "white", lwd = 1),
                                         column_title_gp =  gpar(col = "black", fontsize = 18,fontface="bold")
)

promoter_oncoplot_pathway_109_row_percent <- round(promoter_pathway_mat_sigpw_clinical_merge_length[make.names(rownames(promoter_pathway_mat_sigpw))]/109*100)

custom_legend <- Legend(title = "mutaion status",
                        at = c(0, 1), 
                        labels = c("WT","MUT"),
                        direction = "vertical",
                        legend_gp = gpar(fill = c("white","#34495E")),
                        title_gp = gpar(fontsize = 13, fontface = "bold"),
                        labels_gp = gpar(fontsize = 13),
                        title_gap = unit(2, "mm"),
                        row_gap = unit(1, "mm")
                        )

row_109_ha <- rowAnnotation(Percentage = anno_text(paste0(promoter_oncoplot_pathway_109_row_percent, "%"),
                                                   gp = gpar(fontsize = 14),
                                                   just = "right",
                                                   location = 1,
                                                   width = unit(1, "cm"))
)



pdf("Figure_2D.pdf",width = 8,height = 8)

draw(promoter_oncoplot_pathway_109 + row_109_ha,annotation_legend_list = list(custom_legend))

dev.off()


#FigureS2A
pdf("Figure_S2A.pdf",width = 14,height = 14)

oncoplot(TCGA_LUSC_maf,
         genes = exon_gene_sigpw,
         fontSize = 2,
         titleFontSize = 2.5,
         legendFontSize = 2.7,
         annotationFontSize = 2.7,
         groupAnnotationBySize = TRUE,
         SampleNamefontSize = 1,
         anno_height = 2.5,
         legend_height = 5,
         showPct = TRUE,
         drawRowBar = TRUE,
         drawColBar = TRUE,
         logColBar =  FALSE,
         showTitle = TRUE,
         bgCol = "white"
         )

dev.off()


#FigureS2B


exon_oncoplot_pathway_tcga <- Heatmap(tcga_exon_pathway_mat_sigpw,
                                      name = "TCGA exon pathway mutaion",
                                      column_title = "TCGA exon pathway mutaion",
                                      show_column_names = FALSE,
                                      show_heatmap_legend = FALSE,
                                      show_column_dend = FALSE,
                                      show_row_dend = FALSE,
                                      cluster_rows = TRUE,
                                      cluster_columns = TRUE,
                                      row_names_side = "left",
                                      row_names_gp = gpar(fontsize = 16),
                                      col = Heatmap_color,
                                      rect_gp = gpar(col = "white", lwd = 1),
                                      column_title_gp =  gpar(col = "black", fontsize = 18,fontface="bold")
)

exon_oncoplot_pathway_tcga_row_percent <- round(tcga_exon_pathway_mat_sigpw_clinical_merge_length[make.names(rownames(tcga_exon_pathway_mat_sigpw))]/492*100)

custom_legend <- Legend(title = "mutaion status",
                        at = c(0, 1), 
                        labels = c("WT","MUT"),
                        direction = "vertical",
                        legend_gp = gpar(fill = c("white","black")),
                        title_gp = gpar(fontsize = 13, fontface = "bold"),
                        labels_gp = gpar(fontsize = 13),
                        title_gap = unit(2, "mm"),
                        row_gap = unit(1, "mm")
)

row_tcga_ha <- rowAnnotation(Percentage = anno_text(paste0(exon_oncoplot_pathway_tcga_row_percent, "%"),
                                                    gp = gpar(fontsize = 14),
                                                    just = "right",
                                                    location = 1,
                                                    width = unit(1, "cm"))
)



pdf("Figure_S2B.pdf",width = 14,height = 6)

draw(exon_oncoplot_pathway_tcga + row_tcga_ha,annotation_legend_list = list(custom_legend))

dev.off()


#Figure3
#Survival analysis of OS in exon genes of cohort 109

survival_plot_os_109_list <- list()

for (i in exon_gene_sigpw) {
  
  multi_data_exon_temp <- multi_data_exon[,c("OS.months.","Survival.status",i)]
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_temp[,i] <- factor(multi_data_exon_temp[,i],levels = unique(multi_data_exon_temp[,i]))
  multi_data_exon_temp[,i] <- relevel(multi_data_exon_temp[,i],ref="WT")
  

  coxph_model <- coxph(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_exon_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_exon_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_os_109_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                               palette = palette_colors,
                                               size=0.6,
                                               censor.shape= "|",
                                               censor.size= 3,
                                               ggtheme = theme_classic()+ theme(
                                                 axis.line = element_line(size = 0.35,colour = "black"),
                                                 axis.ticks = element_line(size = 0.35),
                                                 axis.ticks.length = unit(0.25, "cm"),
                                                 axis.text  = element_text(size = 15,face = "bold"),
                                                 axis.title = element_text(size = 20),
                                                 axis.title.x = element_text(size = 18),
                                                 axis.title.y = element_text(size = 20),
                                                 legend.direction  = "horizontal",
                                                 legend.background = element_blank(),
                                                 legend.key.size = unit(10,"mm"),
                                                 legend.text =element_text(size = 14),
                                                 legend.title = element_text(size = 15,face = "bold"),
                                                 plot.title  = element_text(size = 20)
                                                 
                                               
                                               ),
                                               tables.theme = theme(text = element_text(size = 15),
                                                                    axis.text  = element_text(size = 15),
                                                                    axis.line = element_line(size = 0.35,colour = "black"),
                                                                    axis.ticks = element_line(size = 0.35),
                                                                    axis.ticks.length = unit(0.25,"mm"),
                                                                    axis.title = element_text(size = 20),
                                                                    plot.title  = element_text(size = 15)),
  
                                               legend.title = "",
                                               ylab = "Overall survival probability",
                                               xlab = "months",
                                               data=multi_data_exon_temp)
  

  survival_plot_os_109_list[[i]]$plot <- survival_plot_os_109_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(multi_data_exon_temp[,"OS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(multi_data_exon_temp[,"OS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(multi_data_exon_temp[,"OS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)

  
}


for (i in names(survival_plot_os_109_list)) {
  pdf(paste0("Figure_3_exon_OS_109_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_os_109_list[[i]],newpage = FALSE)
  dev.off()
}

#Survival analysis of RFS for exon genes in cohort 109


survival_plot_rfs_109_list <- list()

for (i in exon_gene_sigpw) {
  
  multi_data_exon_temp <- multi_data_exon[,c("RFS.months.","Recurrence.status",i)]
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_temp[,i] <- factor(multi_data_exon_temp[,i],levels = unique(multi_data_exon_temp[,i]))
  multi_data_exon_temp[,i] <- relevel(multi_data_exon_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_exon_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_exon_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_rfs_109_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                palette = palette_colors,
                                                size=0.6,
                                                censor.shape= "|",
                                                censor.size= 3,
                                                
                                                ggtheme = theme_classic()+ theme(
                                                  axis.line = element_line(size = 0.35,colour = "black"),
                                                  axis.ticks = element_line(size = 0.35),
                                                  axis.ticks.length = unit(0.25, "cm"),
                                                  axis.text  = element_text(size = 15,face = "bold"),
                                                  axis.title = element_text(size = 20),
                                                  axis.title.x = element_text(size = 18),
                                                  axis.title.y = element_text(size = 20),
                                                  legend.direction  = "horizontal",
                                                  legend.background = element_blank(),
                                                  legend.key.size = unit(10,"mm"),
                                                  legend.text =element_text(size = 14),
                                                  legend.title = element_text(size = 15,face = "bold"),
                                                  plot.title  = element_text(size = 20)
                                                  
                                                  
                                                ),
                                                tables.theme = theme(text = element_text(size = 15),
                                                                     axis.text  = element_text(size = 15),
                                                                     axis.line = element_line(size = 0.35,colour = "black"),
                                                                     axis.ticks = element_line(size = 0.35),
                                                                     axis.ticks.length = unit(0.25,"mm"),
                                                                     axis.title = element_text(size = 20),
                                                                     plot.title  = element_text(size = 15)),
                                                
                                                legend.title = "",
                                                ylab = "Recurrence survival probability",
                                                xlab = "months",
                                                data=multi_data_exon_temp)
  
  
  survival_plot_rfs_109_list[[i]]$plot <- survival_plot_rfs_109_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(multi_data_exon_temp[,"RFS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(multi_data_exon_temp[,"RFS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(multi_data_exon_temp[,"RFS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
  
  
}


for (i in names(survival_plot_rfs_109_list)) {
  pdf(paste0("Figure_3_exon_RFS_109_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_rfs_109_list[[i]],newpage = FALSE)
  dev.off()
}

#OS survival analysis of exon pathways in cohort 109


survival_plot_os_109_pathway_list <- list()

for (i in make.names(rownames(exon_pathway_mat_sigpw))) {
  
  exon_pathway_mat_sigpw_clinical_merge_temp <- exon_pathway_mat_sigpw_clinical_merge[,c("OS.months.","Survival.status",i)]
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = exon_pathway_mat_sigpw_clinical_merge_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = exon_pathway_mat_sigpw_clinical_merge_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_os_109_pathway_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                       palette = palette_colors,
                                                       legend.title = "",
                                                       size=0.6,
                                                       censor.shape= "|",
                                                       censor.size= 3,
                                                       
                                                       ggtheme = theme_classic()+ theme(
                                                         axis.line = element_line(size = 0.35,colour = "black"),
                                                         axis.ticks = element_line(size = 0.35),
                                                         axis.ticks.length = unit(0.25, "cm"),
                                                         axis.text  = element_text(size = 15,face = "bold"),
                                                         axis.title = element_text(size = 20),
                                                         axis.title.x = element_text(size = 18),
                                                         axis.title.y = element_text(size = 20),
                                                         legend.direction  = "horizontal",
                                                         legend.background = element_blank(),
                                                         legend.key.size = unit(10,"mm"),
                                                         legend.text =element_text(size = 14),
                                                         legend.title = element_text(size = 15,face = "bold"),
                                                         plot.title  = element_text(size = 20)
                                                         
                                                         
                                                       ),
                                                       tables.theme = theme(text = element_text(size = 15),
                                                                            axis.text  = element_text(size = 15),
                                                                            axis.line = element_line(size = 0.35,colour = "black"),
                                                                            axis.ticks = element_line(size = 0.35),
                                                                            axis.ticks.length = unit(0.25,"mm"),
                                                                            axis.title = element_text(size = 20),
                                                                            plot.title  = element_text(size = 15)),
                                                       
                                                       ylab = "Overall survival probability",
                                                       xlab = "months",
                                                       data=exon_pathway_mat_sigpw_clinical_merge_temp)
  
  
  survival_plot_os_109_pathway_list[[i]]$plot <- survival_plot_os_109_pathway_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(exon_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(exon_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(exon_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
  
  
}


for (i in names(survival_plot_os_109_pathway_list)) {
  pdf(paste0("Figure_3_exon_OS_109_pathway_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_os_109_pathway_list[[i]],newpage = FALSE)
  dev.off()
}

#Survival analysis of RFS in the 109 cohort exon pathway

survival_plot_rfs_109_pathway_list <- list()

for (i in make.names(rownames(exon_pathway_mat_sigpw))) {
  
  exon_pathway_mat_sigpw_clinical_merge_temp <- exon_pathway_mat_sigpw_clinical_merge[,c("RFS.months.","Recurrence.status",i)]
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = exon_pathway_mat_sigpw_clinical_merge_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = exon_pathway_mat_sigpw_clinical_merge_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_rfs_109_pathway_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                palette = palette_colors,
                                                size=0.6,
                                                legend.title = "",
                                                censor.shape= "|",
                                                censor.size= 3,
                                                
                                                ggtheme = theme_classic()+ theme(
                                                  axis.line = element_line(size = 0.35,colour = "black"),
                                                  axis.ticks = element_line(size = 0.35),
                                                  axis.ticks.length = unit(0.25, "cm"),
                                                  axis.text  = element_text(size = 15,face = "bold"),
                                                  axis.title = element_text(size = 20),
                                                  axis.title.x = element_text(size = 18),
                                                  axis.title.y = element_text(size = 20),
                                                  legend.direction  = "horizontal",
                                                  legend.background = element_blank(),
                                                  legend.key.size = unit(10,"mm"),
                                                  legend.text =element_text(size = 14),
                                                  legend.title = element_text(size = 15,face = "bold"),
                                                  plot.title  = element_text(size = 20)
                                                  
                                                  
                                                ),
                                                tables.theme = theme(text = element_text(size = 15),
                                                                     axis.text  = element_text(size = 15),
                                                                     axis.line = element_line(size = 0.35,colour = "black"),
                                                                     axis.ticks = element_line(size = 0.35),
                                                                     axis.ticks.length = unit(0.25,"mm"),
                                                                     axis.title = element_text(size = 20),
                                                                     plot.title  = element_text(size = 15)),
                                                
                                                ylab = "Recurrence survival probability",
                                                xlab = "months",
                                                data=exon_pathway_mat_sigpw_clinical_merge_temp)
  
  
  survival_plot_rfs_109_pathway_list[[i]]$plot <- survival_plot_rfs_109_pathway_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(exon_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(exon_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(exon_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
  
  
}


for (i in names(survival_plot_rfs_109_pathway_list)) {
  pdf(paste0("Figure_3_exon_RFS_109_pathway_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_rfs_109_pathway_list[[i]],newpage = FALSE)
  dev.off()
}




#Survival analysis of OS in the promoter genes of cohort 109

survival_plot_os_109_promoter_list <- list()

for (i in promoter_gene_sigpw) {
  
  multi_data_promoter_temp <- multi_data_promoter[,c("OS.months.","Survival.status",i)]
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%1,i] <- "MUT"
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%0,i] <- "WT"
  multi_data_promoter_temp[,i] <- factor(multi_data_promoter_temp[,i],levels = unique(multi_data_promoter_temp[,i]))
  multi_data_promoter_temp[,i] <- relevel(multi_data_promoter_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_promoter_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_promoter_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_os_109_promoter_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                               palette = palette_colors,
                                               size=0.6,
                                               legend.title = "",
                                               censor.shape= "|",
                                               censor.size= 3,
                                               
                                               ggtheme = theme_classic()+ theme(
                                                 axis.line = element_line(size = 0.35,colour = "black"),
                                                 axis.ticks = element_line(size = 0.35),
                                                 axis.ticks.length = unit(0.25, "cm"),
                                                 axis.text  = element_text(size = 15,face = "bold"),
                                                 axis.title = element_text(size = 20),
                                                 axis.title.x = element_text(size = 18),
                                                 axis.title.y = element_text(size = 20),
                                                 legend.direction  = "horizontal",
                                                 legend.background = element_blank(),
                                                 legend.key.size = unit(10,"mm"),
                                                 legend.text =element_text(size = 14),
                                                 legend.title = element_text(size = 15,face = "bold"),
                                                 plot.title  = element_text(size = 20)
                                                 
                                                 
                                               ),
                                               tables.theme = theme(text = element_text(size = 15),
                                                                    axis.text  = element_text(size = 15),
                                                                    axis.line = element_line(size = 0.35,colour = "black"),
                                                                    axis.ticks = element_line(size = 0.35),
                                                                    axis.ticks.length = unit(0.25,"mm"),
                                                                    axis.title = element_text(size = 20),
                                                                    plot.title  = element_text(size = 15)),
                                               
                                               ylab = "Overall survival probability",
                                               xlab = "months",
                                               data=multi_data_promoter_temp)
  
  
  survival_plot_os_109_promoter_list[[i]]$plot <- survival_plot_os_109_promoter_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(multi_data_promoter_temp[,"OS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(multi_data_promoter_temp[,"OS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(multi_data_promoter_temp[,"OS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
    # 
  
}




for (i in names(survival_plot_os_109_promoter_list)) {
  pdf(paste0("Figure_3_promoter_OS_109_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_os_109_promoter_list[[i]],newpage = FALSE)
  dev.off()
}

#Survival analysis of RFS for promoter genes in cohort 109


survival_plot_rfs_109_promoter_list <- list()

for (i in promoter_gene_sigpw) {
  
  multi_data_promoter_temp <- multi_data_promoter[,c("RFS.months.","Recurrence.status",i)]
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%1,i] <- "MUT"
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%0,i] <- "WT"
  multi_data_promoter_temp[,i] <- factor(multi_data_promoter_temp[,i],levels = unique(multi_data_promoter_temp[,i]))
  multi_data_promoter_temp[,i] <- relevel(multi_data_promoter_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_promoter_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_promoter_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_rfs_109_promoter_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                         palette = palette_colors,
                                                         size=0.6,
                                                         legend.title = "",
                                                         censor.shape= "|",
                                                         censor.size= 3,
                                                         
                                                         ggtheme = theme_classic()+ theme(
                                                           axis.line = element_line(size = 0.35,colour = "black"),
                                                           axis.ticks = element_line(size = 0.35),
                                                           axis.ticks.length = unit(0.25, "cm"),
                                                           axis.text  = element_text(size = 15,face = "bold"),
                                                           axis.title = element_text(size = 20),
                                                           axis.title.x = element_text(size = 18),
                                                           axis.title.y = element_text(size = 20),
                                                           legend.direction  = "horizontal",
                                                           legend.background = element_blank(),
                                                           legend.key.size = unit(10,"mm"),
                                                           legend.text =element_text(size = 14),
                                                           legend.title = element_text(size = 15,face = "bold"),
                                                           plot.title  = element_text(size = 20)
                                                           
                                                           
                                                         ),
                                                         tables.theme = theme(text = element_text(size = 15),
                                                                              axis.text  = element_text(size = 15),
                                                                              axis.line = element_line(size = 0.35,colour = "black"),
                                                                              axis.ticks = element_line(size = 0.35),
                                                                              axis.ticks.length = unit(0.25,"mm"),
                                                                              axis.title = element_text(size = 20),
                                                                              plot.title  = element_text(size = 15)),
                                                         
                                                         ylab = "Recurrence survival probability",
                                                         xlab = "months",
                                                         data=multi_data_promoter_temp)
  
  
  survival_plot_rfs_109_promoter_list[[i]]$plot <- survival_plot_rfs_109_promoter_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(multi_data_promoter_temp[,"RFS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(multi_data_promoter_temp[,"RFS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(multi_data_promoter_temp[,"RFS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
  
  
}


for (i in names(survival_plot_rfs_109_promoter_list)) {
  pdf(paste0("Figure_3_promoter_RFS_109_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_rfs_109_promoter_list[[i]],newpage = FALSE)
  dev.off()
}



#Survival analysis of OS in the 109 cohort promoter pathway


survival_plot_os_109_promoter_pathway_list <- list()

for (i in make.names(rownames(promoter_pathway_mat_sigpw))) {
  
  promoter_pathway_mat_sigpw_clinical_merge_temp <- promoter_pathway_mat_sigpw_clinical_merge[,c("OS.months.","Survival.status",i)]
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(promoter_pathway_mat_sigpw_clinical_merge_temp[,i]))
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = promoter_pathway_mat_sigpw_clinical_merge_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = promoter_pathway_mat_sigpw_clinical_merge_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_os_109_promoter_pathway_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                        palette = palette_colors,
                                                        size=0.6,
                                                        legend.title = "",
                                                        censor.shape= "|",
                                                        censor.size= 3,
                                                        ggtheme = theme_classic()+ theme(
                                                          axis.line = element_line(size = 0.35,colour = "black"),
                                                          axis.ticks = element_line(size = 0.35),
                                                          axis.ticks.length = unit(0.25, "cm"),
                                                          axis.text  = element_text(size = 15,face = "bold"),
                                                          axis.title = element_text(size = 20),
                                                          axis.title.x = element_text(size = 18),
                                                          axis.title.y = element_text(size = 20),
                                                          legend.direction  = "horizontal",
                                                          legend.background = element_blank(),
                                                          legend.key.size = unit(10,"mm"),
                                                          legend.text =element_text(size = 14),
                                                          legend.title = element_text(size = 15,face = "bold"),
                                                          plot.title  = element_text(size = 20)
                                                          
                                                          
                                                        ),
                                                        tables.theme = theme(text = element_text(size = 15),
                                                                             axis.text  = element_text(size = 15),
                                                                             axis.line = element_line(size = 0.35,colour = "black"),
                                                                             axis.ticks = element_line(size = 0.35),
                                                                             axis.ticks.length = unit(0.25,"mm"),
                                                                             axis.title = element_text(size = 20),
                                                                             plot.title  = element_text(size = 15)),
                                                        
                                                        ylab = "Overall survival probability",
                                                        xlab = "months",
                                                        data=promoter_pathway_mat_sigpw_clinical_merge_temp)
  
  
  survival_plot_os_109_promoter_pathway_list[[i]]$plot <- survival_plot_os_109_promoter_pathway_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(promoter_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(promoter_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(promoter_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
    # 
  
}

for (i in names(survival_plot_os_109_promoter_pathway_list)) {
  pdf(paste0("Figure_3_promoter_OS_109_pathway_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_os_109_promoter_pathway_list[[i]],newpage = FALSE)
  dev.off()
}




#Survival analysis of RFS in the 109 cohort exon pathway


survival_plot_rfs_109_promoter_pathway_list <- list()

for (i in make.names(rownames(promoter_pathway_mat_sigpw))) {
  
  promoter_pathway_mat_sigpw_clinical_merge_temp <- promoter_pathway_mat_sigpw_clinical_merge[,c("RFS.months.","Recurrence.status",i)]
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(promoter_pathway_mat_sigpw_clinical_merge_temp[,i]))
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = promoter_pathway_mat_sigpw_clinical_merge_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = promoter_pathway_mat_sigpw_clinical_merge_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_rfs_109_promoter_pathway_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                                 palette = palette_colors,
                                                                 size=0.6,
                                                                 legend.title = "",
                                                                 censor.shape= "|",
                                                                 censor.size= 3,
                                                                 ggtheme = theme_classic()+ theme(
                                                                   axis.line = element_line(size = 0.35,colour = "black"),
                                                                   axis.ticks = element_line(size = 0.35),
                                                                   axis.ticks.length = unit(0.25, "cm"),
                                                                   axis.text  = element_text(size = 15,face = "bold"),
                                                                   axis.title = element_text(size = 20),
                                                                   axis.title.x = element_text(size = 18),
                                                                   axis.title.y = element_text(size = 20),
                                                                   legend.direction  = "horizontal",
                                                                   legend.background = element_blank(),
                                                                   legend.key.size = unit(10,"mm"),
                                                                   legend.text =element_text(size = 14),
                                                                   legend.title = element_text(size = 15,face = "bold"),
                                                                   plot.title  = element_text(size = 20)
                                                                   
                                                                   
                                                                 ),
                                                                 tables.theme = theme(text = element_text(size = 15),
                                                                                      axis.text  = element_text(size = 15),
                                                                                      axis.line = element_line(size = 0.35,colour = "black"),
                                                                                      axis.ticks = element_line(size = 0.35),
                                                                                      axis.ticks.length = unit(0.25,"mm"),
                                                                                      axis.title = element_text(size = 20),
                                                                                      plot.title  = element_text(size = 15)),
                                                                 ylab = "Recurrence survival probability",
                                                                 xlab = "months",
                                                                 data=promoter_pathway_mat_sigpw_clinical_merge_temp)
  
  
  survival_plot_rfs_109_promoter_pathway_list[[i]]$plot <- survival_plot_rfs_109_promoter_pathway_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(promoter_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(promoter_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(promoter_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
    # 
  
}

for (i in names(survival_plot_rfs_109_promoter_pathway_list)) {
  pdf(paste0("Figure_3_promoter_RFS_109_pathway_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_rfs_109_promoter_pathway_list[[i]],newpage = FALSE)
  dev.off()
}



#Survival analysis of OS in exon genes of the TCGA cohort


multi_data_exon_tcga <- na.omit(TCGA_Clinical.tidy.lusc.maf[,c(TCGA_exon_gene_sigpw,"OS.months.","Survival.status")])

survival_plot_os_tcga_list <- list()

for (i in TCGA_exon_gene_sigpw) {
  
  multi_data_exon_tcga_temp <- multi_data_exon_tcga[,c("OS.months.","Survival.status",i)]
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_tcga_temp[,i] <- factor(multi_data_exon_tcga_temp[,i],levels = unique(multi_data_exon_tcga_temp[,i]))
  multi_data_exon_tcga_temp[,i] <- relevel(multi_data_exon_tcga_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_exon_tcga_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_exon_tcga_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_os_tcga_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                palette = palette_colors,
                                                legend.title = "",
                                                size=0.6,
                                                censor.shape= "|",
                                                censor.size= 3,
                                                ggtheme = theme_classic()+ theme(
                                                  axis.line = element_line(size = 0.35,colour = "black"),
                                                  axis.ticks = element_line(size = 0.35),
                                                  axis.ticks.length = unit(0.25, "cm"),
                                                  axis.text  = element_text(size = 15,face = "bold"),
                                                  axis.title = element_text(size = 20),
                                                  axis.title.x = element_text(size = 18),
                                                  axis.title.y = element_text(size = 20),
                                                  legend.direction  = "horizontal",
                                                  legend.background = element_blank(),
                                                  legend.key.size = unit(10,"mm"),
                                                  legend.text =element_text(size = 14),
                                                  legend.title = element_text(size = 15,face = "bold"),
                                                  plot.title  = element_text(size = 20)
                                                  
                                                  
                                                ),
                                                tables.theme = theme(text = element_text(size = 15),
                                                                     axis.text  = element_text(size = 15),
                                                                     axis.line = element_line(size = 0.35,colour = "black"),
                                                                     axis.ticks = element_line(size = 0.35),
                                                                     axis.ticks.length = unit(0.25,"mm"),
                                                                     axis.title = element_text(size = 20),
                                                                     plot.title  = element_text(size = 15)),
                                                ylab = "Overall survival probability",
                                                xlab = "months",
                                                data=multi_data_exon_tcga_temp)
  
  
  survival_plot_os_tcga_list[[i]]$plot <- survival_plot_os_tcga_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(multi_data_exon_tcga_temp[,"OS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(multi_data_exon_tcga_temp[,"OS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(multi_data_exon_tcga_temp[,"OS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
    # 
  
}


for (i in names(survival_plot_os_tcga_list)) {
  pdf(paste0("Figure_3_exon_OS_tcga_",i,".pdf"),width = 5,height = 5.5)
  print(survival_plot_os_tcga_list[[i]],newpage = FALSE)
  dev.off()
}




#Survival analysis of RFS in exon genes of the TCGA cohort



multi_data_exon_tcga <- na.omit(TCGA_Clinical.tidy.lusc.maf[,c(TCGA_exon_gene_sigpw,"RFS.months.","Recurrence.status")])


survival_plot_rfs_tcga_list <- list()

for (i in TCGA_exon_gene_sigpw) {
  
  multi_data_exon_tcga_temp <- multi_data_exon_tcga[,c("RFS.months.","Recurrence.status",i)]
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_tcga_temp[,i] <- factor(multi_data_exon_tcga_temp[,i],levels = unique(multi_data_exon_tcga_temp[,i]))
  multi_data_exon_tcga_temp[,i] <- relevel(multi_data_exon_tcga_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_exon_tcga_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_exon_tcga_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_rfs_tcga_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                 palette = palette_colors,
                                                 legend.title = "",
                                                 size=0.6,
                                                 censor.shape= "|",
                                                 censor.size= 3,
                                                 ggtheme = theme_classic()+ theme(
                                                   axis.line = element_line(size = 0.35,colour = "black"),
                                                   axis.ticks = element_line(size = 0.35),
                                                   axis.ticks.length = unit(0.25, "cm"),
                                                   axis.text  = element_text(size = 15,face = "bold"),
                                                   axis.title = element_text(size = 20),
                                                   axis.title.x = element_text(size = 18),
                                                   axis.title.y = element_text(size = 20),
                                                   legend.direction  = "horizontal",
                                                   legend.background = element_blank(),
                                                   legend.key.size = unit(10,"mm"),
                                                   legend.text =element_text(size = 14),
                                                   legend.title = element_text(size = 15,face = "bold"),
                                                   plot.title  = element_text(size = 20)
                                                   
                                                   
                                                 ),
                                                 tables.theme = theme(text = element_text(size = 15),
                                                                      axis.text  = element_text(size = 15),
                                                                      axis.line = element_line(size = 0.35,colour = "black"),
                                                                      axis.ticks = element_line(size = 0.35),
                                                                      axis.ticks.length = unit(0.25,"mm"),
                                                                      axis.title = element_text(size = 20),
                                                                      plot.title  = element_text(size = 15)),
                                                 ylab = "Recurrence survival probability",
                                                 xlab = "months",
                                                 data=multi_data_exon_tcga_temp)
  
  
  survival_plot_rfs_tcga_list[[i]]$plot <- survival_plot_rfs_tcga_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(multi_data_exon_tcga_temp[,"RFS.months."])/100, y = 0.15, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(multi_data_exon_tcga_temp[,"RFS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(multi_data_exon_tcga_temp[,"RFS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
  
  
}


for (i in names(survival_plot_rfs_tcga_list)) {
  pdf(paste0("Figure_3_exon_RFS_tcga_",i,".pdf"),width = 5,height = 5)
  print(survival_plot_rfs_tcga_list[[i]],newpage = FALSE)
  dev.off()
}


#Survival analysis of OS in exon pathways of the TCGA cohort
survival_plot_os_tcga_pathway_list <- list()

for (i in make.names(colnames(tcga_exon_pathway_mat_sigpw_clinical)[2:ncol(tcga_exon_pathway_mat_sigpw_clinical)])) {
  
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp <- na.omit(tcga_exon_pathway_mat_sigpw_clinical_merge[,c("OS.months.","Survival.status",i)])
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = tcga_exon_pathway_mat_sigpw_clinical_merge_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = tcga_exon_pathway_mat_sigpw_clinical_merge_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_os_tcga_pathway_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                        palette = palette_colors,
                                                        legend.title = "",
                                                        size=0.6,
                                                        censor.shape= "|",
                                                        censor.size= 3,
                                                        ggtheme = theme_classic()+ theme(
                                                          axis.line = element_line(size = 0.35,colour = "black"),
                                                          axis.ticks = element_line(size = 0.35),
                                                          axis.ticks.length = unit(0.25, "cm"),
                                                          axis.text  = element_text(size = 15,face = "bold"),
                                                          axis.title = element_text(size = 20),
                                                          axis.title.x = element_text(size = 18),
                                                          axis.title.y = element_text(size = 20),
                                                          legend.direction  = "horizontal",
                                                          legend.background = element_blank(),
                                                          legend.key.size = unit(10,"mm"),
                                                          legend.text =element_text(size = 14),
                                                          legend.title = element_text(size = 15,face = "bold"),
                                                          plot.title  = element_text(size = 20)
                                                          
                                                          
                                                        ),
                                                        tables.theme = theme(text = element_text(size = 15),
                                                                             axis.text  = element_text(size = 15),
                                                                             axis.line = element_line(size = 0.35,colour = "black"),
                                                                             axis.ticks = element_line(size = 0.35),
                                                                             axis.ticks.length = unit(0.25,"mm"),
                                                                             axis.title = element_text(size = 20),
                                                                             plot.title  = element_text(size = 15)),
                                                        ylab = "Overall survival probability",
                                                        xlab = "months",
                                                        data=tcga_exon_pathway_mat_sigpw_clinical_merge_temp)
  
  
  survival_plot_os_tcga_pathway_list[[i]]$plot <- survival_plot_os_tcga_pathway_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.30, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,"OS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
  
  
}


for (i in names(survival_plot_os_tcga_pathway_list)) {
  pdf(paste0("Figure_3_exon_OS_tcga_pathway_",i,".pdf"),width = 5,height = 5)
  print(survival_plot_os_tcga_pathway_list[[i]],newpage = FALSE)
  dev.off()
}



#Survival analysis of RFS in exon pathways of the TCGA cohort

survival_plot_rfs_tcga_pathway_list <- list()

for (i in make.names(colnames(tcga_exon_pathway_mat_sigpw_clinical)[2:ncol(tcga_exon_pathway_mat_sigpw_clinical)])) {
  
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp <- na.omit(tcga_exon_pathway_mat_sigpw_clinical_merge[,c("RFS.months.","Recurrence.status",i)])
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  
  
  coxph_model <- coxph(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = tcga_exon_pathway_mat_sigpw_clinical_merge_temp, na.action = na.omit)
  
  coxph_summary <- summary(coxph_model)
  pvalue <- coxph_summary$coefficients[, 5]
  HR <- coxph_summary$coefficients[, 2]
  up95 <- coxph_summary$conf.int[, 4]
  low95 <- coxph_summary$conf.int[, 3]
  fit_temp <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = tcga_exon_pathway_mat_sigpw_clinical_merge_temp)
  names(fit_temp$strata) <- sub(".*=", "", names(fit_temp$strata))
  palette_colors <- c("red","blue")
  names(palette_colors) <- c("MUT","WT")
  survival_plot_rfs_tcga_pathway_list[[i]] <- ggsurvplot(fit_temp, pval = FALSE, risk.table = TRUE, risk.table.y.text.col = TRUE,
                                                         palette = palette_colors,
                                                         legend.title = "",
                                                         size=0.6,
                                                         censor.shape= "|",
                                                         censor.size= 3,
                                                         ggtheme = theme_classic()+ theme(
                                                           axis.line = element_line(size = 0.35,colour = "black"),
                                                           axis.ticks = element_line(size = 0.35),
                                                           axis.ticks.length = unit(0.25, "cm"),
                                                           axis.text  = element_text(size = 15,face = "bold"),
                                                           axis.title = element_text(size = 20),
                                                           axis.title.x = element_text(size = 18),
                                                           axis.title.y = element_text(size = 20),
                                                           legend.direction  = "horizontal",
                                                           legend.background = element_blank(),
                                                           legend.key.size = unit(10,"mm"),
                                                           legend.text =element_text(size = 14),
                                                           legend.title = element_text(size = 15,face = "bold"),
                                                           plot.title  = element_text(size = 20)
                                                           
                                                           
                                                         ),
                                                         tables.theme = theme(text = element_text(size = 15),
                                                                              axis.text  = element_text(size = 15),
                                                                              axis.line = element_line(size = 0.35,colour = "black"),
                                                                              axis.ticks = element_line(size = 0.35),
                                                                              axis.ticks.length = unit(0.25,"mm"),
                                                                              axis.title = element_text(size = 20),
                                                                              plot.title  = element_text(size = 15)),
                                                         ylab = "Recurrence survival probability",
                                                         xlab = "months",
                                                         data=tcga_exon_pathway_mat_sigpw_clinical_merge_temp)
  
  
  survival_plot_rfs_tcga_pathway_list[[i]]$plot <- survival_plot_rfs_tcga_pathway_list[[i]]$plot +
    ggtitle(i) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggplot2::annotate("text", x = max(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.30, label = paste("p = ", round(pvalue, 3)), color = "black", size = 6, hjust = 0) #+
    # ggplot2::annotate("text", x = max(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.20, label = paste("HR = ", round(HR, 3)), color = "black", size = 4.4, hjust = 0) +
    # ggplot2::annotate("text", x = max(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,"RFS.months."])/100, y = 0.10, label = paste("(", "95%CI:", round(low95, 3), "-", round(up95, 3), ")", sep = ""), color = "black", size = 4.4, hjust = 0)
    # 
  
}


for (i in names(survival_plot_rfs_tcga_pathway_list)) {
  pdf(paste0("Figure_3_exon_RFS_tcga_pathway_",i,".pdf"),width = 5,height = 5)
  print(survival_plot_rfs_tcga_pathway_list[[i]],newpage = FALSE)
  dev.off()
}



#table1


clinical_109_table1 <- list_109[[1]][,c(1:11,which(colnames(list_109[[1]])%in%c("FAT3","TLE1","SMAD2","WNT10B")))]
colnames(clinical_109_table1)[1:11] <- c("Patient_ID","Gender","Age","Smoking_history","Differentiation","Tumor_size","Stage",
                                         "Recurrence_status","RFS(months)","Survival_status","OS(months)")

colnames(clinical_109_table1)[which(colnames(clinical_109_table1) %in% "FAT3")] <- "FAT3_promoter"
colnames(clinical_109_table1)[which(colnames(clinical_109_table1) %in% "TLE1")] <- "TLE1_promoter"
colnames(clinical_109_table1)[which(colnames(clinical_109_table1) %in% "SMAD2")] <- "SMAD2_promoter"
colnames(clinical_109_table1)[which(colnames(clinical_109_table1) %in% "WNT10B")] <- "WNT10B_promoter"

clinical_109_table1[,"FAT4_exon"] <- list_109_exon[[1]][rownames(clinical_109_table1),"FAT4"]
clinical_109_table1[,"DCHS1_exon"] <- list_109_exon[[1]][rownames(clinical_109_table1),"DCHS1"]

clinical_109_table1[,"TGF.Beta"] <- promoter_pathway_mat_sigpw_clinical_merge[rownames(clinical_109_table1),"TGF.Beta"]
clinical_109_table1[,"Hippo"] <- promoter_pathway_mat_sigpw_clinical_merge[rownames(clinical_109_table1),"Hippo"]




clinical_109_table1[clinical_109_table1$Smoking_history %in% 1,"Smoking_history"] <- "smoker"
clinical_109_table1[clinical_109_table1$Smoking_history %in% 0,"Smoking_history"] <- "non-smoker"

clinical_109_table1[clinical_109_table1$Differentiation %in% 1,"Differentiation"] <- "high"
clinical_109_table1[clinical_109_table1$Differentiation %in% 2,"Differentiation"] <- "moderate or low"

clinical_109_table1[clinical_109_table1$Tumor_size %in% 1,"Tumor_size"] <- "> 3cm"
clinical_109_table1[clinical_109_table1$Tumor_size %in% 0,"Tumor_size"] <- "<= 3cm"


clinical_109_table1[clinical_109_table1$Survival_status %in% 1,"Survival_status"] <- "Dead"
clinical_109_table1[clinical_109_table1$Survival_status %in% 0,"Survival_status"] <- "Alive"

clinical_109_table1[clinical_109_table1$Recurrence_status %in% 1,"Recurrence_status"] <- "Recurrence"
clinical_109_table1[clinical_109_table1$Recurrence_status %in% 0,"Recurrence_status"] <- "No Recurrence"


clinical_109_table1[clinical_109_table1$SMAD2_promoter %in% 1,"SMAD2_promoter"] <- "SMAD2_promoter_MUT"
clinical_109_table1[clinical_109_table1$SMAD2_promoter %in% 0,"SMAD2_promoter"] <- "SMAD2_promoter_WT"

clinical_109_table1[clinical_109_table1$FAT3_promoter %in% 1,"FAT3_promoter"] <- "FAT3_promoter_MUT"
clinical_109_table1[clinical_109_table1$FAT3_promoter %in% 0,"FAT3_promoter"] <- "FAT3_promoter_WT"


clinical_109_table1[clinical_109_table1$TLE1_promoter %in% 1,"TLE1_promoter"] <- "TLE1_promoter_MUT"
clinical_109_table1[clinical_109_table1$TLE1_promoter %in% 0,"TLE1_promoter"] <- "TLE1_promoter_WT"

clinical_109_table1[clinical_109_table1$WNT10B_promoter %in% 1,"WNT10B_promoter"] <- "WNT10B_promoter_MUT"
clinical_109_table1[clinical_109_table1$WNT10B_promoter %in% 0,"WNT10B_promoter"] <- "WNT10B_promoter_WT"

colnames(clinical_109_table1) <- make.names(colnames(clinical_109_table1))


table1_all <- table1::table1(~ Gender + Age + Smoking_history + Differentiation +
                 Tumor_size + Stage + Recurrence_status + Survival_status
               ,data=clinical_109_table1)

library(xlsx)



table1_FAT3_promoter <- table1::table1(~ Gender + Age + Smoking_history + Differentiation +
                 Tumor_size + Stage + Recurrence_status + Survival_status | FAT3_promoter ,
               extra.col	= list("P-value"=pvalue),data=clinical_109_table1)


table1_SMAD2_promoter <- table1::table1(~ Gender + Age + Smoking_history + Differentiation +
                 Tumor_size + Stage + Recurrence_status + Survival_status | SMAD2_promoter ,
               extra.col	= list("P-value"=pvalue),data=clinical_109_table1)



table1_TLE1_promoter <- table1::table1(~ Gender + Age + Smoking_history + Differentiation +
                 Tumor_size + Stage + Recurrence_status + Survival_status | TLE1_promoter ,
               extra.col	= list("P-value"=pvalue),data=clinical_109_table1)



table1_WNT10B_promoter <- table1::table1(~ Gender + Age + Smoking_history + Differentiation +
                 Tumor_size + Stage + Recurrence_status + Survival_status | WNT10B_promoter ,
               extra.col	= list("P-value"=pvalue),data=clinical_109_table1)




multicox_109_pathway_exon_genes2 <- multicox_109_pathway_exon_genes
multicox_109_pathway_promoter_genes2 <- multicox_109_pathway_promoter_genes
exon_pathway_multicox_merge_survial2 <- exon_pathway_multicox_merge_survial
promoter_pathway_multicox_merge_survial2 <- promoter_pathway_multicox_merge_survial
tcga_exon_gene_multicox_rfs_os_merge_survial2 <- tcga_exon_gene_multicox_rfs_os_merge_survial
exon_tcga_pathway_multicox_os_rfs_merge_survial2 <- exon_tcga_pathway_multicox_os_rfs_merge_survial

multicox_109_pathway_exon_genes2[,2:9] <- lapply(multicox_109_pathway_exon_genes2[,2:9],signif,digits = 3)
multicox_109_pathway_promoter_genes2[,2:9] <- lapply(multicox_109_pathway_promoter_genes2[,2:9],signif,digits = 3)
exon_pathway_multicox_merge_survial2[,2:9] <- lapply(exon_pathway_multicox_merge_survial2[,2:9],signif,digits = 3)
promoter_pathway_multicox_merge_survial2[,2:9] <- lapply(promoter_pathway_multicox_merge_survial2[,2:9],signif,digits = 3)
tcga_exon_gene_multicox_rfs_os_merge_survial2[,2:9] <- lapply(tcga_exon_gene_multicox_rfs_os_merge_survial2[,2:9],signif,digits = 3)
exon_tcga_pathway_multicox_os_rfs_merge_survial2[,2:9] <- lapply(exon_tcga_pathway_multicox_os_rfs_merge_survial2[,2:9],signif,digits = 3)



clinical_109_table1 <- list_109[[1]][,c(1:11,which(colnames(list_109[[1]])%in%multicox_109_pathway_promoter_genes$Gene))]
colnames(clinical_109_table1)[1:11] <- c("Patient_ID","Gender","Age","Smoking_history","Differentiation","Tumor_size","Stage",
                                         "Recurrence_status","RFS(months)","Survival_status","OS(months)")

colnames(clinical_109_table1)[12:ncol(clinical_109_table1)] <- paste0(colnames(clinical_109_table1)[12:ncol(clinical_109_table1)],"_promoter")
clinical_109_table1[,ncol(clinical_109_table1)+1:length(multicox_109_pathway_exon_genes$Gene)] <- list_109_exon[[1]][rownames(clinical_109_table1),multicox_109_pathway_exon_genes$Gene]

colnames(clinical_109_table1)[(ncol(clinical_109_table1)-length(multicox_109_pathway_exon_genes$Gene)+1):ncol(clinical_109_table1)] <- 
  paste0(colnames(clinical_109_table1)[(ncol(clinical_109_table1)-length(multicox_109_pathway_exon_genes$Gene)+1):ncol(clinical_109_table1)],"_exon")

#
clinical_109_table1[,ncol(clinical_109_table1)+1:length(promoter_pathway_multicox_merge_survial2$Gene)] <- promoter_pathway_mat_sigpw_clinical_merge[rownames(clinical_109_table1),promoter_pathway_multicox_merge_survial2$Gene]

colnames(clinical_109_table1)[(ncol(clinical_109_table1)-length(promoter_pathway_multicox_merge_survial2$Gene)+1):ncol(clinical_109_table1)] <- 
  paste0(colnames(clinical_109_table1)[(ncol(clinical_109_table1)-length(promoter_pathway_multicox_merge_survial2$Gene)+1):ncol(clinical_109_table1)],"_pathway_promoter")

clinical_109_table1[,ncol(clinical_109_table1)+1:length(exon_pathway_multicox_merge_survial2$Gene)] <- exon_pathway_mat_sigpw_clinical_merge[rownames(clinical_109_table1),exon_pathway_multicox_merge_survial2$Gene]

colnames(clinical_109_table1)[(ncol(clinical_109_table1)-length(exon_pathway_multicox_merge_survial2$Gene)+1):ncol(clinical_109_table1)] <- 
  paste0(colnames(clinical_109_table1)[(ncol(clinical_109_table1)-length(exon_pathway_multicox_merge_survial2$Gene)+1):ncol(clinical_109_table1)],"_pathway_exon")


# clinical_109_table1[,"TGF_Beta_promoter"] <- promoter_pathway_mat_sigpw_clinical_merge[rownames(clinical_109_table1),"TGF.Beta"]
# clinical_109_table1[,"Hippo_promoter"] <- promoter_pathway_mat_sigpw_clinical_merge[rownames(clinical_109_table1),"Hippo"]
clinical_109_table1[,c(12:ncol(clinical_109_table1))][clinical_109_table1[,c(12:ncol(clinical_109_table1))] == 0] <- "WT"

clinical_109_table1[,c(12:ncol(clinical_109_table1))][clinical_109_table1[,c(12:ncol(clinical_109_table1))] == 1] <- "MUT"

clinical_109_table1[,c(12:ncol(clinical_109_table1))] <- lapply(clinical_109_table1[,c(12:ncol(clinical_109_table1))], factor)

clinical_109_table1 <- clinical_109_table1[,c(9,11,c(12:ncol(clinical_109_table1)))]

tableOne_list <- list()
tableOne_list2 <- list()
tableOne_list2_out <- data.frame()

for (i in colnames(clinical_109_table1)[3:ncol(clinical_109_table1)]) {
  
clinical_109_table1_2 <- clinical_109_table1[,c("RFS(months)","OS(months)",i)]

tableOne_list[[i]]  <- CreateTableOne(data = clinical_109_table1_2,
                           strata = colnames(clinical_109_table1_2)[3],
                           factorVars =  colnames(clinical_109_table1_2)[3])

tableOne_list2[[i]] <- t(print(tableOne_list[[i]],nonnormal = c("RFS(months)","OS(months)"),showAllLevels = TRUE)[c(2:3),c(2:3)])
tableOne_list2[[i]] <- data.frame(name1=rep(i,ncol(tableOne_list2[[i]])),name2=rownames(tableOne_list2[[i]]),tableOne_list2[[i]],check.names = FALSE)
rownames(tableOne_list2[[i]]) <- NULL
tableOne_list2_out <- rbind(tableOne_list2_out,tableOne_list2[[i]])

}

tableOne_list2_out <- data.frame(flag=paste0(tableOne_list2_out$name1,"_",tableOne_list2_out$name2),tableOne_list2_out,check.names = FALSE)




library(tableone)
tableOne_all  <- CreateTableOne(data = clinical_109_table1[,3:ncol(clinical_109_table1)])


rownames(print(tableOne_all,nonnormal = c("RFS(months)","OS(months)"),showAllLevels = TRUE))

tableOne_all_out <- print(tableOne_all,nonnormal = c("RFS(months)","OS(months)"),showAllLevels = TRUE)
tableOne_all_out <- tableOne_all_out[-1,]

for (i in seq(1, nrow(tableOne_all_out), by = 2)) {
  if (i + 1 <= nrow(tableOne_all_out)) {
    rownames(tableOne_all_out)[i + 1] <- rownames(tableOne_all_out)[i]
  }
}

tableOne_all_out <- data.frame(name=gsub("\ \\(%\\)","",rownames(tableOne_all_out)),tableOne_all_out)
rownames(tableOne_all_out) <- NULL
tableOne_all_out <- data.frame(flag=paste0(tableOne_all_out$name,"_",tableOne_all_out$level),tableOne_all_out)

tableOne_out <- merge(tableOne_all_out[,c("flag","Overall")],tableOne_list2_out,by="flag")




multicox_109_pathway_exon_genes2 <- do.call(rbind, lapply(1:nrow(multicox_109_pathway_exon_genes2), function(i) {
  multicox_109_pathway_exon_genes2 <- rbind(multicox_109_pathway_exon_genes2[i, ], multicox_109_pathway_exon_genes2[i, ])
  multicox_109_pathway_exon_genes2 <- data.frame(flag=paste0(multicox_109_pathway_exon_genes2$Gene,"_exon"),
                                                 OS_CI = paste0(multicox_109_pathway_exon_genes2$OS_HR," (",
                                                                multicox_109_pathway_exon_genes2$OS_lower_95," - ",
                                                                multicox_109_pathway_exon_genes2$OS_upper_95,
                                                                ")"),
                                                 OS_P_value=multicox_109_pathway_exon_genes2$OS_P_value,
                                                 RFS_CI = paste0(multicox_109_pathway_exon_genes2$RFS_HR," (",
                                                                 multicox_109_pathway_exon_genes2$RFS_lower_95," - ",
                                                                 multicox_109_pathway_exon_genes2$RFS_upper_95,
                                                                 ")"),
                                                 RFS_P_value=multicox_109_pathway_exon_genes2$RFS_P_value
  )
  multicox_109_pathway_exon_genes2[1, 1] <- paste0(multicox_109_pathway_exon_genes2[1, 1], "_WT")  # 第一列添加"_WT"
  multicox_109_pathway_exon_genes2[1, c(2,3,4,5)] <- c("1 [reference]",NA,"1 [reference]",NA)
  multicox_109_pathway_exon_genes2[2, 1] <- paste0(multicox_109_pathway_exon_genes2[2, 1], "_MUT")  # 第一列添加"_WT"
  
  return(multicox_109_pathway_exon_genes2)
  
  
}))

multicox_109_pathway_promoter_genes2 <- do.call(rbind, lapply(1:nrow(multicox_109_pathway_promoter_genes2), function(i) {
  multicox_109_pathway_promoter_genes2 <- rbind(multicox_109_pathway_promoter_genes2[i, ], multicox_109_pathway_promoter_genes2[i, ])
  multicox_109_pathway_promoter_genes2 <- data.frame(flag=paste0(multicox_109_pathway_promoter_genes2$Gene,"_promoter"),
                                                     OS_CI = paste0(multicox_109_pathway_promoter_genes2$OS_HR," (",
                                                                    multicox_109_pathway_promoter_genes2$OS_lower_95," - ",
                                                                    multicox_109_pathway_promoter_genes2$OS_upper_95,
                                                                    ")"),
                                                     OS_P_value=multicox_109_pathway_promoter_genes2$OS_P_value,
                                                     RFS_CI = paste0(multicox_109_pathway_promoter_genes2$RFS_HR," (",
                                                                     multicox_109_pathway_promoter_genes2$RFS_lower_95," - ",
                                                                     multicox_109_pathway_promoter_genes2$RFS_upper_95,
                                                                     ")"),
                                                     RFS_P_value=multicox_109_pathway_promoter_genes2$RFS_P_value
  )
  multicox_109_pathway_promoter_genes2[1, 1] <- paste0(multicox_109_pathway_promoter_genes2[1, 1], "_WT")  # 第一列添加"_WT"
  multicox_109_pathway_promoter_genes2[1, c(2,3,4,5)] <- c("1 [reference]",NA,"1 [reference]",NA)
  multicox_109_pathway_promoter_genes2[2, 1] <- paste0(multicox_109_pathway_promoter_genes2[2, 1], "_MUT")  # 第一列添加"_WT"
  
  return(multicox_109_pathway_promoter_genes2)
  
  
}))

exon_pathway_multicox_merge_survial2 <- do.call(rbind, lapply(1:nrow(exon_pathway_multicox_merge_survial2), function(i) {
  exon_pathway_multicox_merge_survial2 <- rbind(exon_pathway_multicox_merge_survial2[i, ], exon_pathway_multicox_merge_survial2[i, ])
  exon_pathway_multicox_merge_survial2 <- data.frame(flag=paste0(exon_pathway_multicox_merge_survial2$Gene,"_pathway_exon"),
                                                     OS_CI = paste0(exon_pathway_multicox_merge_survial2$OS_HR," (",
                                                                    exon_pathway_multicox_merge_survial2$OS_lower_95," - ",
                                                                    exon_pathway_multicox_merge_survial2$OS_upper_95,
                                                                    ")"),
                                                     OS_P_value=exon_pathway_multicox_merge_survial2$OS_P_value,
                                                     RFS_CI = paste0(exon_pathway_multicox_merge_survial2$RFS_HR," (",
                                                                     exon_pathway_multicox_merge_survial2$RFS_lower_95," - ",
                                                                     exon_pathway_multicox_merge_survial2$RFS_upper_95,
                                                                     ")"),
                                                     RFS_P_value=exon_pathway_multicox_merge_survial2$RFS_P_value
  )
  exon_pathway_multicox_merge_survial2[1, 1] <- paste0(exon_pathway_multicox_merge_survial2[1, 1], "_WT")  # 第一列添加"_WT"
  exon_pathway_multicox_merge_survial2[1, c(2,3,4,5)] <- c("1 [reference]",NA,"1 [reference]",NA)
  exon_pathway_multicox_merge_survial2[2, 1] <- paste0(exon_pathway_multicox_merge_survial2[2, 1], "_MUT")  # 第一列添加"_WT"
  
  return(exon_pathway_multicox_merge_survial2)
  
  
}))

promoter_pathway_multicox_merge_survial2 <- do.call(rbind, lapply(1:nrow(promoter_pathway_multicox_merge_survial2), function(i) {
  promoter_pathway_multicox_merge_survial2 <- rbind(promoter_pathway_multicox_merge_survial2[i, ], promoter_pathway_multicox_merge_survial2[i, ])
  promoter_pathway_multicox_merge_survial2 <- data.frame(flag=paste0(promoter_pathway_multicox_merge_survial2$Gene,"_pathway_promoter"),
                                                         OS_CI = paste0(promoter_pathway_multicox_merge_survial2$OS_HR," (",
                                                                        promoter_pathway_multicox_merge_survial2$OS_lower_95," - ",
                                                                        promoter_pathway_multicox_merge_survial2$OS_upper_95,
                                                                        ")"),
                                                         OS_P_value=promoter_pathway_multicox_merge_survial2$OS_P_value,
                                                         RFS_CI = paste0(promoter_pathway_multicox_merge_survial2$RFS_HR," (",
                                                                         promoter_pathway_multicox_merge_survial2$RFS_lower_95," - ",
                                                                         promoter_pathway_multicox_merge_survial2$RFS_upper_95,
                                                                         ")"),
                                                         RFS_P_value=promoter_pathway_multicox_merge_survial2$RFS_P_value
  )
  promoter_pathway_multicox_merge_survial2[1, 1] <- paste0(promoter_pathway_multicox_merge_survial2[1, 1], "_WT")  # 第一列添加"_WT"
  promoter_pathway_multicox_merge_survial2[1, c(2,3,4,5)] <- c("1 [reference]",NA,"1 [reference]",NA)
  promoter_pathway_multicox_merge_survial2[2, 1] <- paste0(promoter_pathway_multicox_merge_survial2[2, 1], "_MUT")  # 第一列添加"_WT"
  
  return(promoter_pathway_multicox_merge_survial2)
  
  
}))


tcga_exon_gene_multicox_rfs_os_merge_survial2 <- do.call(rbind, lapply(1:nrow(tcga_exon_gene_multicox_rfs_os_merge_survial2), function(i) {
  tcga_exon_gene_multicox_rfs_os_merge_survial2 <- rbind(tcga_exon_gene_multicox_rfs_os_merge_survial2[i, ], tcga_exon_gene_multicox_rfs_os_merge_survial2[i, ])
  tcga_exon_gene_multicox_rfs_os_merge_survial2 <- data.frame(flag=paste0(tcga_exon_gene_multicox_rfs_os_merge_survial2$Gene,"_exon"),
                                                              OS_CI = paste0(tcga_exon_gene_multicox_rfs_os_merge_survial2$OS_HR," (",
                                                                             tcga_exon_gene_multicox_rfs_os_merge_survial2$OS_lower_95," - ",
                                                                             tcga_exon_gene_multicox_rfs_os_merge_survial2$OS_upper_95,
                                                                             ")"),
                                                              OS_P_value=tcga_exon_gene_multicox_rfs_os_merge_survial2$OS_P_value,
                                                              RFS_CI = paste0(tcga_exon_gene_multicox_rfs_os_merge_survial2$RFS_HR," (",
                                                                              tcga_exon_gene_multicox_rfs_os_merge_survial2$RFS_lower_95," - ",
                                                                              tcga_exon_gene_multicox_rfs_os_merge_survial2$RFS_upper_95,
                                                                              ")"),
                                                              RFS_P_value=tcga_exon_gene_multicox_rfs_os_merge_survial2$RFS_P_value
  )
  tcga_exon_gene_multicox_rfs_os_merge_survial2[1, 1] <- paste0(tcga_exon_gene_multicox_rfs_os_merge_survial2[1, 1], "_WT")  # 第一列添加"_WT"
  tcga_exon_gene_multicox_rfs_os_merge_survial2[1, c(2,3,4,5)] <- c("1 [reference]",NA,"1 [reference]",NA)
  tcga_exon_gene_multicox_rfs_os_merge_survial2[2, 1] <- paste0(tcga_exon_gene_multicox_rfs_os_merge_survial2[2, 1], "_MUT")  # 第一列添加"_WT"
  
  return(tcga_exon_gene_multicox_rfs_os_merge_survial2)
  
  
}))


exon_tcga_pathway_multicox_os_rfs_merge_survial2 <- do.call(rbind, lapply(1:nrow(exon_tcga_pathway_multicox_os_rfs_merge_survial2), function(i) {
  exon_tcga_pathway_multicox_os_rfs_merge_survial2 <- rbind(exon_tcga_pathway_multicox_os_rfs_merge_survial2[i, ], exon_tcga_pathway_multicox_os_rfs_merge_survial2[i, ])
  exon_tcga_pathway_multicox_os_rfs_merge_survial2 <- data.frame(flag=paste0(exon_tcga_pathway_multicox_os_rfs_merge_survial2$Gene,"_pathway_exon"),
                                                                 OS_CI = paste0(exon_tcga_pathway_multicox_os_rfs_merge_survial2$OS_HR," (",
                                                                                exon_tcga_pathway_multicox_os_rfs_merge_survial2$OS_lower_95," - ",
                                                                                exon_tcga_pathway_multicox_os_rfs_merge_survial2$OS_upper_95,
                                                                                ")"),
                                                                 OS_P_value=exon_tcga_pathway_multicox_os_rfs_merge_survial2$OS_P_value,
                                                                 RFS_CI = paste0(exon_tcga_pathway_multicox_os_rfs_merge_survial2$RFS_HR," (",
                                                                                 exon_tcga_pathway_multicox_os_rfs_merge_survial2$RFS_lower_95," - ",
                                                                                 exon_tcga_pathway_multicox_os_rfs_merge_survial2$RFS_upper_95,
                                                                                 ")"),
                                                                 RFS_P_value=exon_tcga_pathway_multicox_os_rfs_merge_survial2$RFS_P_value
  )
  exon_tcga_pathway_multicox_os_rfs_merge_survial2[1, 1] <- paste0(exon_tcga_pathway_multicox_os_rfs_merge_survial2[1, 1], "_WT")  # 第一列添加"_WT"
  exon_tcga_pathway_multicox_os_rfs_merge_survial2[1, c(2,3,4,5)] <- c("1 [reference]",NA,"1 [reference]",NA)
  exon_tcga_pathway_multicox_os_rfs_merge_survial2[2, 1] <- paste0(exon_tcga_pathway_multicox_os_rfs_merge_survial2[2, 1], "_MUT")  # 第一列添加"_WT"
  
  return(exon_tcga_pathway_multicox_os_rfs_merge_survial2)
  
  
}))


tableOne_out_2 <- rbind(multicox_109_pathway_exon_genes2,multicox_109_pathway_promoter_genes2,
                                          exon_pathway_multicox_merge_survial2,promoter_pathway_multicox_merge_survial2)



fit_os_109_list <- list()
year_survival_os_109_df <- data.frame()

for (i in exon_gene_sigpw) {
  
  multi_data_exon_temp <- multi_data_exon[,c("OS.months.","Survival.status",i)]
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_temp[,i] <- factor(multi_data_exon_temp[,i],levels = unique(multi_data_exon_temp[,i]))
  multi_data_exon_temp[,i] <- relevel(multi_data_exon_temp[,i],ref="WT")
  fit_os_109_list[[i]] <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_exon_temp)
  year_survival_os_109_df <- rbind(year_survival_os_109_df,
                                   data.frame(flag=gsub("=","_exon_",as.vector(summary(fit_os_109_list[[i]],times = 12)$strata)),
                                              `OS 1 Year Survival`=  summary(fit_os_109_list[[i]],times = 12)$surv,
                                              `OS 2 Year Survival`=  summary(fit_os_109_list[[i]],times = 24)$surv,check.names = FALSE)
                                   )
  
}


fit_rfs_109_list <- list()
year_survival_rfs_109_df <- data.frame()

for (i in exon_gene_sigpw) {
  
  multi_data_exon_temp <- multi_data_exon[,c("RFS.months.", "Recurrence.status",i)]
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_temp[multi_data_exon_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_temp[,i] <- factor(multi_data_exon_temp[,i],levels = unique(multi_data_exon_temp[,i]))
  multi_data_exon_temp[,i] <- relevel(multi_data_exon_temp[,i],ref="WT")
  fit_rfs_109_list[[i]] <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_exon_temp)
  year_survival_rfs_109_df <- rbind(year_survival_rfs_109_df,
                                   data.frame(flag=gsub("=","_exon_",as.vector(summary(fit_rfs_109_list[[i]],times = 12)$strata)),
                                              `RFS 1 Year Survival`=  summary(fit_rfs_109_list[[i]],times = 12)$surv,
                                              `RFS 2 Year Survival`=  summary(fit_rfs_109_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}



fit_os_109_pathway_list <- list()
year_survival_os_109_pathway_df <- data.frame()

for (i in make.names(rownames(exon_pathway_mat_sigpw))) {
  
  exon_pathway_mat_sigpw_clinical_merge_temp <- exon_pathway_mat_sigpw_clinical_merge[,c("OS.months.","Survival.status",i)]
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  fit_os_109_pathway_list[[i]] <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = exon_pathway_mat_sigpw_clinical_merge_temp)
  year_survival_os_109_pathway_df <- rbind(year_survival_os_109_pathway_df,
                                           data.frame(flag=gsub("=","_pathway_exon_",as.vector(summary(fit_os_109_pathway_list[[i]],times = 12)$strata)),
                                                      `OS 1 Year Survival`=  summary(fit_os_109_pathway_list[[i]],times = 12)$surv,
                                                      `OS 2 Year Survival`=  summary(fit_os_109_pathway_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}


fit_rfs_109_pathway_list <- list()
year_survival_rfs_109_pathway_df <- data.frame()

for (i in make.names(rownames(exon_pathway_mat_sigpw))) {
  
  exon_pathway_mat_sigpw_clinical_merge_temp <- exon_pathway_mat_sigpw_clinical_merge[,c("RFS.months.", "Recurrence.status",i)]
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  exon_pathway_mat_sigpw_clinical_merge_temp[exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  fit_rfs_109_pathway_list[[i]] <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = exon_pathway_mat_sigpw_clinical_merge_temp)
  year_survival_rfs_109_pathway_df <- rbind(year_survival_rfs_109_pathway_df,
                                    data.frame(flag=gsub("=","_pathway_exon_",as.vector(summary(fit_rfs_109_pathway_list[[i]],times = 12)$strata)),
                                               `RFS 1 Year Survival`=  summary(fit_rfs_109_pathway_list[[i]],times = 12)$surv,
                                               `RFS 2 Year Survival`=  summary(fit_rfs_109_pathway_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}



fit_os_109_promoter_list <- list()
year_survival_os_promoter_109_df <- data.frame()

for (i in promoter_gene_sigpw) {
  
  multi_data_promoter_temp <- multi_data_promoter[,c("OS.months.","Survival.status",i)]
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%1,i] <- "MUT"
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%0,i] <- "WT"
  multi_data_promoter_temp[,i] <- factor(multi_data_promoter_temp[,i],levels = unique(multi_data_promoter_temp[,i]))
  multi_data_promoter_temp[,i] <- relevel(multi_data_promoter_temp[,i],ref="WT")
  fit_os_109_promoter_list[[i]] <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_promoter_temp)
  year_survival_os_promoter_109_df <- rbind(year_survival_os_promoter_109_df,
                                            data.frame(flag=gsub("=","_promoter_",as.vector(summary(fit_os_109_promoter_list[[i]],times = 12)$strata)),
                                                       `OS 1 Year Survival`=  summary(fit_os_109_promoter_list[[i]],times = 12)$surv,
                                                       `OS 2 Year Survival`=  summary(fit_os_109_promoter_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}


fit_rfs_109_promoter_list <- list()
year_survival_rfs_promoter_109_df <- data.frame()

for (i in promoter_gene_sigpw) {
  
  multi_data_promoter_temp <- multi_data_promoter[,c("RFS.months.", "Recurrence.status",i)]
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%1,i] <- "MUT"
  multi_data_promoter_temp[multi_data_promoter_temp[,i]%in%0,i] <- "WT"
  multi_data_promoter_temp[,i] <- factor(multi_data_promoter_temp[,i],levels = unique(multi_data_promoter_temp[,i]))
  multi_data_promoter_temp[,i] <- relevel(multi_data_promoter_temp[,i],ref="WT")
  fit_rfs_109_promoter_list[[i]] <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_promoter_temp)
  year_survival_rfs_promoter_109_df <- rbind(year_survival_rfs_promoter_109_df,
                                             data.frame(flag=gsub("=","_promoter_",as.vector(summary(fit_rfs_109_promoter_list[[i]],times = 12)$strata)),
                                                        `RFS 1 Year Survival`=  summary(fit_rfs_109_promoter_list[[i]],times = 12)$surv,
                                                        `RFS 2 Year Survival`=  summary(fit_rfs_109_promoter_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}



fit_os_109_promoter_pathway_list <- list()
year_survival_os_promoter_109_pathway_df <- data.frame()

for (i in make.names(rownames(promoter_pathway_mat_sigpw))) {
  
  promoter_pathway_mat_sigpw_clinical_merge_temp <- promoter_pathway_mat_sigpw_clinical_merge[,c("OS.months.","Survival.status",i)]
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(promoter_pathway_mat_sigpw_clinical_merge_temp[,i]))
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  fit_os_109_promoter_pathway_list[[i]] <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = promoter_pathway_mat_sigpw_clinical_merge_temp)
  year_survival_os_promoter_109_pathway_df <- rbind(year_survival_os_promoter_109_pathway_df,
                                           data.frame(flag=gsub("=","_pathway_promoter_",as.vector(summary(fit_os_109_promoter_pathway_list[[i]],times = 12)$strata)),
                                                      `OS 1 Year Survival`=  summary(fit_os_109_promoter_pathway_list[[i]],times = 12)$surv,
                                                      `OS 2 Year Survival`=  summary(fit_os_109_promoter_pathway_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}


fit_rfs_109_promoter_pathway_list <- list()
year_survival_rfs_promoter_109_pathway_df <- data.frame()

for (i in make.names(rownames(promoter_pathway_mat_sigpw))) {
  
  promoter_pathway_mat_sigpw_clinical_merge_temp <- promoter_pathway_mat_sigpw_clinical_merge[,c("RFS.months.", "Recurrence.status",i)]
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[promoter_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(promoter_pathway_mat_sigpw_clinical_merge_temp[,i]))
  promoter_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(promoter_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  fit_rfs_109_promoter_pathway_list[[i]] <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = promoter_pathway_mat_sigpw_clinical_merge_temp)
  year_survival_rfs_promoter_109_pathway_df <- rbind(year_survival_rfs_promoter_109_pathway_df,
                                                     data.frame(flag=gsub("=","_pathway_promoter_",as.vector(summary(fit_rfs_109_promoter_pathway_list[[i]],times = 12)$strata)),
                                                                `RFS 1 Year Survival`=  summary(fit_rfs_109_promoter_pathway_list[[i]],times = 12)$surv,
                                                                `RFS 2 Year Survival`=  summary(fit_rfs_109_promoter_pathway_list[[i]],times = 24)$surv,check.names = FALSE)
  )
  
}

year_survival_109_df <- rbind(merge(year_survival_os_109_df,year_survival_rfs_109_df,by="flag"),
                              merge(year_survival_os_109_pathway_df,year_survival_rfs_109_pathway_df,by="flag"),
                              merge(year_survival_os_promoter_109_df,year_survival_rfs_promoter_109_df,by="flag"),
                              merge(year_survival_os_promoter_109_pathway_df,year_survival_rfs_promoter_109_pathway_df,by="flag"))


tableOne_out_3 <- merge(tableOne_out,tableOne_out_2,by="flag")
tableOne_out_4 <- merge(tableOne_out_3,year_survival_109_df,by="flag")
tableOne_out_4$pathway <- "FALSE"
tableOne_out_4[tableOne_out_4$flag%in%year_survival_rfs_promoter_109_pathway_df$flag,"pathway"] <- "TRUE"
tableOne_out_5 <- tableOne_out_4[,c("name1","name2","Overall","RFS(months) (median [IQR])",
                                    "RFS 1 Year Survival","RFS 2 Year Survival","RFS_CI",
                                    "RFS_P_value","Overall","OS(months) (median [IQR])",
                                    "OS 1 Year Survival","OS 2 Year Survival",
                                    "OS_CI","OS_P_value","pathway")]

tableOne_out_5$`RFS 1 Year Survival` <- signif(tableOne_out_5$`RFS 1 Year Survival`,digits = 3)
tableOne_out_5$`RFS 2 Year Survival` <- signif(tableOne_out_5$`RFS 2 Year Survival`,digits = 3)

tableOne_out_5$`OS 1 Year Survival` <- signif(tableOne_out_5$`OS 1 Year Survival`,digits = 3)
tableOne_out_5$`OS 2 Year Survival` <- signif(tableOne_out_5$`OS 2 Year Survival`,digits = 3)

tableOne_out_5$Overall <- gsub("\\(","\r\\(",tableOne_out_5$Overall)

empty_row <- data.frame()


tableOne_out_6 <- do.call(rbind, lapply(1:nrow(tableOne_out_5), function(i) {
  if (i %% 2 == 1) {
    return(rbind(c(tableOne_out_5[i,1],rep(NA,14)), tableOne_out_5[i, , drop = FALSE]))
  } else {
    return(tableOne_out_5[i, , drop = FALSE])
  }
}))


multi_data_exon_tcga <- na.omit(TCGA_Clinical.tidy.lusc.maf[,c(TCGA_exon_gene_sigpw,"OS.months.","Survival.status")])

fit_os_tcga_exon_list <- list()
year_survival_os_exon_tcga_df <- data.frame()

for (i in TCGA_exon_gene_sigpw) {
  
  multi_data_exon_tcga_temp <- multi_data_exon_tcga[,c("OS.months.","Survival.status",i)]
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_tcga_temp[,i] <- factor(multi_data_exon_tcga_temp[,i],levels = unique(multi_data_exon_tcga_temp[,i]))
  multi_data_exon_tcga_temp[,i] <- relevel(multi_data_exon_tcga_temp[,i],ref="WT")
  fit_os_tcga_exon_list[[i]] <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = multi_data_exon_tcga_temp)
  year_survival_os_exon_tcga_df <- rbind(year_survival_os_exon_tcga_df,
                                         data.frame(flag=gsub("=","_exon_",as.vector(summary(fit_os_tcga_exon_list[[i]],times = 12)$strata)),
                                                    `OS 1 Year Survival`=  summary(fit_os_tcga_exon_list[[i]],times = 12,extend = TRUE)$surv,
                                                    `OS 2 Year Survival`=  summary(fit_os_tcga_exon_list[[i]],times = 24,extend = TRUE)$surv,
                                                    check.names = FALSE)
  )
  
}

multi_data_exon_tcga <- na.omit(TCGA_Clinical.tidy.lusc.maf[,c(TCGA_exon_gene_sigpw,"RFS.months.", "Recurrence.status")])

fit_rfs_tcga_exon_list <- list()
year_survival_rfs_exon_tcga_df <- data.frame()

for (i in TCGA_exon_gene_sigpw) {
  
  multi_data_exon_tcga_temp <- multi_data_exon_tcga[,c("RFS.months.", "Recurrence.status",i)]
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%1,i] <- "MUT"
  multi_data_exon_tcga_temp[multi_data_exon_tcga_temp[,i]%in%0,i] <- "WT"
  multi_data_exon_tcga_temp[,i] <- factor(multi_data_exon_tcga_temp[,i],levels = unique(multi_data_exon_tcga_temp[,i]))
  multi_data_exon_tcga_temp[,i] <- relevel(multi_data_exon_tcga_temp[,i],ref="WT")
  fit_rfs_tcga_exon_list[[i]] <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = multi_data_exon_tcga_temp)
  year_survival_rfs_exon_tcga_df <- rbind(year_survival_rfs_exon_tcga_df,
                                          data.frame(flag=gsub("=","_exon_",as.vector(summary(fit_rfs_tcga_exon_list[[i]],times = 12)$strata)),
                                                     `RFS 1 Year Survival`=  summary(fit_rfs_tcga_exon_list[[i]],times = 12,extend = TRUE)$surv,
                                                     `RFS 2 Year Survival`=  summary(fit_rfs_tcga_exon_list[[i]],times = 24,extend = TRUE)$surv
                                                     ,check.names = FALSE
                                                     )
  )
  
}



fit_os_tcga_exon_pathway_list <- list()
year_survival_os_exon_tcga_pathway_df <- data.frame()

for (i in make.names(colnames(tcga_exon_pathway_mat_sigpw_clinical)[2:ncol(tcga_exon_pathway_mat_sigpw_clinical)])) {
  
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp <- na.omit(tcga_exon_pathway_mat_sigpw_clinical_merge[,c("OS.months.","Survival.status",i)])
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  fit_os_tcga_exon_pathway_list[[i]] <- survfit(as.formula(paste("Surv(OS.months., Survival.status) ~", i)), data = tcga_exon_pathway_mat_sigpw_clinical_merge_temp)
  year_survival_os_exon_tcga_pathway_df <- rbind(year_survival_os_exon_tcga_pathway_df,
                                                 data.frame(flag=gsub("=","_pathway_exon_",as.vector(summary(fit_os_tcga_exon_pathway_list[[i]],times = 12)$strata)),
                                                            `OS 1 Year Survival`=  summary(fit_os_tcga_exon_pathway_list[[i]],times = 12,extend = TRUE)$surv,
                                                            `OS 2 Year Survival`=  summary(fit_os_tcga_exon_pathway_list[[i]],times = 24,extend = TRUE)$surv,
                                                            check.names = FALSE)
  )
  
}


fit_rfs_tcga_exon_pathway_list <- list()
year_survival_rfs_exon_tcga_pathway_df <- data.frame()

for (i in make.names(colnames(tcga_exon_pathway_mat_sigpw_clinical)[2:ncol(tcga_exon_pathway_mat_sigpw_clinical)])) {
  
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp <- na.omit(tcga_exon_pathway_mat_sigpw_clinical_merge[,c("RFS.months.", "Recurrence.status",i)])
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%1,i] <- "MUT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]%in%0,i] <- "WT"
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- factor(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],levels = unique(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i]))
  tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i] <- relevel(tcga_exon_pathway_mat_sigpw_clinical_merge_temp[,i],ref="WT")
  fit_rfs_tcga_exon_pathway_list[[i]] <- survfit(as.formula(paste("Surv(RFS.months., Recurrence.status) ~", i)), data = tcga_exon_pathway_mat_sigpw_clinical_merge_temp)
  year_survival_rfs_exon_tcga_pathway_df <- rbind(year_survival_rfs_exon_tcga_pathway_df,
                                                  data.frame(flag=gsub("=","_pathway_exon_",as.vector(summary(fit_rfs_tcga_exon_pathway_list[[i]],times = 12)$strata)),
                                                             `RFS 1 Year Survival`=  summary(fit_rfs_tcga_exon_pathway_list[[i]],times = 12,extend = TRUE)$surv,
                                                             `RFS 2 Year Survival`=  summary(fit_rfs_tcga_exon_pathway_list[[i]],times = 24,extend = TRUE)$surv,
                                                             check.names = FALSE)
  )
  
}


clinical_TCGA_table1 <- TCGA_Clinical.tidy.lusc.maf[,c("OS.months.","RFS.months.",TCGA_exon_gene_sigpw)]
tcga_exon_pathway_mat_sigpw_clinical2 <- tcga_exon_pathway_mat_sigpw_clinical
colnames(tcga_exon_pathway_mat_sigpw_clinical2) <- paste0(colnames(tcga_exon_pathway_mat_sigpw_clinical2),"_pathway")
rownames(tcga_exon_pathway_mat_sigpw_clinical2) <- tcga_exon_pathway_mat_sigpw_clinical2$patients
clinical_TCGA_table1[,ncol(clinical_TCGA_table1)+1:(ncol(tcga_exon_pathway_mat_sigpw_clinical2)-1)] <- tcga_exon_pathway_mat_sigpw_clinical2[rownames(clinical_TCGA_table1),c(2:ncol(tcga_exon_pathway_mat_sigpw_clinical2))]
clinical_TCGA_table1[,c(3:ncol(clinical_TCGA_table1))][clinical_TCGA_table1[,c(3:ncol(clinical_TCGA_table1))] == 1] <- "MUT"
clinical_TCGA_table1[,c(3:ncol(clinical_TCGA_table1))][clinical_TCGA_table1[,c(3:ncol(clinical_TCGA_table1))] == 0] <- "WT"

colnames(clinical_TCGA_table1) <- make.names(colnames(clinical_TCGA_table1))

TCGA_tableOne_all  <- CreateTableOne(data = clinical_TCGA_table1)

TCGA_tableOne_all_out <- print(TCGA_tableOne_all,nonnormal = c("RFS.months.","OS.months."),showAllLevels = TRUE)
TCGA_tableOne_all_out <- TCGA_tableOne_all_out[-1,]

for (i in seq(1, nrow(TCGA_tableOne_all_out), by = 2)) {
  if (i + 1 <= nrow(TCGA_tableOne_all_out)) {
    rownames(TCGA_tableOne_all_out)[i + 1] <- rownames(TCGA_tableOne_all_out)[i]
  }
}

TCGA_tableOne_all_out <- data.frame(name=gsub("\ \\(%\\)","",rownames(TCGA_tableOne_all_out)),TCGA_tableOne_all_out)
rownames(TCGA_tableOne_all_out) <- NULL
TCGA_tableOne_all_out <- data.frame(flag=paste0(TCGA_tableOne_all_out$name,"_exon_",TCGA_tableOne_all_out$level),TCGA_tableOne_all_out)






TCGA_tableOne_list <- list()
TCGA_tableOne_list2 <- list()
TCGA_tableOne_list2_out <- data.frame()

for (i in colnames(clinical_TCGA_table1)[3:ncol(clinical_TCGA_table1)]) {
  
  clinical_TCGA_table1_2 <- clinical_TCGA_table1[,c("RFS.months.","OS.months.",i)]
  
  TCGA_tableOne_list[[i]]  <- CreateTableOne(data = clinical_TCGA_table1_2,
                                        strata = colnames(clinical_TCGA_table1_2)[3],
                                        factorVars =  colnames(clinical_TCGA_table1_2)[3])
  
  TCGA_tableOne_list2[[i]] <- t(print(TCGA_tableOne_list[[i]],nonnormal = c("RFS.months.","OS.months."),showAllLevels = TRUE)[c(2:3),c(2:3)])
  TCGA_tableOne_list2[[i]] <- data.frame(name1=rep(i,ncol(TCGA_tableOne_list2[[i]])),name2=rownames(TCGA_tableOne_list2[[i]]),TCGA_tableOne_list2[[i]],check.names = FALSE)
  rownames(TCGA_tableOne_list2[[i]]) <- NULL
  TCGA_tableOne_list2_out <- rbind(TCGA_tableOne_list2_out,TCGA_tableOne_list2[[i]])
  
}

TCGA_tableOne_list2_out <- data.frame(flag=paste0(TCGA_tableOne_list2_out$name1,"_exon_",TCGA_tableOne_list2_out$name2),TCGA_tableOne_list2_out,check.names = FALSE)


TCGA_tableOne_out2 <- merge(TCGA_tableOne_all_out[,c("flag","Overall")],TCGA_tableOne_list2_out,by="flag")



TCGA_tableOne_out <- rbind(tcga_exon_gene_multicox_rfs_os_merge_survial2,exon_tcga_pathway_multicox_os_rfs_merge_survial2)
year_survival_TCGA_df <- rbind(merge(year_survival_os_exon_tcga_df,year_survival_rfs_exon_tcga_df,by="flag"),
                            merge(year_survival_os_exon_tcga_pathway_df,year_survival_rfs_exon_tcga_pathway_df,by="flag")
)

TCGA_tableOne_out3 <- merge(TCGA_tableOne_out2,TCGA_tableOne_out,by="flag")
TCGA_tableOne_out4 <-  merge(TCGA_tableOne_out3,year_survival_TCGA_df,by="flag")

TCGA_tableOne_out4$pathway <- "FALSE"
TCGA_tableOne_out4[TCGA_tableOne_out4$flag %in% year_survival_rfs_exon_tcga_pathway_df$flag,"pathway"]<- "TRUE"
TCGA_tableOne_out4$name1 <- paste0(TCGA_tableOne_out4$name1,"_exon")

TCGA_tableOne_out5 <- TCGA_tableOne_out4[,c("name1","name2","Overall","RFS.months. (median [IQR])",
                                            "RFS 1 Year Survival","RFS 2 Year Survival","RFS_CI",
                                            "RFS_P_value","Overall","OS.months. (median [IQR])",
                                            "OS 1 Year Survival","OS 2 Year Survival",
                                            "OS_CI","OS_P_value","pathway")]

TCGA_tableOne_out5$`RFS 1 Year Survival` <- signif(TCGA_tableOne_out5$`RFS 1 Year Survival`,digits = 3)
TCGA_tableOne_out5$`RFS 2 Year Survival` <- signif(TCGA_tableOne_out5$`RFS 2 Year Survival`,digits = 3)

TCGA_tableOne_out5$`OS 1 Year Survival` <- signif(TCGA_tableOne_out5$`OS 1 Year Survival`,digits = 3)
TCGA_tableOne_out5$`OS 2 Year Survival` <- signif(TCGA_tableOne_out5$`OS 2 Year Survival`,digits = 3)

empty_row <- data.frame()


TCGA_tableOne_out6 <- do.call(rbind, lapply(1:nrow(TCGA_tableOne_out5), function(i) {
  if (i %% 2 == 1) {
    return(rbind(c(TCGA_tableOne_out5[i,1],rep(NA,14)), TCGA_tableOne_out5[i, , drop = FALSE]))
  } else {
    return(TCGA_tableOne_out5[i, , drop = FALSE])
  }
}))




clinical_109_table1[clinical_109_table1$Stage%in%c("I"),"Stage"] <- 1
clinical_109_table1[clinical_109_table1$Stage%in%c("II"),"Stage"] <- 2
clinical_109_table1[clinical_109_table1$Stage%in%c("III"),"Stage"] <- 3
clinical_109_table1[clinical_109_table1$Stage%in%c("IV"),"Stage"] <- 4
clinical_109_table1$Stage <- as.numeric(clinical_109_table1$Stage)

clinical_109_table1[clinical_109_table1$Gender=="Male","Gender"] <- 1
clinical_109_table1[clinical_109_table1$Gender=="Female","Gender"] <- 0
clinical_109_table1$Gender <- as.numeric(clinical_109_table1$Gender)

colnames(clinical_109_table1) <- make.names(colnames(clinical_109_table1))


cox_results_109_exon_OS <- list()
cli_factor <- c("Age", "Stage","Smoking_history")
for (i in c("FAT4_exon","DCHS1_exon")) {
  formula_string <- paste0("Surv(`OS.months.`, `Survival_status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  cox_results_109_exon_OS[[i]] <- coxph(as.formula(formula_string), data = clinical_109_table1, na.action = na.omit)
}


cox_results_109_exon_RFS <- list()
cli_factor <- c("Age", "Stage","Smoking_history")
for (i in c("FAT4_exon","DCHS1_exon")) {
  formula_string <- paste0("Surv(`RFS.months.`, `Recurrence_status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  cox_results_109_exon_RFS[[i]] <- coxph(as.formula(formula_string), data = clinical_109_table1, na.action = na.omit)
}



cox_results_109_promoter_OS <- list()

for (i in c("TLE1_promoter","FAT3_promoter","SMAD2_promoter")) {
  formula_string <- paste0("Surv(`OS.months.`, `Survival_status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  cox_results_109_promoter_OS[[i]] <- coxph(as.formula(formula_string), data = clinical_109_table1, na.action = na.omit)
}

cox_results_109_promoter_RFS <- list()

for (i in c("WNT10B_promoter")) {
  formula_string <- paste0("Surv(`RFS.months.`, `Recurrence_status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  cox_results_109_promoter_RFS[[i]] <- coxph(as.formula(formula_string), data = clinical_109_table1, na.action = na.omit)
}


cox_results_109_promoter_pathway_OS <- list()

for (i in c("TGF.Beta","Hippo")) {
  formula_string <- paste0("Surv(`OS.months.`, `Survival_status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  cox_results_109_promoter_pathway_OS[[i]] <- coxph(as.formula(formula_string), data = clinical_109_table1, na.action = na.omit)
}



#
TCGA_Clinical.tidy.lusc.maf2 <- TCGA_Clinical.tidy.lusc.maf
TCGA_Clinical.tidy.lusc.maf2 <- TCGA_Clinical.tidy.lusc.maf2[,c("Age", "Stage","Smoking_history","Survival.status","OS.months.","Recurrence.status",
                                                                "RFS.months.","ERBB4","IRS1","NOTCH2","FBXW7","NTRK2","NOTCH3")]
colnames(TCGA_Clinical.tidy.lusc.maf2) <- c("Age", "Stage","Smoking_history","Survival.status","OS.months.","Recurrence.status",
                                            "RFS.months.","ERBB4_exon","IRS1_exon","NOTCH2_exon","FBXW7_exon","NTRK2_exon","NOTCH3_exon")



TCGA_cox_results_109_exon_OS <- list()
cli_factor <- c("Age", "Stage","Smoking_history")
for (i in c("ERBB4_exon","IRS1_exon","NOTCH2_exon")) {
  formula_string <- paste0("Surv(`OS.months.`, `Survival.status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  TCGA_cox_results_109_exon_OS[[i]] <- coxph(as.formula(formula_string), data = TCGA_Clinical.tidy.lusc.maf2, na.action = na.omit)
}


TCGA_cox_results_109_exon_RFS <- list()
cli_factor <- c("Age", "Stage","Smoking_history")
for (i in c("FBXW7_exon","NTRK2_exon","NOTCH3_exon")) {
  formula_string <- paste0("Surv(`RFS.months.`, `Recurrence.status`) ~ ", "`",paste(cli_factor, collapse = "` + `"),"` + `", i,"`")
  TCGA_cox_results_109_exon_RFS[[i]] <- coxph(as.formula(formula_string), data = TCGA_Clinical.tidy.lusc.maf2, na.action = na.omit)
}

rownames(tcga_exon_gene_multicox_os_survial) <- tcga_exon_gene_multicox_os_survial$Gene
rownames(tcga_exon_gene_multicox_rfs_survial) <- tcga_exon_gene_multicox_rfs_survial$Gene

tcga_exon_gene_multicox_os_survial_select <- tcga_exon_gene_multicox_os_survial[c("ERBB4","IRS1","NOTCH2"),]
tcga_exon_gene_multicox_rfs_survial_select <- tcga_exon_gene_multicox_rfs_survial[c("FBXW7","NTRK2"),]



tcga_exon_gene_multicox_os_survial_select[,"95%CI"] <- paste0(signif(tcga_exon_gene_multicox_os_survial_select$OS_HR,digits = 3),"(",
                                                          signif(tcga_exon_gene_multicox_os_survial_select$OS_lower_95,digits = 3),"-",
                                                          signif(tcga_exon_gene_multicox_os_survial_select$OS_upper_95,digits = 3)
                                                          ,")")

tcga_exon_gene_multicox_os_survial_select$OS_P_value <- signif(tcga_exon_gene_multicox_os_survial_select$OS_P_value,digits = 3)

tcga_exon_gene_multicox_os_survial_select$Gene <- paste0(tcga_exon_gene_multicox_os_survial_select$Gene,"_exon")

tcga_exon_gene_multicox_rfs_survial_select[,"95%CI"] <- paste0(signif(tcga_exon_gene_multicox_rfs_survial_select$RFS_HR,digits = 3),"(",
                                                              signif(tcga_exon_gene_multicox_rfs_survial_select$RFS_lower_95,digits = 3),"-",
                                                              signif(tcga_exon_gene_multicox_rfs_survial_select$RFS_upper_95,digits = 3)
                                                              ,")")

tcga_exon_gene_multicox_rfs_survial_select$RFS_P_value <- signif(tcga_exon_gene_multicox_rfs_survial_select$RFS_P_value,digits = 3)

tcga_exon_gene_multicox_rfs_survial_select$Gene <- paste0(tcga_exon_gene_multicox_rfs_survial_select$Gene,"_exon")



#forestplot
rownames(multicox_OS_pathway_exon_genes) <- multicox_OS_pathway_exon_genes$Gene
rownames(multicox_RFS_pathway_exon_genes) <- multicox_RFS_pathway_exon_genes$Gene
rownames(multicox_OS_pathway_promoter_genes) <- multicox_OS_pathway_promoter_genes$Gene
rownames(multicox_RFS_pathway_promoter_genes) <- multicox_RFS_pathway_promoter_genes$Gene
rownames(promoter_pathway_multicox_os_survial) <- promoter_pathway_multicox_os_survial$Gene
rownames(promoter_pathway_multicox_rfs_survial) <- promoter_pathway_multicox_rfs_survial$Gene

multicox_OS_pathway_exon_genes_select <- multicox_OS_pathway_exon_genes[c("FAT4"),]
multicox_RFS_pathway_exon_genes_select <- multicox_RFS_pathway_exon_genes[c("FAT4"),]
multicox_OS_pathway_promoter_genes_select <- multicox_OS_pathway_promoter_genes[c("FAT3","TLE1","SMAD2"),]
multicox_RFS_pathway_promoter_genes_select <- multicox_RFS_pathway_promoter_genes[c("WNT10B"),]
promoter_pathway_multicox_os_survial_select <- promoter_pathway_multicox_os_survial[c("TGF.Beta","Hippo"),]


multicox_OS_pathway_exon_genes_select[,"95%CI"] <- paste0(signif(multicox_OS_pathway_exon_genes_select$OS_HR,digits = 3),"(",
                                                          signif(multicox_OS_pathway_exon_genes_select$OS_lower_95,digits = 3),"-",
                                                          signif(multicox_OS_pathway_exon_genes_select$OS_upper_95,digits = 3)
                                                          ,")")

multicox_OS_pathway_exon_genes_select$OS_P_value <- signif(multicox_OS_pathway_exon_genes_select$OS_P_value,digits = 3)

multicox_OS_pathway_exon_genes_select$Gene <- paste0(multicox_OS_pathway_exon_genes_select$Gene,"_exon")

multicox_RFS_pathway_exon_genes_select[,"95%CI"] <- paste0(signif(multicox_RFS_pathway_exon_genes_select$RFS_HR,digits = 3),"(",
                                                           signif(multicox_RFS_pathway_exon_genes_select$RFS_lower_95,digits = 3),"-",
                                                           signif(multicox_RFS_pathway_exon_genes_select$RFS_upper_95,digits = 3)
                                                           ,")")

multicox_RFS_pathway_exon_genes_select$RFS_P_value <- signif(multicox_RFS_pathway_exon_genes_select$RFS_P_value,digits = 3)

multicox_RFS_pathway_exon_genes_select$Gene <- paste0(multicox_RFS_pathway_exon_genes_select$Gene,"_exon")


multicox_OS_pathway_promoter_genes_select[,"95%CI"] <- paste0(signif(multicox_OS_pathway_promoter_genes_select$OS_HR,digits = 3),"(",
                                                              signif(multicox_OS_pathway_promoter_genes_select$OS_lower_95,digits = 3),"-",
                                                              signif(multicox_OS_pathway_promoter_genes_select$OS_upper_95,digits = 3)
                                                              ,")")

multicox_OS_pathway_promoter_genes_select$OS_P_value <- signif(multicox_OS_pathway_promoter_genes_select$OS_P_value,digits = 3)

multicox_OS_pathway_promoter_genes_select$Gene <- paste0(multicox_OS_pathway_promoter_genes_select$Gene,"_promoter")


multicox_RFS_pathway_promoter_genes_select[,"95%CI"] <- paste0(signif(multicox_RFS_pathway_promoter_genes_select$RFS_HR,digits = 3),"(",
                                                               signif(multicox_RFS_pathway_promoter_genes_select$RFS_lower_95,digits = 3),"-",
                                                               signif(multicox_RFS_pathway_promoter_genes_select$RFS_upper_95,digits = 3)
                                                               ,")")

multicox_RFS_pathway_promoter_genes_select$RFS_P_value <- signif(multicox_RFS_pathway_promoter_genes_select$RFS_P_value,digits = 3)

multicox_RFS_pathway_promoter_genes_select$Gene <- paste0(multicox_RFS_pathway_promoter_genes_select$Gene,"_promoter")


promoter_pathway_multicox_os_survial_select[,"95%CI"] <- paste0(signif(promoter_pathway_multicox_os_survial_select$OS_HR,digits = 3),"(",
                                                                signif(promoter_pathway_multicox_os_survial_select$OS_lower_95,digits = 3),"-",
                                                                signif(promoter_pathway_multicox_os_survial_select$OS_upper_95,digits = 3)
                                                                ,")")
promoter_pathway_multicox_os_survial_select$OS_P_value <- signif(promoter_pathway_multicox_os_survial_select$OS_P_value,digits = 3)

promoter_pathway_multicox_os_survial_select$Gene <- paste0(promoter_pathway_multicox_os_survial_select$Gene,"_promoter")


multicox_OS_merge <- rbind(multicox_OS_pathway_exon_genes_select,
                           multicox_OS_pathway_promoter_genes_select,
                           promoter_pathway_multicox_os_survial_select)

multicox_OS_merge <- rbind(c("Gene",NA,NA,NA,"P_value","HR(95%CI)"),multicox_OS_merge)

multicox_OS_merge[,2:4] <- lapply(multicox_OS_merge[,2:4],as.numeric)


multicox_RFS_merge <- rbind(multicox_RFS_pathway_exon_genes_select,
                            multicox_RFS_pathway_promoter_genes_select)

multicox_RFS_merge <- rbind(c("Gene",NA,NA,NA,"P_value","HR(95%CI)"),multicox_RFS_merge)

multicox_RFS_merge[,2:4] <- lapply(multicox_RFS_merge[,2:4],as.numeric)


TCGA_multicox_OS_merge <- tcga_exon_gene_multicox_os_survial_select

TCGA_multicox_OS_merge <- rbind(c("Gene",NA,NA,NA,"P_value","HR(95%CI)"),TCGA_multicox_OS_merge)

TCGA_multicox_OS_merge[,2:4] <- lapply(TCGA_multicox_OS_merge[,2:4],as.numeric)


TCGA_multicox_RFS_merge <- tcga_exon_gene_multicox_rfs_survial_select

TCGA_multicox_RFS_merge <- rbind(c("Gene",NA,NA,NA,"P_value","HR(95%CI)"),TCGA_multicox_RFS_merge)

TCGA_multicox_RFS_merge[,2:4] <- lapply(TCGA_multicox_RFS_merge[,2:4],as.numeric)






is_summary <- c(TRUE, rep(FALSE, nrow(multicox_OS_merge) - 1))
library(forestplot)

pdf("forestplot_OS.pdf",height = 5)

forestplot::forestplot(multicox_OS_merge,
                       labeltext = c(Gene,`95%CI`,OS_P_value),
                       mean = OS_HR,
                       lower = OS_lower_95,
                       upper=  OS_upper_95,
                       graph.pos = 3,
                       boxsize = 0.15,
                       xlab = "Hazard Ratio (HR)",
                       clip = c(0, 10),
                       zero = 1,
                       new_page=FALSE,
                       txt_gp = fpTxtGp(
                         xlab = gpar(fontsize = 20), 
                         ticks = gpar(fontsize = 20),
                         label = gpar(fontsize = 12),
                         summary = gpar(fontsize = 20, fontface = "bold")
                       ),
                       is_summary = is_summary,
                       col = fpColors(box = "darkblue", line = "darkblue", summary = "darkred"),
                       lineheight = "auto", 
                       title = "Overall survival"
)

dev.off()

pdf("forestplot_RFS.pdf",height = 3)

forestplot::forestplot(multicox_RFS_merge,
                       labeltext = c(Gene,`95%CI`,RFS_P_value),
                       mean = RFS_HR,
                       lower = RFS_lower_95,
                       upper=  RFS_upper_95,
                       graph.pos = 3,
                       boxsize = 0.15,
                       xlab = "Hazard Ratio (HR)",
                       clip = c(0, 10),
                       zero = 1,
                       new_page=FALSE,
                       txt_gp = fpTxtGp(
                         xlab = gpar(fontsize = 20), 
                         ticks = gpar(fontsize = 20),
                         label = gpar(fontsize = 12),
                         summary = gpar(fontsize = 20, fontface = "bold")
                       ),
                       is_summary = is_summary,
                       col = fpColors(box = "darkblue", line = "darkblue", summary = "darkred"), 
                       lineheight = "auto", 
                       title = "Recurrence free survival"
)

dev.off()



pdf("TCGA_forestplot_OS.pdf",height = 3)

forestplot::forestplot(TCGA_multicox_OS_merge,
                       labeltext = c(Gene,`95%CI`,OS_P_value),
                       mean = OS_HR,
                       lower = OS_lower_95,
                       upper=  OS_upper_95,
                       graph.pos = 3,
                       boxsize = 0.15,
                       xlab = "Hazard Ratio (HR)",
                       clip = c(0, 10),
                       zero = 1,
                       new_page=FALSE,
                       txt_gp = fpTxtGp(
                         xlab = gpar(fontsize = 20), 
                         ticks = gpar(fontsize = 20),
                         label = gpar(fontsize = 12),
                         summary = gpar(fontsize = 20, fontface = "bold")
                       ),
                       is_summary = is_summary,
                       col = fpColors(box = "darkblue", line = "darkblue", summary = "darkred"),
                       lineheight = "auto", 
                       title = "Overall survival"
)

dev.off()

pdf("TCGA_forestplot_RFS.pdf",height = 3)

forestplot::forestplot(TCGA_multicox_RFS_merge,
                       labeltext = c(Gene,`95%CI`,RFS_P_value),
                       mean = RFS_HR,
                       lower = RFS_lower_95,
                       upper=  RFS_upper_95,
                       graph.pos = 3,
                       boxsize = 0.15,
                       xlab = "Hazard Ratio (HR)",
                       clip = c(0, 10),
                       zero = 1,
                       new_page=FALSE,
                       txt_gp = fpTxtGp(
                         xlab = gpar(fontsize = 20), 
                         ticks = gpar(fontsize = 20),
                         label = gpar(fontsize = 12),
                         summary = gpar(fontsize = 20, fontface = "bold")
                       ),
                       is_summary = is_summary,
                       col = fpColors(box = "darkblue", line = "darkblue", summary = "darkred"), 
                       lineheight = "auto",
                       title = "Recurrence free survival"
                       )

dev.off()

cox_results_tcga_exon_OS <- list()
cli_factor_tcga <- c("Age", "Stage","Smoking_history")
for (i in c("ERBB4","IRS1","NOTCH2")) {
  formula_string <- paste0("Surv(`OS.months.`, `Survival.status`) ~ ", "`",paste(cli_factor_tcga, collapse = "` + `"),"` + `", i,"`")
  cox_results_tcga_exon_OS[[i]] <- coxph(as.formula(formula_string), data = TCGA_Clinical.tidy.lusc.maf, na.action = na.omit)
}

cox_results_tcga_exon_RFS <- list()
cli_factor_tcga <- c("Age", "Stage","Smoking_history")
for (i in c("FBXW7","NOTCH3","NTRK2")) {
  formula_string <- paste0("Surv(`RFS.months.`, `Recurrence.status`) ~ ", "`",paste(cli_factor_tcga, collapse = "` + `"),"` + `", i,"`")
  cox_results_tcga_exon_RFS[[i]] <- coxph(as.formula(formula_string), data = TCGA_Clinical.tidy.lusc.maf, na.action = na.omit)
}




T <- function(x){
  A = matrix(0,nrow(x),ncol(x))
  for(i in 1:nrow(x)){
    for(j in 1:ncol(x)) A[i,j] = sum(x[i,])*sum(x[,j])/sum(x)
  }
  A
}

pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  TT <- T(table(y, g))
  nn <- sum(table(y, g))
  
  if(sum(TT) >=5 & nn>=40 ){
    p <- chisq.test(table(y, g),correct = FALSE)$p.value
  }
  if(sum(TT) <5 & sum(TT)>=1 & nn>=40 ){
    p <- chisq.test(table(y, g),correct = TRUE)$p.value
  }
  if(sum(TT)<1 | nn<40 ){
    p <- fisher.test(table(y, g))$p.value
  }
  # sub("<", "&lt;", format.pval(p, digits=4, eps=0.0001))
  p
}


exon_pathway_multicox_merge_survial <- merge(exon_pathway_multicox_os_survial,exon_pathway_multicox_rfs_survial,by="Gene")
exon_pathway_multicox_merge_survial$cohort <- "109"
exon_pathway_multicox_merge_survial$Mutation_Number <- as.vector(unlist(exon_pathway_mat_sigpw_clinical_merge_length[exon_pathway_multicox_merge_survial$Gene]))








