
##################################################################
#############################s############### COAD ############################
##################################################################


# Raw.Data.Path<-"E:/Dropbox/Work_temp/IO_nanoStrings_CRC/TCGA"
# setwd(Raw.Data.Path)

##### View survival_COAD_survival.txt  #####
COAD.Clinical<-read.table(file="data/survival_COAD_survival.txt",sep="\t",head=TRUE,quote = "\"")
dim(COAD.Clinical)
colnames(COAD.Clinical)
head(COAD.Clinical)
colSums(is.na(COAD.Clinical[colSums(is.na(COAD.Clinical)) > 0]))



##### View TCGA.COAD.sampleMap_COAD_clinicalMatrix  #####
COAD.Clinical2<-read.table(file="data/TCGA.COAD.sampleMap_COAD_clinicalMatrix",sep="\t",head=TRUE,quote = "\"")
dim(COAD.Clinical2)
colnames(COAD.Clinical2)
head(COAD.Clinical2[,1:10])
colSums(is.na(COAD.Clinical2[colSums(is.na(COAD.Clinical2)) > 0]))


# Column categories

# 1. 基本病人資料
basic_info <- c(
  "sampleID", "patient_id", "bcr_patient_barcode", "bcr_sample_barcode", "bcr_followup_barcode",
  "gender", "age_at_initial_pathologic_diagnosis", "days_to_birth", "days_to_death",
  "days_to_last_followup", "vital_status", "height", "weight", "initial_weight",
  "year_of_initial_pathologic_diagnosis", "system_version", "form_completion_date"
)

# 2. 分子與基因相關資訊
molecular_gene <- c(
  "MSI_updated_Oct62011", "microsatellite_instability", "CIMP",
  "braf_gene_analysis_performed", "braf_gene_analysis_result",
  "kras_gene_analysis_performed", "kras_mutation_found", "kras_mutation_codon",
  "hypermutation", "loss_expression_of_mismatch_repair_proteins_by_ihc",
  "loss_expression_of_mismatch_repair_proteins_by_ihc_result",
  "non_silent_mutation", "silent_mutation", "total_mutation",
  "non_silent_rate_per_Mb", "silent_rate_per_Mb",
  "number_of_abnormal_loci", "number_of_loci_tested"
)

# 3. 臨床診斷與病理
clinical_pathology <- c(
  "primary_site", "tumor_tissue_site", "anatomic_neoplasm_subdivision",
  "histological_type", "pathologic_T", "pathologic_N", "pathologic_M",
  "pathologic_stage", "residual_tumor",
  "residual_disease_post_new_tumor_event_margin_status",
  "circumferential_resection_margin", "perineural_invasion_present",
  "lymphatic_invasion", "venous_invasion",
  "lymph_node_examined_count", "number_of_lymphnodes_positive_by_he",
  "number_of_lymphnodes_positive_by_ihc", "non_nodal_tumor_deposits", "pathology_report_file_name", "primary_lymph_node_presentation_assessment"
)

# 4. 治療與追蹤
treatment_followup <- c(
  "radiation_therapy", "additional_radiation_therapy", "additional_pharmaceutical_therapy",
  "postoperative_rx_tx", "preoperative_pretreatment_cea_level",
  "primary_therapy_outcome_success", "followup_treatment_success",
  "person_neoplasm_cancer_status",
  "new_tumor_event_after_initial_treatment", "new_neoplasm_event_type",
  "new_tumor_event_additional_surgery_procedure",
  "days_to_new_tumor_event_after_initial_treatment",
  "days_to_new_tumor_event_additional_surgery_procedure",
  "history_of_neoadjuvant_treatment", "lost_follow_up",
  "followup_case_report_form_submission_reason", "site_of_additional_surgery_new_tumor_event_mets"
)

# 5. 其他病史與檢查
history_other <- c(
  "history_of_colon_polyps", "colon_polyps_present", "synchronous_colon_cancer_present",
  "number_of_first_degree_relatives_with_cancer_diagnosis",
  "icd_10", "icd_o_3_histology", "icd_o_3_site",
  "informed_consent_verified",
  "AWG_MLH1_silencing", "AWG_cancer_type_Oct62011", "CDE_ID_3226963", "disease_code"
)

# 6. TCGA 整合與多組學欄位
tcga_multiomics <- c(
  grep("^X_", colnames(COAD.Clinical2), value = TRUE),  # 所有 X_ 開頭的欄位
  "project_code", "X_cohort", "X_primary_disease", "X_primary_site"
)

# 7. 樣本/檢體資訊
sample_info <- c(
  "days_to_collection", "days_to_initial_pathologic_diagnosis", "intermediate_dimension", "longest_dimension", "shortest_dimension",
  "is_ffpe", "oct_embedded", "sample_type", "sample_type_id",
  "tissue_prospective_collection_indicator", "tissue_retrospective_collection_indicator",
  "tissue_source_site", "vial_number"
)


# === Split to 6 data.frame ===
COAD_basic      <- COAD.Clinical2[, basic_info]
COAD_molecular  <- COAD.Clinical2[, molecular_gene]
COAD_pathology  <- COAD.Clinical2[, clinical_pathology]
COAD_treatment  <- COAD.Clinical2[, treatment_followup]
COAD_history    <- COAD.Clinical2[, history_other]
COAD_multiomics <- COAD.Clinical2[, tcga_multiomics]




#### TCGA.COAD.sampleMap_HiSeqV2 AND TCGA.COAD.sampleMap_GAV2 ####
TCGA.temp2<-read.table(file="data/TCGA.COAD.sampleMap_HiSeqV2.data",sep="\t",head=TRUE,quote = "\"")
head(TCGA.temp2[,1:10])
TCGA.temp<-read.table(file="data/TCGA.COAD.sampleMap_GAV2.data",sep="\t",head=TRUE,quote = "\"")
head(TCGA.temp[,1:10])
dim(TCGA.temp)
dim(TCGA.temp2)



# TCGA.temp2[,-1]<-head(TCGA.temp2[,1:10])
# match(colnames(TCGA.temp),"TCGA.AG.3611.01")


# Convert to FPKM
# matrix.temp<-2^TCGA.temp[,-1]-1

# Normalization FPKM
# TCGA.temp2<-log2(matrix.temp*10^6/t(matrix(colSums(matrix.temp),512,60483))+1)
# head(TCGA.temp2[,1:10])

# Match tumor cells by patient
COAD.T.index <- match(unique(paste0(COAD.Clinical$X_PATIENT, "-01")), COAD.Clinical$sample, nomatch = 0)
sum(COAD.T.index != 0)  # Count how many tumor samples were successfully matched

# Match normal cells by patient
COAD.N.index <- match(unique(paste0(COAD.Clinical$X_PATIENT, "-11")), COAD.Clinical$sample, nomatch = 0)
sum(COAD.N.index != 0)  # Count how many normal samples were successfully matched

# Find patients who have both tumor and normal cells
COAD.matchT.index <- match(unique(paste0(COAD.Clinical$X_PATIENT[COAD.N.index], "-01")), COAD.Clinical$sample, nomatch = 0)
cbind(COAD.Clinical$sample[COAD.matchT.index], COAD.Clinical$sample[COAD.N.index])  # Combine matched tumor-normal sample pairs

# Extract tumor patient clinical data
COAD.Clinical.T <- COAD.Clinical[COAD.T.index,]

# Identify the patient's tumor cells using GAV2
common.names <- intersect(substr(colnames(TCGA.temp), 1, 15), make.names(COAD.Clinical.T$sample))
common.names
length(common.names)  # Count common names

# Identify the patient's tumor cells using HiSeqV2
common.names2 <- intersect(substr(colnames(TCGA.temp2), 1, 15), make.names(COAD.Clinical.T$sample))
common.names2
length(common.names2)  # Count common names

# Identify the patient's normal cells using HiSeqV2
common.namesN <- intersect(substr(colnames(TCGA.temp2), 1, 15), make.names(COAD.Clinical$sample[COAD.N.index]))
common.namesN
length(common.namesN)  # Count common names



#### Analysis ####
