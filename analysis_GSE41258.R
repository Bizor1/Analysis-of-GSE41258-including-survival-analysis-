library(GEOquery)
library(DESeq2)
library(IRanges)
library(limma)
library(ggplot2)
library(Biobase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(KEGG.db)
library(KEGGREST)
library(affy)
library(hgu133plus2.db)  
library(pheatmap)
library(STRINGdb)
library(TCGAbiolinks)
library(xml2)
library(dplyr)
library(XML)
library(AnnotationDbi)
library(biomaRt)
library(jsonlite)
library(survival)
library(survminer)
library(enrichplot)
library(DESeq2)

install.packages("survival") 
install.packages("survminer") 
install.packages("xml2")
BiocManager::install("STRINGdb")
BiocManager::install("TCGAbiolinks")
BiocManager::install("hgu133plus2.db",force = TRUE)


#Gene Expresion analysis


if(!requireNamespace("clusterProfiler", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("clusterProfiler")
}
library(clusterProfiler)

# Load the data
GSE41258 <- getGEO("GSE41258", GSEMatrix = TRUE)[[1]]
exprs_GSE41258 <- exprs(GSE41258)
pdata_GSE41258 <- pData(GSE41258)

# Filter out NA values
pdata_GSE41258 <- pdata_GSE41258[!is.na(pdata_GSE41258$`tissue:ch1`),]
pdata_GSE41258GSE41258 <- getGEO("GSE41258", GSEMatrix = TRUE)[[1]]
#Error in open.connection(x, "rb") : HTTP error 404
# Subset for Normal Colon and Primary Tumor
selected_samples <- pdata_GSE41258[pdata_GSE41258$`tissue:ch1` %in% c("Normal Colon", "Primary Tumor"),]
selected_samples
group <- factor(selected_samples$`tissue:ch1`, levels = c("Normal Colon", "Primary Tumor"))
group
unique(pdata_GSE41258$'tissue:ch1')
# Create design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Ensure column names match between expression data and selected samples
exprs_GSE41258_filtered <- exprs_GSE41258[, colnames(exprs_GSE41258) %in% rownames(selected_samples)]

# Fit model
fit <- lmFit(exprs_GSE41258_filtered, design)
fit<-eBayes(fit)

topGenes=topTable(fit,coef = 2)


dim(exprs_GSE41258_filtered)

dim(design)
nrow(selected_samples)
print(colnames(design))



# Define all relevant tissue pairs
tissue_pairs <- list(
  c("Normal Colon", "Primary Tumor"),
  c("Normal Liver", "Liver Metastasis"),
  c("Normal Lung", "Lung Metastasis"),
  c("Polyp", "Polyp, high grade")
)

# Initialize a list to store results
results <- list()

# Loop over each pair
for (pair in tissue_pairs) {
  tissue1 <- pair[1]
  tissue2 <- pair[2]
  
  # Subset the pdata for the two tissues
  selected_samples <- pdata_GSE41258[pdata_GSE41258$`tissue:ch1` %in% c(tissue1, tissue2),]
  
  # Create group factor for the pair
  group <- factor(selected_samples$`tissue:ch1`, levels = c(tissue1, tissue2))
  
  # Create design matrix with valid names
  design <- model.matrix(~0 + group)
  valid_names <- make.names(levels(group))
  colnames(design) <- valid_names
  
  # Ensure column names match between expression data and selected samples
  selected_sample_ids <- rownames(selected_samples)
  exprs_GSE41258_filtered <- exprs_GSE41258[, colnames(exprs_GSE41258) %in% selected_sample_ids]
  
  # Perform the linear model fit
  if (ncol(exprs_GSE41258_filtered) > 0 && nrow(design) > 0) {
    fit <- lmFit(exprs_GSE41258_filtered, design)
    
    # Define dynamic contrast name
    contrast_name1 <- valid_names[1]
    contrast_name2 <- valid_names[2]
    
    # Define contrast matrix using exact column names
    contrast_matrix <- makeContrasts(
      contrast = paste0(contrast_name1, "-", contrast_name2),
      levels = design
    )
    
    # Apply contrasts and fit
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # Store results with dynamic name
    results[[paste0(contrast_name1, "_vs_", contrast_name2)]] <- fit2
    print(paste("Model fitting completed for", tissue1, "vs", tissue2))
  } else {
    print(paste("Error: No data for", tissue1, "vs", tissue2))
  }
}

# Results will store the fitted models and contrasts for all tissue pairs
# Extract results from the fit object
results_list <- lapply(results, function(fit2) {
  topTable(fit2, number = Inf)  # Extract all results
})
# Extract results from the fit object
results_list <- lapply(results, function(fit2) {
  topTable(fit2, number = Inf)  # Extract all results
})

results_list
# Apply multiple testing correction
results_adjusted <- lapply(results_list, function(res) {
  res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")  # Adjust p-values using Benjamini-Hochberg method
  res
})
results_adjusted
# Example for one pair
res <- results_adjusted[[2]]
res
ggplot(res, aes(x = logFC, y = -log10(P.Value), color = adj.P.Val < 0.05)) +
  geom_point() +
  labs(title = "Volcano Plot", x = "Log Fold Change", y = "-Log10 P-Value")


results[Polyp_vs_Polyp..high.grade]
res

# Sort the results based on adjusted p-value

sorted_results<- res[order(res$adj.P.Val),]


sorted_results
# Extract the top 10 most significant genes
top10_genes <- head(sorted_results, 10)
significant_genes <- sorted_results[sorted_results$adj.P.Val < 0.01, ]
significant_genes 
significant_genes$probe_id<-rownames(significant_genes)
significant_genes

upregulated_genes <- significant_genes[significant_genes$logFC > 0 & significant_genes$adj.P.Val < 1e-20, ]
upregulated_genes
downregulated_genes <- significant_genes[significant_genes$logFC < 0 & significant_genes$adj.P.Val <1e-10 , ]
downregulated_genes
top10_genes
top10_genes$probe_id <- rownames(top10_genes)
#alternative
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
print(ensembl)

probe_ids <- rownames(top10_genes)

probe_ids<-rownames(significant_genes)


probe_ids



annotations <- getBM(attributes = c("affy_hg_u133_plus_2", "hgnc_symbol", "gene_biotype"),
                     filters = "affy_hg_u133_plus_2", 
                     values = probe_ids, 
                     mart = ensembl)
print(annotations)


print(top10_genes)
print(rownames(top10_genes))
print(unique(rownames(significant_genes)))

keytypes(hgu133plus2.db)
available_probes <- keys(hgu133plus2.db, keytype = "PROBEID")
print(available_probes[1:30])  # Show first 10 available probes
gene_annotations <- select(hgu133plus2.db, 
                           keys = rownames(top10_genes), 
                           columns = c("SYMBOL", "GENENAME"), 
                           keytype = "PROBEID")
gene_annotations <- select(hgu133plus2.db, 
                           keys = rownames(significant_genes), 
                           columns = c("SYMBOL", "GENENAME"), 
                           keytype = "PROBEID")


annotated_top10<-merge(top10_genes, annotations, by.x = "probe_id", by.y = "affy_hg_u133_plus_2")
annotated_significant_genes<-merge(significant_genes, annotations, by.x = "probe_id", by.y = "affy_hg_u133_plus_2")

print(annotated_top10)
print(annotated_significant_genes)



gene_sygene_sygene_symbols<-annotated_significant_genes$hgnc_symbol

gene_logfcs=annotated_significant_genes$logFC
gene_logfcs
length(gene_logfcs)

length(gene_symbols)


valid_gene_symbols <- gene_symbols[gene_symbols != ""]
unique(valid_gene_symbols)

entrez_ids<-bitr(gene_symbols,fromType = "SYMBOL",toType="ENTREZID",OrgDb = "org.Hs.eg.db")
ensembl_ids<-bitr(gene_symbols,fromType = "SYMBOL",toType="ENSEMBL",OrgDb = "org.Hs.eg.db")
entrez_genes_list<-entrez_ids$ENTREZID
ensembl_genes_list<-ensembl_ids$ENSEMBL
length(entrez_genes_list)
length(ensembl_genes_list)
ensembl_genes_list


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
probe_ids <- annotated_significant_genes$probe_id  # The probe IDs you want to map

ensembl_ids <- getBM(filters = "affy_hg_u133_plus_2", 
                     attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene_id"), 
                     values = probe_ids, 
                     mart = ensembl)

merged_data <- merge(annotated_significant_genes, ensembl_ids, by.x = "probe_id", by.y = "affy_hg_u133_plus_2", all.x = TRUE)
merged_data

valid_data <- merged_data[!is.na(merged_data$ensembl_gene_id), ]
valid_data
ordered_data <- valid_data[order(-valid_data$logFC), ]
ordered_data
final_data <- ordered_data[, c("probe_id", "ensembl_gene_id", "entrezgene_id", "logFC")]
final_data

vector1 <- final_data$logFC
names(vector1) <- final_data$ensembl_gene_id
vector1

entrez_ids <- bitr(unique(final_data$ensembl_gene_id), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids
write.csv(entrez_ids[1:1000,], file = "final_significant_genes_with_entrez.csv", row.names = FALSE)

gse <- gseGO(vector1,
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)
gse
fit <- gseaplot(gse, geneSetID = 1)
fit
gse

as.data.frame(gse)


kegg_enrichment<-enrichKEGG(gene=entrez_genes_list,organism='hsa',pvalueCutoff=0.05)



kegg_enrichment

barplot(kegg_enrichment,showCategory=10,tite="Top 10 kegg PATHWAYS")


#gene set enrichment analysis
annotated_significant_genes

analysis_gene_list<-sorted_results$logFC
names(analysis_gene_list)

gse <- gseGO(gene_list,
             ont = "BP",
             keyType = "ENSEMBL",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)




#survival analysis
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  data.type = "Clinical Supplement"
                )
GDCdownload(query)
files <- getResults(query)
files

xml_dir <- "C:/Users/bizor/Documents/GDCdata/TCGA-COAD/Clinical/Clinical_Supplement/"


xml_files <- list.files(path = xml_dir, pattern = "\\.xml$", full.names = TRUE, recursive = TRUE)
xml_files


inspect_structure <- function(file) {
  xml_data <- xmlParse(file)
  patient_data <- xmlToList(xml_data)$patient
  # Print the structure of the patient data
  str(patient_data)
}

# Inspect the structure of a few XML files
lapply(xml_files[1:5], inspect_structure)



# Function to parse XML and convert to JSON
parse_xml_to_json <- function(file) {
  xml_data <- read_xml(file)
  patient_data <- xml_find_all(xml_data, ".//patient")
  
  # Extract data and convert to a list
  data_list <- lapply(patient_data, function(node) {
    list(
      vital_status = xml_text(xml_find_first(node, ".//vital_status")),
      days_to_death = as.numeric(xml_text(xml_find_first(node, ".//days_to_death"))),
      days_to_last_followup = as.numeric(xml_text(xml_find_first(node, ".//days_to_last_followup")))
    )
  })
  
  # Convert list to JSON
  json_data <- toJSON(data_list, pretty = TRUE)
  return(json_data)
}

# Process all XML files and write to JSON
#xml_files <- list.files(path = "xml_dir", pattern = "\\.xml$", full.names = TRUE)
xml_files
json_list <- lapply(xml_files, parse_xml_to_json)
json_list[1]


# Combine all JSON data into one
combined_json <- paste(unlist(json_list), collapse = ",\n")
combined_json <- paste0("[\n", combined_json, "\n]")

# Write JSON data to a file
writeLines(combined_json, "all_clinical_data.json")

all_clinical_data <- fromJSON("all_clinical_data.json")
print(all_clinical_data)


xml_file <- xml_files[1]  # Use the first XML file as an example
xml_file







xml_file <- xml_files[1]
xml_data <- read_xml(xml_file)
xml_data

patients <- xml_find_all(xml_data, ".//coad:patient")
print(xml_text(patients[1]))
print(as.character(xml_data))

namespaces <- xml_ns(xml_data)
print(namespaces)


patients <- xml_find_all(xml_data, ".//coad:patient", ns = namespaces)
patients


json_dir <- "C:/Users/bizor/Documents/GDCdata/TCGA-COAD/Clinical/Clinical_Supplement/"
json_files <- list.files(path = json_dir, pattern = "\\.json$", full.names = TRUE, recursive = TRUE)
json_files
inspect_json_structure <- function(file) {
  json_data <- tryCatch(fromJSON(file, flatten = TRUE), 
                        error = function(e) { 
                          message("Error reading file: ", file)
                          return(NULL) 
                        })
  
  if (is.null(json_data)) {
    return(NULL)
  }
  
  # Check for the presence of key fields
  has_vital_status <- "vital_status" %in% names(json_data)
  has_days_to_death <- "days_to_death" %in% names(json_data)
  has_days_to_last_followup <- "days_to_last_followup" %in% names(json_data)
  
  # Return structure type
  if (has_vital_status && has_days_to_death && has_days_to_last_followup) {
    return("Survival Analysis")
  } else {
    return("Other")
  }
}

# Inspect all JSON files
file_classifications <- sapply(json_files, inspect_json_structure)

# Combine file names with their classifications
classified_files <- data.frame(
  file = json_files,
  classification = file_classifications,
  stringsAsFactors = FALSE
)

# Save classifications to a CSV file
write.csv(classified_files, "classified_files.csv", row.names = FALSE)
survival_files <- classified_files[classified_files$classification == "Survival Analysis", "file"]
other_files <- classified_files[classified_files$classification == "Other", "file"]
# Save lists to files for review
writeLines(survival_files, "survival_files.txt")
writeLines(other_files, "other_files.txt")


#Function to check for relevant fields in the JSON data
check_fields <- function(data) {
  if (is.null(data)) return("Unknown")
  
  fields <- c(
    "coad:patient$clin_shared$vital_status",
    "coad:patient$clin_shared$days_to_death",
    "coad:patient$clin_shared$days_to_last_followup"
  )
  
  # Check presence of fields
  if (all(fields %in% names(data))) {
    if (!is.null(data$`coad:patient$clin_shared$vital_status`) &
        !is.null(data$`coad:patient$clin_shared$days_to_death`) &
        !is.null(data$`coad:patient$clin_shared$days_to_last_followup`)) {
      return("Survival Data")
    }
  }
  
  return("Other Data")
}

# Function to read JSON and categorize
categorize_json_file <- function(file) {
  data <- read_json_file(file)
  category <- check_fields(data)
  return(data.frame(file = basename(file), category = category, stringsAsFactors = FALSE))
}

# Read all JSON files and categorize
categorized_data <- do.call(rbind, lapply(json_files, categorize_json_file))

# Save the categorized results
write.csv(categorized_data, "categorized_json_files.csv", row.names = FALSE)

# Preview the categorized results
head(categorized_data)





data<-read.csv("C:/Users/bizor/Documents/GDCdata/TCGA-COAD/Clinical/survival_analysis_data.csv")

head(data)


data$days_to_death<- as.numeric((as.character(data$days_to_death)))
data$time <- ifelse(data$vital_status == "Dead", data$days_to_death, data$days_to_followup)

data$event <- ifelse(data$vital_status == "Dead", 1, 0)

head(data)

surv_object <- Surv(time = data$time, event = data$event)



print(surv_object)


km_fit<-survfit(surv_object~1)

summary(km_fit)

ggsurvplot(km_fit, data = data, conf.int = TRUE, pval = TRUE, 
           title = "Kaplan-Meier Survival Curve", 
           xlab = "Days", ylab = "Survival Probability")





























# Sample gene list with ENTREZ IDs and log fold changes
sample_gene_list <- c(
  "7157" = 2.5,  # TP53 - Tumor Protein p53
  "5290" = -1.8, # AKT1 - v-akt murine thymoma viral oncogene homolog 1
  "1956" = 0.5,  # MYC - Myelocytomatosis viral oncogene homolog
  "367"  = 3.2,  # HSP90AA1 - Heat Shock Protein 90 Alpha Family Class A Member 1
  "174"  = -0.3, # FOS - FBJ murine osteosarcoma viral oncogene homolog
  "2033" = 1.1,  # JUN - Jun proto-oncogene
  "1026" = -2.4   # EGFR - Epidermal Growth Factor Receptor
)

# Ensure the gene list is sorted in decreasing order
sorted_gene_list <- sort(sample_gene_list, decreasing = TRUE)

# Display the sorted gene list
print(sorted_gene_list)



library(clusterProfiler)
library(org.Hs.eg.db)

# Convert ENTREZ IDs to GO terms
go_genes <- bitr(names(sorted_gene_list), fromType = "ENTREZID", toType = "GO", OrgDb = org.Hs.eg.db)

# Perform GO enrichment analysis
gsea_go <- gseGO(geneList = sorted_gene_list,
                 OrgDb = org.Hs.eg.db,
                 ont = "ALL",  # Biological Process
                 pvalueCutoff = 0.05,
                 scoreType = "pos")

# Check results
head(gsea_go)


