# Analysis-of-GSE41258-including-survival-analysis-
DIFFERENTIAL GENE EXPRESSION, PATHWAY ENRICHMENT AND SURVIVAL ANALYSIS IN COLORECTAL CANCER: INSIGHTS FROM GSE41258 
ABSTRACT
This work explores the molecular and biological underpinnings of colorectal cancer based on the most differentially expressed genes in the GSE41258 expression dataset. The analysis focuses on the differential expression of genes between normal colonic and tumor tissues. Several genes including PLG, TFR2 and SERPINA10 were significantly upregulated whiles others like BTC were significantly downregulated. These upregulated genes had roles in tumor progression, invasion and metastasis. Pathway enrichment revealed alterations in the metabolic, immune and inflammatory related response pathways as a major contributor to tumor progression. Also, survival analysis identified TFR2 and SERPINA10 as poor prognostic markers, making them potential therapeutic targets. This study outlines the ever-growing importance of integrating bioinformatics research into therapeutics while acknowledging the limitations in terms of data and the need for in vivo validation.

 
BACKGROUND 
Colorectal cancer (CRC) is the third most common cancer and the second most common cause of cancer related mortality worldwide (Roshandel et al., 2024). Approximately 9.4% of cancer-related deaths were due to colorectal cancer in 2020.  
CRC occurs is caused by hyper proliferation of glandular epithelial cells. Many factors are related to the development of CRC including genetics, polyposis, chronic inflammation, inflammatory bowel disease, increased body mass index, sedentary lifestyle, cigarrete smoking, alcohol abuse and dietary habits(Wang et al., 2022) .The transition from normal epithelium to dysplasia occurs over time and involves genetic changes and abnormalities. There are three branches of these changes in terms of pathways, and they are: chromosomal instability, Microsatellite Instability(MSI) and CpG island methylation(CIMPG) phenotype (Mármol et al., 2017). 
Chromosomal Instability(CIN) accounts for 80-85% of CRC cases. This instability leads to chromosome segregation, telomere dysfunction and an impaired response to DNA damage. Genes like APC, KRAS,PI3k, and TP53, promote tumor genesis and uncontrolled cell proliferation. MSI, which is defined as somatic alterations in microsatellite sequences due to deletion or insertion of one or more repeat unit, mostly occurs due to defective DNA repair, however, tumors associated with this pathway often have a relatively better prognosis. CIMP often co-occurs with genetic mutations such as BRAF and contributes to tumor development (Simons et al., 2013). 
With the literature being clear on the importance of molecular pathways and gene expression, profiling these genes and knowing the level of expression in CRC has become quintessential to understanding the complexities involved in CRC progression (Marisa et al., 2013). Through the use of high-throughput techniques like micro array analysis, we can identify differentially expressed genes that have key roles to play in the development of primary tumors and metastasis. These will also aid in providing insights into the molecular mechanisms of CRC and serve as biomarkers for diagnosis, prognosis and therapeutic targeting (Gerhold et al., 2002).  
Pathway enrichment analyses, such as the Kyoto Encyclopedia of Genes and Genomes (KEGG) and Gene Set Enrichment Analysis(GSEA) help provide additional insights into understanding and revealing the biological pathways and processes significantly altered in CRC(Kanehisa et al., 2017). This will be very relevant in providing potential avenues for therapeutic interventions. 
In addition to the exploring the differentially expressed genes, survival analysis which remains an important aspect of understanding the clinical relevance of identified genes will be zoomed in on. By correlating gene expression patterns with patient outcomes, we can identify genes that play a role in cancer survival trends even beyond that of CRC (Guastadisegni et al., 2010). 
METHODOLOGY 
Data Acquisition and Preprocessing 
In conducting this research, we searched for datasets in the Gene Expression Omnibus(GEO) database with the keyword colorectal cancer and decided on utilizing the GSE41258 dataset. The dataset comprises gene expression profiles from patients with colonic neoplasms who presented at a medical center between 1992 and 2004. It includes 390 samples representing 4 tissue types namely, primary colon adeno carcinomas, adenomas, metastasis and corresponding normal mucosae. 
The dataset was loaded using the biobase and GQuery packages in R. The expression matrix was extracted and the phenotypic data was filtered to exclude those with missing annotations (Klaus & Reisenauer, 2018). 
 
Differential Gene Expression Analysis 
To identify the differentially expressed genes between tissue types, we performed a differential gene expression analysis with the limma package (Ritchie et al., 2015).  Tissue samples were grouped based on tissue type and a design matrix was created to model comparisons. After modelling eBayes was used to compute the model t-statistics. The most significant genes for each tissue pair were obtained and a Benjamini-Hochberg correction added to adjust the p-values obtained. Comparisons between the Normal Colon and Primary turmor, Normal liver and liver metastasis, normal lung and lung metastasis and finally for polyps and high-grade polyps were made for this research but only that of the normal colon and primary tumor was discussed in this article. 
Gene Annotation and Identification 
Gene annotations of the differentially expressed probe IDs with their gene symbols were made with the hguplus2.db package. This also came with their biological information allowing for a mapping of probe IDs to HGNC symbols (Braschi et al., 2019) and other identifiers such as Ensembl and EntrezIDs(Aken et al., 2016). 
Pathway Enrichment Analysis (KEGG) and Gene Set Enrichment Analysis 
A KEGG pathway enrichment analysis was then performed using the clusterProfiler package. The differentially expressed genes were mapped to KEGG pathways and visualized to identify the overly represented biological processes (Kanehisa et al., 2017). 
GSEA was also conducted on the significantly expressed genes and ranked based on their log fold changes. This helps in the identification of gene sets significantly associated with colorectal cancer progression. 
After enrichment analysis was conducted and significantly expressed genes identified, survival analysis was performed to evaluate the impact of gene expression on patient outcomes. Clinical data was obtained from The Cancer Genome Atlas (TCGA) (Colaprico et al., 2016) and survival information such as “days to last follow-up” and “days to death” along with their “vital status” was used in conjunction with the expression data for the survival analysis.  
A univariate Cox proportional hazard regression model was developed and applied using the survival package in R to evaluate the relationship between the expression of individual genes and patient survival. The Kaplan Meir method which has been proven to be definitive in clinical survival analysis was employed based on the high or low expression of the differentially expressed genes(Benítez-Parejo et al., 2011). 
 
RESULTS 
Differentially Expressed Genes 
There were 22283 genes in the dataset used for the differential gene expression between the normal colon and primary tumor. Out of those, 10081 were differentially expressed and statistically significant. The top 15 most significantly upregulated and the top 9 genes most significantly downregulated are provided below. 
Upregulated genes 
Gene 	probe_id 	logFC 	P.Value 	adj.P.Val 
PLG 	209978_s_at 	72369.43 	1.14E-24 	9.72E-23 
TFR2 	210215_at 	48129.57 	1.20E-22 	7.21E-21 
SERPINA10 	220626_at 	30572.89 	5.46E-26 	6.11E-24 
NAT8 	210289_at 	7614.238 	8.87E-50 	1.98E-45 
ANGPTL3 	219803_at 	5048.767 	3.54E-28 	6.58E-26 
CYP3A4 	208367_x_at 	4435.72 	4.34E-26 	4.99E-24 
SAA4 	207096_at 	4309.519 	1.14E-26 	1.53E-24 
CYP2C9 	216025_x_at 	3894.575 	3.62E-26 	4.29E-24 
DPYS 	206065_s_at 	3785.595 	8.74E-43 	3.24E-39 
TF 	214063_s_at 	3666.145 	9.20E-24 	6.77E-22 
F9 	207218_at 	3631.417 	9.72E-29 	2.10E-26 
HPR 	208471_at 	3580.712 	5.79E-33 	3.58E-30 
CYP2C9 	220017_x_at 	3542.987 	8.48E-26 	8.83E-24 
CPS1 	204920_at 	3447.76 	1.63E-26 	2.12E-24 
Table1.0 
Downregulated genes  
Gene 	probe_id 	logFC 	P.Value 	adj.P.Val 
BTC 	207326_at 	-10.3941 	3.19E-12 	4.59E-11 
BBS9 	209958_s_at 	-15.3672 	5.94E-14 	1.10E-12 
SEC24D 	215641_at 	-16.909 	1.90E-14 	3.73E-13 
SLC25A12 	203339_at 	-19.5637 	5.73E-12 	7.91E-11 
CSHL1 	208295_x_at 	-24.4496 	2.18E-12 	3.20E-11 
TMC7 	220021_at 	-24.9164 	6.39E-14 	1.18E-12 
SMYD5 	209516_at 	-30.1066 	5.21E-12 	7.25E-11 
LPAR2 	206723_s_at 	-32.0751 	9.84E-14 	1.76E-12 
DCHS2 	220373_at 	-35.6755 	2.93E-15 	6.47E-14 
Table 1.1 
 
 
GO biological Process Enrichment 
The Gene Ontology analysis of the differentially expressed genes found several biological, molecular and cellular processes that were significantly enriched and these included glucuronate and uronic acid metabolic processes, cellular glucuronidation and biological processes involving a response to xenobiotics.  
 

GSEA Enrichment Analysis 
In the GSEA ,significant enrichment in immune-related and metabolic processes was observed.  Top enriched processes included myeloid cell activation, leucocyte migration inflammatory response. Metabolic pathways such as the monocarboxylic acid metabolism, lipid localization and small molecule biosynthesis were significantly enriched. 
 
KEGG Pathway Enrichment 
Applying the KEGG pathway analysis revealed significantly enriched pathways among the differentially expressed genes, These pathways included fatty acid  degradation, tyrosine metabolism, glycolysis/gluconeogenesis porphyrin metabolism and steroid hormone biosynthesis. 
 

 
 
Survival Analysis 
Survival analysis was conducted on using the most significantly expressed genes on a TCGA-COAD clinical dataset. Among the analyzed genes several showed statistically significant associations with patient survival. Some of these genes include the C8G gene, the OCT1 gene, TFR2, LP(A) genes. All of these genes exhibited poorer survival rates in the higher expression groups compared to the lower expression groups whiles others like PLG had little to no difference in both expression groups.
 
 
Cox ph model was also conducted and the results displayed . 
The Cox proportional hazards regression analysis for the various genes presents a mix of statistically significant and non-significant results. Here's a brief breakdown: 
  
GENE ID	HAZARD RATIO (HR)	P-VALUE	SIGNIFICANCE
ENSG00000122194	0.99	0.972	Not Significant
ENSG00000129988	1.88	0.0228	Significant
ENSG00000113889	1.07	0.799	Not Significant
ENSG00000198670	0.71	0.197	Not Significant
ENSG00000158104	0.70	0.19	Not Significant
ENSG00000106327	0.79	0.388	Not Significant
ENSG00000187048	0.93	0.779	Not Significant
ENSG00000162365	0.82	0.485	Not Significant
ENSG00000136881	0.98	0.937	Not Significant
ENSG00000135447	0.62	0.0737	Approaching Significance
ENSG00000140093	1.18	0.535	Not Significant
ENSG00000109758	0.77	0.338	Not Significant
ENSG00000163581	1.95	0.0473	Significant
ENSG00000176919	0.61	0.0704	Approaching Significance
ENSG00000118271	0.97	0.922	Not Significant

  Table 2.0
 
GENE ID	COEFFICIENT	EXP(COEF)	SE(COEF)	Z-VALUE	P-VALUE	CONCOEDANCE	LIKELIHOOD RATIO	WALD TEST	LOGRANK SCORE
ENSG00000122194	-0.0096	0.9904	0.2723	-0.035	0.972	0.492	0	0	0
ENSG00000129988	0.6324	1.8821	0.2777	2.277	0.0228	0.553	5.44	5.18	5.36
ENSG00000113889

	0.0676	1.0699	0.2653	0.255	0.799	0.510	0.06	0.06	0.06
ENSG00000198670

	-0.3430	0.7096	0.2656	-1.292	0.197	0.533	1.67	1.67	1.68
ENSG00000158104	-0.3525	0.7030	0.2692	-1.309	0.19	0.553	1.73	1.71	1.73
ENSG00000106327	-0.2331	0.7921	0.2702	-0.863	0.388	0.536	0.75	0.74	0.75
ENSG00000187048	-0.0748	0.9279	0.2662	-0.281	0.779	0.497	0.08	0.08	0.08
ENSG00000162365	-0.1993	0.8193	0.2854	-0.699	0.485	0.533	0.48	0.49	0.49
ENSG00000136881	-0.0213	0.9789	-0.08	0.937	0.512	0.512	0.01	0.01	0.01
ENSG00000135447	-0.4809	0.6182	0.2689	-1.788	0.0737	0.564	3.25	3.2	3.26
ENSG00000140093	0.1654	1.1798	0.2663	0.621	0.535	0.533	0.39	0.39	0.39
Table 2.1
DISCUSSION 
In this exploratory research utilizing GSE41258, among the top differentially expressed and upregulated genes were PLG,TFR2,and SERPINA10, all of which are implicated in CRC progression and has been suggested by a number of research papers to be vital in determining the prognosis of CRC and determining further therapies for patients (Seetoo et al., 2003) .The overexpression of PLG(plasminogen) is particularly important as it plays a crucial role in fibrinolysis, tumor cell invasion and metastasis by enhancing the breakdown if the extracellular matrix(ECM) and facilitating angiogenesis.  The overexpression of this gene is known to be an important driver of metastasis in various cancers.
TFR2(Transferin Receptor 2) which plays a key role in the regulation of iron homeostatis was also significantly upregulated. Dysregulated iron homeostasis in a known key factor in cancer biology (Calzolari et al., 2007)  and has been seen not only in CRC but in most cancer cell lines (CALZOLARI et al., 2009) where increases iron uptake supports the proliferation of cancer cells by enhancing oxidative stress. Enrichment of iron-related pathways further supports and emphasizes the role of metabolic changes in CRC tumor progression.
SERPINA10 which is a member of the Serpin family is involved in coagulation and has been explored as a useful tool in breast cancer prognosis (YU et al.,2024), making it unsurprising that it as found to be differentially expressed in CRC. This goes on to suggest that coagulation-related factors contribute to tumorigenesis by fostering a pro-thrombic environment that supports cancer progression.
In addition to the significantly upregulated genes, there were others that were found to be downregulated lice the BTC(Betacellulin) gene. This gene is involved in epithelial cell growth and differentiation. Several mechanisms have been proposed for its role in different kinds of cancer (Yokoyama et al. 1995) however, the downregulation in CRC may indicate that the inability of cells to differentiate, indicating an impaired growth factor signaling causes a reflex proliferation of nascent cells affecting normal cellular proliferation. Similarly, the downregulation of BBS9 and SEC24D all indicate irregular cell signaling and protein transport important for cellular homeostasis in non-cancerous tissues.
The KEGG pathway analysis revealed metabolic pathways such as fatty acid degradation, glycolysis/gluconeogenesis and steroid hormone biosynthesis as being significantly involved in CRC. The energy production and synthetic pathways are vital for sustaining the high metabolic demands of tumor cells and this is clearly understood as the Warburg effect in cancer literature Chen et al., 2007; Wu & Zhao, 2013) . 
In the GO biological process analysis, cellular glucuronidation, response to xenobiotics and processes involved in detoxification and metabolic reprogramming were enriched. These processes maybe be enriched because tumor cells are trying to overcome the unsuitable microenvironment, they find themselves in or maybe trying to survive chemotherapy and other therapeutic interventions in which most primary sources of the tissues used to obtain the dataset were undergoing. Glucuronidation in particular has been identified as a useful protective mechanism for tumor cells (Cummings et al., 2004).
In the GSEA some previous findings were reinforced. Gene sets related to immune response and metabolic processes were enriched. Gene sets related to immune-related pathways such as myeloid cell activation and leucocyte migration were highly enriched, highlighting the role immune cell infiltration plays in CRC (Ge et al., 2019). This adds up to findings suggesting that immune evasion and inflammation are central to CRC progression (Grasso et al., 2018; Schmitt & Greten, 2021). Lipid localization, monocarboxylic acid metabolism and small molecule biosynthesis were all significantly enriched providing evidence to what was previously found in the GO about CRC cells undergoing metabolic programming to ensure sustained growth and tumor survival.
To investigate the clinical significance of the differentially expressed genes survival analysis was conducted. The Cox proportional hazard regression revealed that TFR2(HR=1.88,P=0.0228) and SERPINA10(HR=1.95,p=0.0473) were significantly associated with poor patient outcomes. This indicates that their overexpression contributes to an increase in the hazard or danger posed by CRC and serves as an unfavorable prognostic marker manifesting in aggressive disease phenotypes. The value for BBS9(HR=0.61,P=0.0704) which was approaching significance was indicative that the gene may have protective effects. This makes sense given that it was significantly downregulated. PLG on the other hand was a bit conflicting given that it was over expressed but with an HR=0.62, it indicated that it sis not increase the hazard risk in CRC. This may be due to the fact that other genes that assist in its fuction were not found to be significantly overexpressed in CRC and hence did not have any impact on the survival of patients in their absence. 
This study is in no way exhaustive and has several limitations. The use of microarray data offers a relatively lower insight as compared to data obtained from RNA-sequencing. Also this analysis was based on transcriptomic data solely; integrating additional layers of molecular data such as DNA methylation or proteomic profiles could provide a more comprehensive view of CRC’s molecular landscape. In vivo studies need to be conducted to validate the findings of this study
Conclusion 
In summary we identified several key differentially expressed genes and enriched pathways in colorectal cancer based on the GSE41258 microarray, that provided insights into the molecular mechanisms possibly underlying CRC progression and patient survival. We were also able to identify and associate certain genes with their impact on survival outcomes in high and low expression groups as well as highlight new potential therapeutic targets and prognostic biomarkers.  
 
REFERENCES
 
 
1.	 Aken, B.L. et al. (2016) ‘The Ensembl gene annotation system’, Database, 2016, p. baw093. Available at: https://doi.org/10.1093/database/baw093. 
2.	Benítez-Parejo, N., Rodríguez del Águila, M.M. and Pérez-Vicente, S. (2011) ‘Survival analysis and Cox regression’, Allergologia et Immunopathologia, 39(6), pp. 362–373. Available at: https://doi.org/10.1016/j.aller.2011.07.007.
3.	 Braschi, B. et al. (2019) ‘Genenames.org: the HGNC and VGNC resources in 2019’, Nucleic Acids Research, 47(D1), pp. D786–D792. Available at: https://doi.org/10.1093/nar/gky930.
4.	 Calzolari, A. et al. (2007) ‘Transferrin receptor 2 is frequently expressed in human cancer cell lines’, Blood Cells, Molecules, and Diseases, 39(1), pp. 82–91. Available at: https://doi.org/10.1016/j.bcmd.2007.02.003. 
5.	CALZOLARI, A. et al. (2009) ‘Regulation of transferrin receptor 2 in human cancer cell lines’, Blood Cells, Molecules, and Diseases, 42(1), pp. 5–13. Available at: https://doi.org/10.1016/j.bcmd.2008.10.001.
6.	 Chen, Z. et al. (2007) ‘The Warburg effect and its cancer therapeutic implications’, Journal of Bioenergetics and Biomembranes, 39(3), pp. 267–274. Available at: https://doi.org/10.1007/s10863-007-9086-x. 
7.	Colaprico, A. et al. (2016) ‘TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data’, Nucleic Acids Research, 44(8), pp. e71–e71. Available at: https://doi.org/10.1093/nar/gkv1507. 
8.	Cummings, J. et al. (2004) ‘Glucuronidation as a mechanism of intrinsic drug resistance in colon cancer cells: contribution of drug transport proteins’, Biochemical Pharmacology, 67(1), pp. 31–39. Available at: https://doi.org/10.1016/j.bcp.2003.07.019. 
9.	Fletcher, R. et al. (2018) ‘Colorectal cancer prevention: Immune modulation taking the stage’, Biochimica et Biophysica Acta (BBA) - Reviews on Cancer, 1869(2), pp. 138–148. Available at: https://doi.org/10.1016/j.bbcan.2017.12.002. 
10.	Gerhold, D.L., Jensen, R. V. and Gullans, S.R. (2002) ‘Better therapeutics through microarrays’, Nature Genetics, 32(S4), pp. 547–552. Available at: https://doi.org/10.1038/ng1042. 
11.	Grasso, C.S. et al. (2018) ‘Genetic Mechanisms of Immune Evasion in Colorectal Cancer’, Cancer Discovery, 8(6), pp. 730–749. Available at: https://doi.org/10.1158/2159-8290.CD-17-1327.
12.	Guastadisegni, C. et al. (2010) ‘Microsatellite instability as a marker of prognosis and response to therapy: A meta-analysis of colorectal cancer survival data’, European Journal of Cancer, 46(15), pp. 2788–2798. Available at: https://doi.org/10.1016/j.ejca.2010.05.009. 
13.	Guglietta, S. and Rescigno, M. (2016) ‘Hypercoagulation and complement: Connected players in tumor development and metastases’, Seminars in Immunology, 28(6), pp. 578–586. Available at: https://doi.org/10.1016/j.smim.2016.10.011. 
14.	Kanehisa, M. et al. (2017) ‘KEGG: new perspectives on genomes, pathways, diseases and drugs’, Nucleic Acids Research, 45(D1), pp. D353–D361. Available at: https://doi.org/10.1093/nar/gkw1092. 
15.	Klaus, B. and Reisenauer, S. (2018) ‘An end to end workflow for differential gene expression using Affymetrix microarrays’, F1000Research, 5, p. 1384. Available at: https://doi.org/10.12688/f1000research.8967.2. 
16.	Marisa, L. et al. (2013) ‘Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value’, PLoS Medicine, 10(5), p. e1001453. Available at: https://doi.org/10.1371/journal.pmed.1001453. 
17.	Peltier, J. et al. (2016) ‘Quantitative proteomic analysis exploring progression of colorectal cancer: Modulation of the serpin family’, Journal of Proteomics, 148, pp. 139–148. Available at: https://doi.org/10.1016/j.jprot.2016.07.031. 
18.	Ritchie, M.E. et al. (2015) ‘limma powers differential expression analyses for RNA-sequencing and microarray studies’, Nucleic Acids Research, 43(7), pp. e47–e47. Available at: https://doi.org/10.1093/nar/gkv007. 
19.	Roshandel, G., Ghasemi-Kebria, F. and Malekzadeh, R. (2024) ‘Colorectal Cancer: Epidemiology, Risk Factors, and Prevention’, Cancers, 16(8), p. 1530. Available at: https://doi.org/10.3390/cancers16081530.
20.	 Schmitt, M. and Greten, F.R. (2021) ‘The inflammatory pathogenesis of colorectal cancer’, Nature Reviews Immunology, 21(10), pp. 653–667. Available at: https://doi.org/10.1038/s41577-021-00534-x.
21.	 Seetoo, D. et al. (2003) ‘Quantitative expression of protein markers of plasminogen activation system in prognosis of colorectal cancer’, Journal of Surgical Oncology, 82(3), pp. 184–193. Available at: https://doi.org/10.1002/jso.10210.
22.	 Simons, C.C.J.M. et al. (2013) ‘A novel classification of colorectal tumors based on microsatellite instability, the CpG island methylator phenotype and chromosomal instability: Implications for prognosis’, Annals of Oncology, 24(8), pp. 2048–2056. Available at: https://doi.org/10.1093/annonc/mdt076. 
23.	Song, K. et al. (2019) ‘Identification of genes with universally upregulated or downregulated expressions in colorectal cancer’, Journal of Gastroenterology and Hepatology, 34(5), pp. 880–889. Available at: https://doi.org/10.1111/jgh.14529. 
24.	Stepchenko, A.G. et al. (2022) ‘Suppression of OCT-1 in Metastatic Breast Cancer Cells Reduces Tumor Metastatic Potential, Hypoxia Resistance, and Drug Resistance’, Life, 12(9), p. 1435. Available at: https://doi.org/10.3390/life12091435. 
25.	Wang, Guanglin et al. (2022) ‘Uncovering potential genes in colorectal cancer based on integrated and DNA methylation analysis in the gene expression omnibus database’, BMC Cancer, 22(1). Available at: https://doi.org/10.1186/s12885-022-09185-0.
26.	 Willier, S., Butt, E. and Grunewald, T.G.P. (2013) ‘Lysophosphatidic acid (LPA) signalling in cell migration and cancer invasion: A focussed review and analysis of LPA receptor gene expression on the basis of more than 1700 cancer microarrays’, Biology of the Cell, 105(8), pp. 317–333. Available at: https://doi.org/10.1111/boc.201300011.
27.	 Wu, W. and Zhao, S. (2013) ‘Metabolic changes in cancer: beyond the Warburg effect’, Acta Biochimica et Biophysica Sinica, 45(1), pp. 18–26. Available at: https://doi.org/10.1093/abbs/gms104.  
 
 
 
 
 
 
 

