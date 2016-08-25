
fileName = '~/Documents/Public-Datasets/2016-8-24_Fatty-Liver-Disease_Gene-Expression_Metabolomics/Original-Data/MetaboLights_ MTBLS174/MTBLS174_20160823_115437/m_hna_fld_metabolite_profiling_NMR_spectroscopy_v2_maf.tsv'
metabolomics <- read.delim(file = fileName, sep = "\t", stringsAsFactors = FALSE)

fileName = '~/Documents/Public-Datasets/2016-8-24_Fatty-Liver-Disease_Gene-Expression_Metabolomics/Processed-Data/MetaboLights/m_hna_fld_metabolite_profiling_NMR_spectroscopy_v2_maf.txt'
write.table(metabolomics, fileName, sep = "\t", row.names = FALSE)