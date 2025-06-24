workdir = "~/Bureau/tdefreitas_genobioinfo/PROJET_H2AV_2025/"

LIST_QUANTIF_K36=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF_K36.RDS"))
ZSCORE_H3K36me3_WT3_GB_N_WT = LIST_QUANTIF_K36$ZSCORE_H3K36me3_WT3_GB_N_WT
ZSCORE_H3K36me3_SEA4_GB_N_WT = LIST_QUANTIF_K36$ZSCORE_H3K36me3_SEA4_GB_N_WT

LIST_QUANTIF=readRDS(paste0(workdir,"DATA/LIST_FEATURES/LIST_QUANTIF.RDS"))
ZSCORE_H2AV_WT3_GB_N_WT = LIST_QUANTIF$ZSCORE_H2AV_WT3_GB_N_WT

ZSCORE_RAD51_WT3_GB_N_WT = LIST_QUANTIF$ZSCORE_RAD51_WT3_GB_N_WT

ZSCORE_H2AV_GB_WT_N = LIST_QUANTIF$ZSCORE_H2AV_GB_WT_N

getNameList = function(Vec, topdown = "top", prct = 10){
  Vec = Vec[order(Vec, decreasing=T)]
  if(topdown %in% "top"){
    GN = names(Vec[Vec > quantile(Vec, (100-prct)/100)])
  }
  if(topdown %in% "down"){
    GN = names(Vec[Vec < quantile(Vec, (prct)/100)])
  }
  if(topdown %in% "mid"){
    tmp1 = names(Vec[Vec < quantile(Vec, (100/2-prct/2)/100)])
    tmp2 = names(Vec[Vec < quantile(Vec, (100/2-prct/2+prct)/100)])
    GN = tmp2[tmp2 %ni% tmp1]
  }
  return(GN)
}

LIST_ZSCORE_H2AV_WT3_GB_N_WT = getNameList(ZSCORE_H2AV_WT3_GB_N_WT, topdown = "top", prct = 10)

LIST_ZSCORE_H3K36me3_WT3_GB_N_WT = getNameList(ZSCORE_H3K36me3_WT3_GB_N_WT, topdown = "down", prct = 10)

# LIST_ZSCORE_H3K36me3_SEA4_GB_N_WT = getNameList(ZSCORE_H3K36me3_SEA4_GB_N_WT, topdown = "top", prct = 25)

LIST_ZSCORE_RAD51_WT3_GB_N_WT= getNameList(ZSCORE_RAD51_WT3_GB_N_WT, topdown = "down", prct = 10)

Reduce(intersect, list(LIST_ZSCORE_H2AV_WT3_GB_N_WT,LIST_ZSCORE_H3K36me3_WT3_GB_N_WT,LIST_ZSCORE_RAD51_WT3_GB_N_WT))

LIST_ZSCORE_H2AV_GB_WT_N = getNameList(ZSCORE_H2AV_GB_WT_N, topdown = "down", prct = 5)
LIST_ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse = getNameList(ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse, topdown = "top", prct = 5)
LIST_ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse = getNameList(ZSCORE_PROFMAT_Q_DRIP_HLuc_vs_HLuc_RNAse, topdown = "down", prct = 75)

Reduce(intersect, list(LIST_ZSCORE_H2AV_GB_WT_N,LIST_ZSCORE_PROFMAT_Q_DRIP_Nelf_vs_Nelf_RNAse))

# Reduce(intersect, list(LIST_ZSCORE_H2AV_WT3_GB_N_WT,LIST_ZSCORE_H3K36me3_WT3_GB_N_WT,LIST_ZSCORE_RAD51_WT3_GB_N_WT))
# [1] "FBgn0035147.1" "FBgn0040235.1" "FBgn0031022.1" "FBgn0035153.1" "FBgn0033926.1" "FBgn0037610.1" "FBgn0264267.1" "FBgn0035140.1" "FBgn0035149.1" "FBgn0037653.1" "FBgn0039691.1"
# [12] "FBgn0033813.1" "FBgn0028687.1" "FBgn0011787.1" "FBgn0028692.1" "FBgn0003062.1" "FBgn0267778.1" "FBgn0020626.1"

