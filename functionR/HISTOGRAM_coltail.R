# Nuc Positionning prof
# Cuvier''s Team
# Depierre*
# 2018 16 avril

# plot HISTO with colored tail corresponding to a given percent

print('USAGE : plotHISTO_coltail(HIST_data = vec(), THqtl = 0.05, VarName = "name", info = "info", outDir = paste0(workdir, "David_M2_K9K27DE/integration_DE_DHMG/HISTOGRAM/"))')

plotHISTO_coltail = function(HIST_data, THqtl = 0.05, VarName = NULL, info = NULL, outDir = paste0(workdir, "David_M2_K9K27DE/integration_DE_DHMG/HISTOGRAM/")){
  histgm <- hist(HIST_data, breaks=100, plot=FALSE)
  mycol = ifelse(histgm$breaks < quantile(HIST_data, THqtl), "red", ifelse (histgm$breaks > quantile(HIST_data, 1-THqtl), "blue", "gray50"))
  v1 = round(quantile(HIST_data, THqtl), 2)
  v2 = round(quantile(HIST_data, 1-THqtl), 2)
  pdf(paste0(outDir, "HISTOGRAM_",VarName ,"_", THqtl,"tail_",info,".pdf"), width = 16, height = 10 )
  plot(histgm, col = mycol, border=T, main=paste0("Histogram of ", VarName, " / red > ",v1," / blue < ",v2," / tail=",THqtl))
  dev.off()
}
