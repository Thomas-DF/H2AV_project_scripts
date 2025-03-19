

# Heterochromatin spread / Single Cell RNA seq project/ Nuc_Positioning
# Cuvier''s Team
# Depierre*
# 2018


print("USAGE :     List_NOOUTLIER = removeOutliers(list_MTX = List_1COLmtx, THpct = 0.05)")
print("USAGE :     PROFList_NOOUTLIER = removeOutliers(list_prof = List_PROF, THpct = 0.05) -> ONLY LOW PROF ARE REMOVED")

removeOutliers = function(list_MTX, THpct = 0.05){
	list_MTX = lapply(list_MTX, function(MTX){
		if(ncol(MTX)==1){
			MTX.temp = MTX[MTX[,1] >= quantile(MTX, THpct),,drop=F]
			MTX.temp = MTX.temp[MTX.temp[,1] <= quantile(MTX, 1-THpct),,drop=F]
			return(MTX.temp)
		}
	})
	list_MTX_gn = lapply(list_MTX, rownames)
	comRownames_NOoutlier = Reduce(intersect,list_MTX_gn)
	list_MTX = lapply(list_MTX, function(feat1){feat1 = feat1[rownames(feat1) %in% comRownames_NOoutlier,,drop=F]})
	return(list_MTX)
}

removeOutliersProf = function(list_prof,THpct = 0.05 ){
		list_prof = lapply(list_prof, function(prof){
		if(ncol(prof)>1000){
			sum_prof = apply(prof, 1, sum)
			prof = prof[sum_prof >= quantile(sum_prof, THpct),,drop=F]
		}
		return(prof)
	})
	list_prof_gn = lapply(list_prof, rownames)
	comRownames_NOoutlier = Reduce(intersect,list_prof_gn)
	list_prof = lapply(list_prof, function(feat1){feat1 = feat1[rownames(feat1) %in% comRownames_NOoutlier,,drop=F]})
	return(list_prof)

}
