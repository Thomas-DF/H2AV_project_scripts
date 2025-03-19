# Heterochromatin k9k27 variation project
# DEPIERRE 2018

#
# print("Function that takes as input list of profiles (list of matrix) and a list of genes names and return list of average profil (list of vector) only for the selected genes ------------------
# 		 2 STEP : 1/Selection select_prof() 2/Normalisation NORMselect_prof()")
# # Function that takes as input list of profiles (list of matrix) and a list of genes names and return list of average profil (list of vector) only for the selected genes
# select_prof = function(prof_list, gene_names, avg=FALSE, Nsplit=0){
# 	if(Nsplit != 0){
# 			gene_names = split(gene_names, ceiling(seq_along(gene_names)/ceiling(length(gene_names)/Nsplit)))
# 			prof_list_selctd = lapply(prof_list, function(prof){
# 				lapply(gene_names, function(gene_names_sub){prof = prof[rownames(prof) %in% gene_names_sub,]})
# 			})
# 			if(avg==FALSE){
# 				return(prof_list_selctd)
# 			}else{
# 				avg_prof_list = lapply(prof_list_selctd, function(prof_by_qrtl){
# 					lapply(prof_by_qrtl,apply,2,mean,na.rm=T)
# 				})
# 				return(avg_prof_list)
# 			}
# 	}else{
# 		prof_list_selctd = lapply(prof_list, function(prof){prof = prof[rownames(prof) %in% gene_names,]})
# 		if(avg==FALSE){
# 			return(prof_list_selctd)
# 		}else{
# 			avg_prof_list = lapply(prof_list_selctd,apply,2,mean,na.rm=T)
# 			return(avg_prof_list)
# 		}
# 	}
# }


print("Function that takes as input list of profiles (list of matrix) and a list of genes names and return list of average profil (list of vector) only for the selected genes ------------------
		 2 STEP : 1/Selection select_prof() 2/Normalisation NORMselect_prof()")
# Function that takes as input list of profiles (list of matrix) and a list of genes names and return list of average profil (list of vector) only for the selected genes
select_prof = function(prof_list, gene_names, avg=FALSE, Nsplit=0){
	if(Nsplit != 0){
			gene_names = split(gene_names, ceiling(seq_along(gene_names)/ceiling(length(gene_names)/Nsplit)))
			prof_list_selctd = lapply(prof_list, function(prof){
				lapply(gene_names, function(gene_names_sub){prof = prof[rownames(prof) %in% gene_names_sub,]})
			})
			if(avg==FALSE){
				return(prof_list_selctd)
			}
			if(avg=="mean"){
				avg_prof_list = lapply(prof_list_selctd, function(prof_by_qrtl){
					lapply(prof_by_qrtl,apply,2,mean,na.rm=T)
				})
				return(avg_prof_list)
			}
			if(avg=="median"){
				avg_prof_list = lapply(prof_list_selctd, function(prof_by_qrtl){
					lapply(prof_by_qrtl,apply,2,median,na.rm=T)
				})
				return(avg_prof_list)
			}
	}else{
		prof_list_selctd = lapply(prof_list, function(prof){prof = prof[rownames(prof) %in% gene_names,]})
		if(avg==FALSE){
			return(prof_list_selctd)
		}else{
			avg_prof_list = lapply(prof_list_selctd,apply,2,mean,na.rm=T)
			return(avg_prof_list)
		}
	}
}




NORMselect_prof_distalTSS = function(output_select_prof, start=1, stop=1000){
	NormFact = mean(output_select_prof[[1]][[1]][start:stop])
	lapply(output_select_prof, function(prof_list){
		lapply(prof_list, function(prof){
			prof/mean(prof[start:stop])*NormFact
		})
	})
}


NORMselect_prof = function(output_select_prof){
	NormFact = mean(output_select_prof[[1]][[1]])
	lapply(output_select_prof, function(prof_list){
		lapply(prof_list, function(prof){
			prof/mean(prof)*NormFact
		})
	})
}
