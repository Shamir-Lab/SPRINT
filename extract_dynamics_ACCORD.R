# Prepare ACCORD matrix to include dynamic features
preprocess_accord<-function(x,t,bp,dynamics_extraction_method)
{
	bp = bp[,c(1:4)]
	bp = bp[-which(bp[,2]=="EXIT"),]
	bp[,2] = substr(bp[,2],2,3)
	bp[bp[,2]=="LR",2] = 0
	bp_2 = matrix(c(as.numeric(bp[,1]),as.numeric(bp[,2]),as.numeric(bp[,3]),as.numeric(bp[,4])),nrow=nrow(bp),ncol=4)
	bp = bp_2
	colnames(bp) = c("p_name","VISITCODE","SBP","DBP")
	x_m = x
	x_m = cbind(x_m,matrix(rep(0,nrow(x_m)*5),ncol=5))
	colnames(x_m)[(ncol(x_m)-4):ncol(x_m)] = c("pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")
	for(i in 1:nrow(x_m))
	{	
		p = rownames(x_m)[i]
		val = dynamics_extraction_method(p,x_m,bp,t)
		x_m[p,"pp_m"] = val[1]
		x_m[p,"pp_mean"] = val[2]
		x_m[p,"pp_R2"] = val[3]
		x_m[p,"pp_f_pvalue"] = val[4]	
		x_m[p,"pp_max_min"] = val[5]	
	}
	x_m = x_m[which(!is.na(x_m[,"pp_m"])),]
	x_m = x_m[which(!is.na(x_m[,"pp_mean"])),]
	x_m = x_m[which(!is.na(x_m[,"pp_R2"])),]
	x_m = x_m[which(!is.na(x_m[,"pp_f_pvalue"])),]
	return(x_m)
}

# Extract dynamic features from longitudinal data for patient p
calculate_cox_accord<-function(p,x,bp_data,t)
{
	pp = c()
	timepoints = c()
	p_bp = bp_data[which(bp_data[,1] == p),]
	if(is.null(nrow(p_bp)) || nrow(p_bp) < 3)
	{
		return(c(NA,NA,NA,NA,NA))
	}
	p_bp = p_bp[order(p_bp[,"VISITCODE"]),]
	if((p_bp[3,"VISITCODE"] + t)*30 > x[p,"fuyrs_po"])
	{
		return(c(NA,NA,NA,NA,NA))
	}
	i = 1
	while(i <= nrow(p_bp) && (p_bp[i,"VISITCODE"] + t)*30 <= x[p,"fuyrs_po"])
	{
		pp = c(pp,p_bp[i,"SBP"] - p_bp[i,"DBP"])
		timepoints = c(timepoints,p_bp[i,"VISITCODE"])
		i = i + 1
	}
	pp_ols = lm(pp~timepoints)
	pp_m = pp_ols$coefficients["timepoints"]
	pp_mean = mean(pp)
	pp_R2 = summary(pp_ols)$r.squared
	if(!(pp_R2 == 0) && !is.na(pp_R2)){
	pp_f_pvalue = pf(summary(pp_ols)$fstatistic[1],summary(pp_ols)$fstatistic[2],summary(pp_ols)$fstatistic[3],lower.tail=F)
	} else{
	pp_f_pvalue = NA
	}
	pp_max_min = max(pp) - min(pp)

	return(c(pp_m,pp_mean,pp_R2,pp_f_pvalue,pp_max_min))
}