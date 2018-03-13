# This function runs a predictor (defined by "command") and returns it's performance (AUC,concordance as specified)
# y is a vector of known labels
run_predictor<-function(y,command,args,repeats = 50,keep_concordance = F)
{
	rocs = c()
	concordance = c()
	for(i in 1:repeats)
	{
		args[["b"]] = sample(c(1:10),nrow(args[[1]]),replace=T)
		cv_results = do.call(command,args)
		if(keep_concordance)
		{
			concordance = c(concordance,cv_results[length(cv_results)])
			cv_results = cv_results[-length(cv_results)]
		} else { concordance = c(concordance,NA) }
		fg = cv_results[y == 1]; bg = cv_results[y == 0]
		roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
		rocs = c(rocs,roc$auc)
	}
	return(rbind(rocs,concordance))
}
