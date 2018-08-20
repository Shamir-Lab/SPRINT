# This function runs a predictor (defined by "command") and returns it's performance (AUC,concordance as specified)
# y is a vector of known labels
run_predictor<-function(y,command,args,repeats = 50)
{
	rocs = c()
	for(i in 1:repeats)
	{
		cv_results = do.call(command,args)
		if(class(y)=="matrix")
		{
			fg = cv_results[y[,"status"] == 1]; bg = cv_results[y[,"status"] == 0]
		} else {
			fg = cv_results[y == 1]; bg = cv_results[y == 0]
		}
		roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
		rocs = c(rocs,roc$auc)
	}
	return(rbind(rocs))
}
