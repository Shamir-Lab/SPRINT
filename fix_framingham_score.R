# Fix SPRINT framingham score that was calculated wrongly
# Fixed values are as corrected by Warner et al.
fix_framingham_score<-function(x)
{
	new_framingham_score = c()
	for(i in 1:nrow(x))
	{
		if(x[i,"FEMALE"]==1)
		{
			beta = 2.32888*log(x[i,"AGE"]) + 1.20904*log(x[i,"CHR"]) - 0.70833*log(x[i,"HDL"]) + 
    				(2.82263 - (x[i,"NOAGENTS"]*0.06106))*log(x[i,"SBP"]) + 0.52873*x[i,"CURR_SMOKER"]
   			new_framingham_score = c(new_framingham_score,(1 - 0.95012^exp(beta-26.1931))*100)
		} else {
			beta = 3.06117*log(x[i,"AGE"]) + 1.12370*log(x[i,"CHR"]) - 0.93263*log(x[i,"HDL"]) + 
   				(1.99881 - (x[i,"NOAGENTS"]*0.06578))*log(x[i,"SBP"]) + 0.65451*x[i,"CURR_SMOKER"]
   			new_framingham_score = c(new_framingham_score,(1 - 0.88936^exp(beta-23.9802))*100)
		}
	}
	return(new_framingham_score)
}