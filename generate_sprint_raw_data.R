# Preprocess of SPRINT data. All files must be in the currect working directory
generate_raw_data<-function(){
	raw_baseline = read.csv('./baseline.csv',row.names=1,header=T)
	raw_outcomes = read.csv('./outcomes.csv',row.names=1,header=T)
	raw_safety = read.csv('./safety.csv',row.names=1,header=T)
	raw_data = merge(raw_baseline,raw_outcomes,by="row.names")
	rownames(raw_data) = raw_data[,1]
	raw_data = raw_data[,2:58]
	raw_data = merge(raw_data,raw_safety,by="row.names")
	rownames(raw_data) = raw_data[,1]
	raw_data = raw_data[,2:98]
	raw_data = race4_to_numeric(raw_data)
	return(raw_data)
}

race4_to_numeric<-function(x)
{
	race_matrix = matrix(rep(0,nrow(x)*4),nrow = nrow(x), ncol=4)
	colnames(race_matrix) = c("BLACK","WHITE","HISPANIC","OTHER")
	for(i in 1:nrow(x))
	{
		race_matrix[i,x[i,"RACE4"]] = 1
	}
	colnames(race_matrix) = c("ETHNIC_BLACK","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER")
	race_feature_position = which(colnames(x) == "RACE4")
	x = cbind(x[,1:(race_feature_position-1)],race_matrix[,2:4],x[,(race_feature_position+1):ncol(x)])
	return(x)
}
