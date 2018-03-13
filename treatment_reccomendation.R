# This function accepts a SPRINT information matrix (x) and assigns patients to intensive/standard treatments 
# According to thier CV/kidney risks. theta is the treatment cutoff (final predicted value > theta: intensive treatment) 
treatment_reccomendation<-function(x,theta)
{
	cv_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","SMOKE_3CAT","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
					,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
					"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO")
	aki_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","SMOKE_3CAT","ASPIRIN","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
					,"ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER","GLUR","TRR","UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","EGFR_SCREAT")
	x = cbind(x,rep(0,nrow(x)))
	colnames(x)[ncol(x)] = "new_assignment"
	x_1 = x
	x_0 = x
	x_1[,"INTENSIVE"] = 1
	x_0[,"INTENSIVE"] = 0
	y = x[,"EVENT_PRIMARY"]
	time = "T_PRIMARY"
	event = "EVENT_PRIMARY"
	list_of_predictors = cv_static_predictors
	predictors = paste(cv_static_predictors,collapse = "+")		
	preds = c()
	batches = sample(c(1:10),nrow(x),replace=T)
	for(batch in unique(batches))
	{
		test_samples = which(batches==batch)
		tr_x = x[-test_samples,];tr_y=y[-test_samples]
		te_x = x_0[test_samples,];te_y = y[test_samples]
		models=simple_sampling_based_learner(x=as.data.frame(tr_x),y=tr_y,d=NULL,max_num_samples = 300,func=get_coxph_model,reps=10,time,event,predictors)
		batch_predictions_obj = predict(models,as.data.frame(te_x),type="risk")			
		preds = c(preds,batch_predictions_obj)
	}
	preds = preds[rownames(x)]
	x = cbind(x,preds)
	colnames(x)[ncol(x)] = "risk_cv_standard"
	y = x[,"AKI_ERS_EVNT"]
	time = "AKI_ERS_DAYS"
	event = "AKI_ERS_EVNT"
	list_of_predictors = aki_static_predictors
	predictors = paste(aki_static_predictors,collapse = "+")		
	preds = c()
	batches = sample(c(1:10),nrow(x),replace=T)
	for(batch in unique(batches))
	{
		test_samples = which(batches==batch)
		tr_x = x[-test_samples,];tr_y=y[-test_samples]
		te_x = x_1[test_samples,];te_y = y[test_samples]
		models=simple_sampling_based_learner(x=as.data.frame(tr_x),y=tr_y,d=NULL,max_num_samples = 300,func=get_coxph_model,reps=10,time,event,predictors)
		batch_predictions_obj = predict(models,as.data.frame(te_x),type="risk")			
		preds = c(preds,batch_predictions_obj)
	}
	preds = preds[rownames(x)]
	x = cbind(x,preds)
	colnames(x)[ncol(x)] = "risk_aki_intensive"
	kidney_risks = x[,"risk_aki_intensive"]	
	cv_risks = x[,"risk_cv_standard"]
	x_events = x[which(x[,"EVENT_PRIMARY"]==1 | x[,"AKI_ERS_EVNT"]==1),]	
	event_type = rep(1,nrow(x_events))
	event_type[which(x_events[,"EVENT_PRIMARY"]==1)] = 0
	x_events = cbind(x_events,event_type)
	x_assign = cbind(cv_risks,kidney_risks)
	rownames(x_assign) = rownames(x)
	batches = sample(c(1:10),nrow(x),replace=T)
	preds = c()
	for(batch in unique(batches))
	{
		test_samples = which(batches==batch)
		tr_x = x_assign[-test_samples,]
		tr_x = tr_x[which(rownames(tr_x) %in% rownames(x_events)),]
		tr_y = x_events[which(rownames(x_events) %in% rownames(tr_x)),"event_type"]
		te_x = x_assign[test_samples,]
		models = simple_sampling_based_learner(tr_x,tr_y,reps=200)
		batch_predictions_obj = predict(models,as.data.frame(te_x),type="response")
		preds = c(preds,batch_predictions_obj)
	}
	preds = preds[rownames(x)]
	preds = 1-preds
	new_assignment = rep(0,nrow(x))	
	new_assignment[preds > theta] = 1
	return(new_assignment)
}
	
