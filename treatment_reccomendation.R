# This function accepts a SPRINT information matrix (x) and assigns patients to intensive/standard treatments 
# According to thier CV/kidney risks. theta is the treatment cutoff (final predicted value > theta: intensive treatment) 
treatment_reccomendation<-function(x,cv_event,cv_time,adverse_event,adverse_time)
{
	cv_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","FORMER_SMOKER","CURR_SMOKER","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
					,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER","TRR",
					"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO")
	aki_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","FORMER_SMOKER","CURR_SMOKER","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
					,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER","TRR",
					"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO")

	x = cbind(x,rep(0,nrow(x)))
	x_1 = x
	x_0 = x
	x_1[,"INTENSIVE"] = 1
	x_0[,"INTENSIVE"] = 0
	y= cbind(as.integer(x[,"T_PRIMARY"]+1),x[,"EVENT_PRIMARY"])
	colnames(y) = c("time","status")
	rownames(y) = rownames(x)
	if(class(y)=="matrix") events = y[,"status"] else  events = y 
	weights = rep((table(events)[1]/length(events)),length(events))
	weights[events==1] = table(events)[2]/length(events)
	names(weights) = names(events)
	cv = cv.glmnet(as.matrix(x[,cv_static_predictors]),y,family="cox",alpha=1,weights=weights)
	preds = c()
	batches = sample(c(1:10),nrow(x),replace=T)
	for(batch in unique(batches))
	{
		test_samples = which(batches==batch)
		tr_x = x[-test_samples,];tr_y=y[-test_samples,]
		te_x = x_0[test_samples,];te_y = y[test_samples,]
		model = run_coxnet(as.matrix(tr_x[,cv_static_predictors]),tr_y,cv$lambda,weights)
		batch_predictions_obj = predict(model,as.matrix(te_x[,cv_static_predictors]),s=cv$lambda.1se)[,1]
		names(batch_predictions_obj) = rownames(te_x)
		preds = c(preds,batch_predictions_obj)
	}
	preds = preds[rownames(x)]
	risk_cv_standard = preds
	y = cbind(as.integer(x[,"AKI_ERS_DAYS"]+1),x[,"AKI_ERS_EVNT"])
	colnames(y) = c("time","status")
	rownames(y) = rownames(x)
	if(class(y)=="matrix") events = y[,"status"] else  events = y 
	weights = rep((table(events)[1]/length(events)),length(events))
	weights[events==1] = table(events)[2]/length(events)
	names(weights) = names(events)
	cv = cv.glmnet(as.matrix(x[,aki_static_predictors]),y,family="cox",alpha=1,weights=weights)
	preds = c()
	batches = sample(c(1:10),nrow(x),replace=T)
	for(batch in unique(batches))
	{
		test_samples = which(batches==batch)
		tr_x = x[-test_samples,];tr_y=y[-test_samples,]
		te_x = x_1[test_samples,];te_y = y[test_samples,]
		model = run_coxnet(as.matrix(tr_x[,aki_static_predictors]),tr_y,cv$lambda,weights)
		batch_predictions_obj = predict(model,as.matrix(te_x[,aki_static_predictors]),s=cv$lambda.1se)[,1]
		names(batch_predictions_obj) = rownames(te_x)
		preds = c(preds,batch_predictions_obj)
	}
	preds = preds[rownames(x)]
	risk_aki_intensive = preds
	kidney_risks = risk_aki_intensive	
	cv_risks = risk_cv_standard
	x_events = x[which(x[,cv_event]==1 | x[,adverse_event]==1),]	
	event_type = rep(1,nrow(x_events))
	event_type[which(x_events[,ifelse(table(x[,cv_event])[2] > table(x[,adverse_event])[2],cv_event,adverse_event)]==1)] = 0
	x_events = cbind(x_events,event_type)
	x_assign = cbind(cv_risks,kidney_risks)
	rownames(x_assign) = rownames(x)
	batches = sample(c(1:10),nrow(x),replace=T)
	new_assignment = rep(0,nrow(x))
	names(new_assignment) = rownames(x)
	for(batch in unique(batches))
	{
		test_samples = which(batches==batch)
		tr_x = x_assign[-test_samples,]
		tr_x = tr_x[which(rownames(tr_x) %in% rownames(x_events)),]
		tr_y = x_events[which(rownames(x_events) %in% rownames(tr_x)),"event_type"]
		te_x = x_assign[test_samples,]
		models = simple_sampling_based_learner(tr_x,tr_y,func=run_glm)
		preds = 1-predict(models,as.data.frame(tr_x),type="response")
		chosen_theta = 0
		highest_HR = 0
		aki_coef = c()
		cv_coef = c()
		for(j in seq(0.5,0.6,0.01))
		{
			assignment = rep(0,nrow(tr_x))
			names(assignment) = rownames(tr_x)
			assignment[preds > j] = 1
			if(length(table(assignment)) == 1) { next }
			if(table(assignment)[1] < length(assignment)/100 | table(assignment)[2] < length(assignment)/100)
			{
				aki_coef = c(aki_coef,0)
				cv_coef = c(cv_coef,0)
				next 
			}
			cox = coxph(Surv(time=x[rownames(tr_x),"T_PRIMARY"],event=x[rownames(tr_x),"EVENT_PRIMARY"]) ~ assignment,
					data=cbind(x[rownames(tr_x),],assignment))
			cox_aki = coxph(Surv(time=x[rownames(tr_x),"AKI_ERS_DAYS"],event=x[rownames(tr_x),"AKI_ERS_EVNT"]) ~ assignment,
					data=cbind(x[rownames(tr_x),],assignment))
			cv_coef = c(cv_coef,exp(cox$"coefficients"))
			aki_coef = c(aki_coef,exp(cox_aki$"coefficients"))
			#print(j)
			#print(exp(cox$"coefficients"))
			#print(exp(cox_aki$"coefficients"))
			if(is.na(exp(cox_aki$"coefficients"))){ next }
			if(exp(cox$"coefficients") > highest_HR & exp(cox_aki$"coefficients") <=1)
			#if(exp(cox$"coefficients") > highest_HR)
			{
				highest_HR = exp(cox$"coefficients")
				chosen_theta = j
			}
		}
		#df_means = cbind(aki_coef,cv_coef)
		#df_means = df_means[1:min(min(which(is.na(aki_coef))),min(which(is.na(cv_coef))))-1,]
		#rownames(df_means) = alphas[1:min(min(which(is.na(aki_coef))),min(which(is.na(cv_coef))))-1]
		
		#df.means <- melt(df_means, id.vars = "x")
		#colnames(df.means) = c("x","event","value")
		#df.means[,"x"] = c(alphas,alphas)
		#plot = ggplot(data=df.means,aes(x=x,y=value,colour=event)) + scale_colour_manual(values=c("red","blue"))+ geom_point(size=1) + geom_line() + labs(x="",y="HR (intensive vs. standard)")  + theme(axis.title.y = element_text(size=26))
		print(chosen_theta)
		preds = 1-predict(models,as.data.frame(te_x),type="response")
		new_assignment[names(preds)[preds > chosen_theta]] = 1
	}
	return(new_assignment)
}
	
