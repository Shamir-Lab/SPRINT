#The whole process might take ~40 minuts to run
#preprocess
raw_accord_data = generate_accord_raw_data()
x_accord = raw_accord_data
x_accord = impute.knn(as.matrix(x_accord))$data
bp = read.sas7bdat("./bloodpressure.sas7bdat")
x_6m_prior_accord = preprocess_accord(x_accord,6,bp,calculate_cox_accord)
x_12m_prior_accord = preprocess_accord(x_accord,12,bp,calculate_cox_accord)
x_6m_prior_sprint = x_6m_prior
x_12m_prior_sprint = x_12m_prior

#fix SPRINT-ACCORD data discrapencies
pp = x_accord[,"sbp"] - x_accord[,"dbp"]
sub_ckd = rep(0,nrow(x_accord))
sub_ckd[which(x_accord[,"gfr"] < 60)] = 1
hdl_chol_ratio = x_accord[,"hdl"]/x_accord[,"chol"]
x_accord = cbind(x_accord,sub_ckd,pp,hdl_chol_ratio)

pp = x_6m_prior_accord[,"sbp"] - x_6m_prior_accord[,"dbp"]
sub_ckd = rep(0,nrow(x_6m_prior_accord))
sub_ckd[which(x_6m_prior_accord[,"gfr"] < 60)] = 1
hdl_chol_ratio = x_6m_prior_accord[,"hdl"]/x_6m_prior_accord[,"chol"]
x_6m_prior_accord = cbind(x_6m_prior_accord,sub_ckd,pp,hdl_chol_ratio)

pp = x_12m_prior_accord[,"sbp"] - x_12m_prior_accord[,"dbp"]
sub_ckd = rep(0,nrow(x_12m_prior_accord))
sub_ckd[which(x_12m_prior_accord[,"gfr"] < 60)] = 1
hdl_chol_ratio = x_12m_prior_accord[,"hdl"]/x_12m_prior_accord[,"chol"]
x_12m_prior_accord = cbind(x_12m_prior_accord,sub_ckd,pp,hdl_chol_ratio)

#predict ACCORD patients
cv_static_predictors = c("INTENSIVE","N_AGENTS","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"CHR","HDL","UMALCR","BMI","STATIN","PP","FORMER_SMOKER","CURR_SMOKER")
cv_dynamic_predictors = c("INTENSIVE","N_AGENTS","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER","CHR","HDL","TRR",
				"UMALCR","BMI","STATIN","PP","FORMER_SMOKER","CURR_SMOKER","pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")
cv_static_predictors_accord = c("arm","n_agents","gfr","screat","sub_ckd","ETHNIC_BLACK","baseline_age","female","cvd_hx_baseline","ETHNIC_WHITE"
					,"ETHNIC_HISPANIC","ETHNIC_OTHER","chol","hdl","uacr","bmi","statin","pp","curr_smoker","former_smoker")
cv_dynamic_predictors_accord = c("arm","n_agents","gfr","screat","sub_ckd","ETHNIC_BLACK","baseline_age","female","cvd_hx_baseline","ETHNIC_WHITE"
					,"ETHNIC_HISPANIC","ETHNIC_OTHER","chol","hdl","trig","uacr","bmi","statin","pp","curr_smoker","former_smoker",
					"pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")

#calculate static CV model performance

tr_x = x[,cv_static_predictors];tr_y=x[,"EVENT_PRIMARY"]
te_x = x_accord[,cv_static_predictors_accord];te_y = x_accord[,"censor_po"]
colnames(te_x) = colnames(tr_x)
names(tr_y) = rownames(tr_x)
weights = rep((table(tr_y)[1]/length(tr_y)),length(tr_y))
weights[tr_y==1] = table(tr_y)[2]/length(tr_y)
names(weights) = names(tr_y)
cv = cv.glmnet(as.matrix(tr_x),tr_y,family="binomial",alpha=1,weights=weights)
accord_static_rocs = c()
for(i in 1:50)
{
models = simple_sampling_based_learner(x=as.matrix(tr_x),y=tr_y,lambda = cv$lambda,weights=weights)
pred = 0
num_of_empty_models = 0
for(j in 1:length(models))
{
	if(length(models[[j]]$lambda) == 1)
	{
		num_of_empty_models = num_of_empty_models + 1
		next
	}
	curr_pred = te_x %*% coef(models[[j]],s=cv$lambda.1se)[-1]
	adj_model = lrm.fit(y=te_y,offset = curr_pred)
	pred = pred + exp(curr_pred+coef(adj_model)[1])/(1+exp(curr_pred+coef(adj_model)[1]))
}
pred = pred/(length(models) - num_of_empty_models)
cv_results = pred
fg = cv_results[te_y == 1]; bg = cv_results[te_y == 0]
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
accord_static_rocs = c(accord_static_rocs,roc$auc)
}

#calculate dynamic CV model performance
tr_x = x_12m_prior_sprint[,cv_dynamic_predictors];tr_y=x_12m_prior_sprint[,"EVENT_PRIMARY"]
te_x = x_12m_prior_accord[,cv_dynamic_predictors_accord];te_y = x_12m_prior_accord[,"censor_po"]
colnames(te_x) = colnames(tr_x)
names(tr_y) = rownames(tr_x)
weights = rep((table(tr_y)[1]/length(tr_y)),length(tr_y))
weights[tr_y==1] = table(tr_y)[2]/length(tr_y)
names(weights) = rownames(tr_x)
cv = cv.glmnet(as.matrix(tr_x),tr_y,family="binomial",alpha=1,weights=weights)
accord_dynamic_12_t_rocs = c()
for(i in 1:50)
{
models = simple_sampling_based_learner(x=as.matrix(tr_x),y=tr_y,lambda = cv$lambda,weights=weights)
pred = 0
num_of_empty_models = 0
for(j in 1:length(models))
{
	if(length(models[[j]]$lambda) == 1)
	{
		num_of_empty_models = num_of_empty_models + 1
		next
	}
	curr_pred = te_x %*% coef(models[[j]],s=cv$lambda.1se)[-1]
	adj_model = lrm.fit(y=te_y,offset = curr_pred)
	pred = pred + exp(curr_pred+coef(adj_model)[1])/(1+exp(curr_pred+coef(adj_model)[1]))
}
pred = pred/(length(models) - num_of_empty_models)
cv_results = pred
fg = cv_results[te_y == 1]; bg = cv_results[te_y == 0]
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
accord_dynamic_12_t_rocs = c(accord_dynamic_12_t_rocs,roc$auc)
}

tr_x = x_6m_prior_sprint[,cv_dynamic_predictors];tr_y=x_6m_prior_sprint[,"EVENT_PRIMARY"]
te_x = x_6m_prior_accord[,cv_dynamic_predictors_accord];te_y = x_6m_prior_accord[,"censor_po"]
colnames(te_x) = colnames(tr_x)
names(tr_y) = rownames(tr_x)
weights = rep((table(tr_y)[1]/length(tr_y)),length(tr_y))
weights[tr_y==1] = table(tr_y)[2]/length(tr_y)
names(weights) = names(tr_y)
cv = cv.glmnet(as.matrix(tr_x),tr_y,family="binomial",alpha=1,weights=weights)
accord_dynamic_6_t_rocs = c()
for(i in 1:50)
{
models = simple_sampling_based_learner(x=as.matrix(tr_x),y=tr_y,lambda = cv$lambda,weights=weights)
pred = 0
num_of_empty_models = 0
for(j in 1:length(models))
{
	if(length(models[[j]]$lambda) == 1)
	{
		num_of_empty_models = num_of_empty_models + 1
		next
	}
	curr_pred = te_x %*% coef(models[[j]],s=cv$lambda.1se)[-1]
	adj_model = lrm.fit(y=te_y,offset = curr_pred)
	pred = pred + exp(curr_pred+coef(adj_model)[1])/(1+exp(curr_pred+coef(adj_model)[1]))
}
pred = pred/(length(models) - num_of_empty_models)
cv_results = pred[,1]
names(cv_results) = rownames(pred)
fg = cv_results[te_y == 1]; bg = cv_results[te_y == 0]
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
accord_dynamic_6_t_rocs = c(accord_dynamic_6_t_rocs,roc$auc)
}


#plot results
cv_mean_rocs = c(mean(accord_dynamic_6_t_rocs),mean(accord_dynamic_12_t_rocs),mean(accord_static_rocs))
cv_sd_rocs = c(sd(accord_dynamic_6_t_rocs),sd(accord_dynamic_12_t_rocs),sd(accord_static_rocs))


model = c("Dynamic \n Multivariate \n (t=6)","Dynamic \n Multivariate \n (t=12)","Static \n Multivariate")
df = data.frame(model,cv_mean_rocs,cv_sd_rocs)
df[,1] = factor(df$"model",levels = df$"model")
ggplot(df,aes(model,cv_mean_rocs,fill = model)) + geom_col(position = 'dodge') +
coord_cartesian(ylim=c(0.5,0.8)) + ylab("AUC ROC for Primary Outcome") + xlab("Prediction Model") + 
ggtitle("Comparison of Dynamic vs. Static multivariate predictive models in CV event prediction for the ACCORD cohort \n(training set: SPRINT) ") +
theme(legend.position="none",axis.title.y =element_text(size=16),
axis.title.x = element_text(size=20),axis.text.x = element_text(size=14), text = element_text(size=22)) +
geom_errorbar(aes(ymin=cv_mean_rocs-cv_sd_rocs, ymax=cv_mean_rocs+cv_sd_rocs),width=.1,size = 1,position=position_dodge(.9))