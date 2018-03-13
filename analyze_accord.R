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
curr_smoker = rep(0,nrow(x))
curr_smoker[which(x[,"SMOKE_3CAT"]==3)] = 1
former_smoker = rep(0,nrow(x))
former_smoker[which(x[,"SMOKE_3CAT"]==2)] = 1
x = cbind(x,curr_smoker,former_smoker)
pp = x_accord[,"sbp"] - x_accord[,"dbp"]
sub_ckd = rep(0,nrow(x_accord))
sub_ckd[which(x_accord[,"gfr"] < 60)] = 1
hdl_chol_ratio = x_accord[,"hdl"]/x_accord[,"chol"]
x_accord = cbind(x_accord,sub_ckd,pp,hdl_chol_ratio)

curr_smoker = rep(0,nrow(x_6m_prior_sprint))
curr_smoker[which(x_6m_prior_sprint[,"SMOKE_3CAT"]==3)] = 1
former_smoker = rep(0,nrow(x_6m_prior_sprint))
former_smoker[which(x_6m_prior_sprint[,"SMOKE_3CAT"]==2)] = 1
x_6m_prior_sprint = cbind(x_6m_prior_sprint,curr_smoker,former_smoker)
pp = x_6m_prior_accord[,"sbp"] - x_6m_prior_accord[,"dbp"]
sub_ckd = rep(0,nrow(x_6m_prior_accord))
sub_ckd[which(x_6m_prior_accord[,"gfr"] < 60)] = 1
hdl_chol_ratio = x_6m_prior_accord[,"hdl"]/x_6m_prior_accord[,"chol"]
x_6m_prior_accord = cbind(x_6m_prior_accord,sub_ckd,pp,hdl_chol_ratio)

curr_smoker = rep(0,nrow(x_12m_prior_sprint))
curr_smoker[which(x_12m_prior_sprint[,"SMOKE_3CAT"]==3)] = 1
former_smoker = rep(0,nrow(x_12m_prior_sprint))
former_smoker[which(x_12m_prior_sprint[,"SMOKE_3CAT"]==2)] = 1
x_12m_prior_sprint = cbind(x_12m_prior_sprint,curr_smoker,former_smoker)
pp = x_12m_prior_accord[,"sbp"] - x_12m_prior_accord[,"dbp"]
sub_ckd = rep(0,nrow(x_12m_prior_accord))
sub_ckd[which(x_12m_prior_accord[,"gfr"] < 60)] = 1
hdl_chol_ratio = x_12m_prior_accord[,"hdl"]/x_12m_prior_accord[,"chol"]
x_12m_prior_accord = cbind(x_12m_prior_accord,sub_ckd,pp,hdl_chol_ratio)

#predict ACCORD patients
cv_static_predictors = c("INTENSIVE","N_AGENTS","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"CHR","HDL","UMALCR","BMI","STATIN","PP","curr_smoker","former_smoker")
cv_dynamic_predictors = c("INTENSIVE","N_AGENTS","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER","CHR","HDL",
				"UMALCR","BMI","STATIN","PP","curr_smoker","former_smoker","pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")
cv_static_predictors_accord = c("arm","n_agents","gfr","screat","sub_ckd","ETHNIC_BLACK","baseline_age","female","cvd_hx_baseline","ETHNIC_WHITE"
					,"ETHNIC_HISPANIC","ETHNIC_OTHER","chol","hdl","uacr","bmi","statin","pp","curr_smoker","former_smoker")
cv_dynamic_predictors_accord = c("arm","n_agents","gfr","screat","sub_ckd","ETHNIC_BLACK","baseline_age","female","cvd_hx_baseline","ETHNIC_WHITE"
					,"ETHNIC_HISPANIC","ETHNIC_OTHER","chol","hdl","uacr","bmi","statin","pp","curr_smoker","former_smoker",
					"pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")

#calculate static CV model performance
tr_x = x[,cv_static_predictors];tr_y=x[,"EVENT_PRIMARY"]
te_x = x_accord[,cv_static_predictors_accord];te_y = x_accord[,"censor_po"]
colnames(te_x) = colnames(tr_x)
accord_static_rocs = c()
for(i in 1:50)
{
models = simple_sampling_based_learner(tr_x,tr_y)
cv_results = predict(models,as.data.frame(te_x))
fg = cv_results[te_y == 1]; bg = cv_results[te_y == 0]
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
accord_static_rocs = c(accord_static_rocs,roc$auc)
}

#calculate dynamic CV model performance
tr_x = x_12m_prior_sprint[,cv_dynamic_predictors];tr_y=x_12m_prior_sprint[,"EVENT_PRIMARY"]
te_x = x_12m_prior_accord[,cv_dynamic_predictors_accord];te_y = x_12m_prior_accord[,"censor_po"]
colnames(te_x) = colnames(tr_x)
accord_dynamic_12_t_rocs = c()
for(i in 1:50)
{
models = simple_sampling_based_learner(tr_x,tr_y)
cv_results = predict(models,as.data.frame(te_x))
fg = cv_results[te_y == 1]; bg = cv_results[te_y == 0]
roc<-roc.curve(scores.class0 = fg, scores.class1 = bg,curve = TRUE)
accord_dynamic_12_t_rocs = c(accord_dynamic_12_t_rocs,roc$auc)
}

tr_x = x_6m_prior_sprint[,cv_dynamic_predictors];tr_y=x_6m_prior_sprint[,"EVENT_PRIMARY"]
te_x = x_6m_prior_accord[,cv_dynamic_predictors_accord];te_y = x_6m_prior_accord[,"censor_po"]
colnames(te_x) = colnames(tr_x)
accord_dynamic_6_t_rocs = c()
for(i in 1:50)
{
models = simple_sampling_based_learner(tr_x,tr_y)
cv_results = predict(models,as.data.frame(te_x))
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