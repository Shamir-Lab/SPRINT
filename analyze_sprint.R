needed_libs = c('PRROC','impute','survival','sas7bdat','ggplot2')
for(j in 1:length(needed_libs)){
	try({library(needed_libs[j],character.only=T)})
}
needed_scripts = c("./extract_dynamics_ACCORD.R","./extract_dynamics_SPRINT.R","./generate_sprint_raw_data.R","./generate_accord_raw_data.R",
			"./fix_framingham_score.R","./predictors.R","./run_predictor.R","./treatment_reccomendation.R")
sapply(needed_scripts,function(x) source(x))
#preprocess
options(warn=-1)
raw_data = generate_raw_data()
feature_set = c("INTENSIVE","RISK10YRS","INCLUSIONFRS","SBP","DBP","N_AGENTS","NOAGENTS","SMOKE_3CAT","ASPIRIN","EGFR","SCREAT","SUB_CKD","RACE_BLACK",
			"AGE","FEMALE","SUB_CVD","SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
			"CHR","GLUR","HDL","TRR","UMALCR","BMI","STATIN","EVENT_PRIMARY","T_PRIMARY","EVENT_30PERCENTREDUCTION_EGFR","T_30PERCENTREDUCTION_EGFR",
			"AKI_ERS_EVNT","AKI_ERS_DAYS")
x = raw_data[,feature_set]
x = impute.knn(as.matrix(x))$data
x[,"ASPIRIN"] = round(x[,"ASPIRIN"])
x[,"STATIN"] = round(x[,"STATIN"])
PP = x[,"SBP"] - x[,"DBP"]	
HDL_CHR_RATIO = x[,"HDL"]/x[,"CHR"]
EGFR_SCREAT = x[,"EGFR"]/x[,"SCREAT"]
x = cbind(x,PP,HDL_CHR_RATIO,EGFR_SCREAT)
x[,"RISK10YRS"] = fix_framingham_score(x)
x[,"INCLUSIONFRS"] = apply(x,1,function(y) y["RISK10YRS"] > 15)
x[which(x[,"INCLUSIONFRS"]==T),"INCLUSIONFRS"] = 1 
x[which(x[,"INCLUSIONFRS"]==F),"INCLUSIONFRS"] = 0 
#calculate dynamic features
x_6m_prior = preprocess_sprint(x,6,"./bp.csv",calculate_cox_sprint)
x_12m_prior = preprocess_sprint(x,12,"./bp.csv",calculate_cox_sprint)

#run univariate CV model
cv_uni_predictors = c("INTENSIVE")
y_cv_uni = x[,"EVENT_PRIMARY"]
res = run_predictor(y = y_cv_uni,command = run_survival_cv_by_batch,args = list(x,y_cv_uni,
			list_of_predictors = cv_uni_predictors,prediction_algo="coxph",time = "T_PRIMARY",event = "EVENT_PRIMARY"),
			keep_concordance=T)
cv_uni_rocs = res[1,]
cv_uni_concordance = res[2,]


#run static CV model
cv_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","SMOKE_3CAT","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO")
y_cv_static = x[,"EVENT_PRIMARY"]
res = run_predictor(y = y_cv_static,command = run_survival_cv_by_batch,args = list(x,y_cv_static,sample(c(1:10),nrow(x),replace=T),
			list_of_predictors = cv_static_predictors,prediction_algo="coxph",time = "T_PRIMARY",event = "EVENT_PRIMARY"),
			keep_concordance=T)
cv_static_rocs = res[1,]
cv_static_concordance = res[2,]

#run univariate AKI model
aki_uni_predictors = c("INTENSIVE")
y_aki_uni= x[,"AKI_ERS_EVNT"]
res = run_predictor(y = y_aki_uni,command = run_survival_cv_by_batch,args = list(x,y_aki_uni,sample(c(1:10),nrow(x),replace=T),
			list_of_predictors = aki_uni_predictors,prediction_algo="coxph",time = "AKI_ERS_DAYS",event = "AKI_ERS_EVNT"),
			keep_concordance=T)
aki_uni_rocs = res[1,]
aki_uni_concordance = res[2,]

#run static AKI model
aki_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","SMOKE_3CAT","ASPIRIN","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER","GLUR","TRR","UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","EGFR_SCREAT")
y_aki_static = x[,"AKI_ERS_EVNT"]
res = run_predictor(y = y_aki_static,command = run_survival_cv_by_batch,args = list(x,y_aki_static,sample(c(1:10),nrow(x),replace=T),
			list_of_predictors = aki_static_predictors,prediction_algo="coxph",time = "AKI_ERS_DAYS",event = "AKI_ERS_EVNT"),
			keep_concordance=T)
aki_static_rocs = res[1,]
aki_static_concordance = res[2,]

#run static CKD model
ckd_uni_predictors = c("INTENSIVE")
x_ckd_uni = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),]
y_ckd_uni = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),"EVENT_30PERCENTREDUCTION_EGFR"]
res = run_predictor(y = y_ckd_uni,command = run_survival_cv_by_batch,args = list(x_ckd_uni,y_ckd_uni,sample(c(1:10),nrow(x_ckd_uni),replace=T),
			list_of_predictors = ckd_uni_predictors,prediction_algo="coxph",time = "T_30PERCENTREDUCTION_EGFR",event = "EVENT_30PERCENTREDUCTION_EGFR"),
			keep_concordance=T)
ckd_uni_rocs = res[1,]
ckd_uni_concordance = res[2,]


#run static CKD model
ckd_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","SMOKE_3CAT","ASPIRIN","RACE_BLACK","AGE","FEMALE","SUB_CVD","ETHNIC_WHITE_RACE",
					"ETHNIC_HISPANIC","ETHNIC_OTHER","TRR","UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","EGFR_SCREAT")
x_ckd_static = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),]
y_ckd_static = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),"EVENT_30PERCENTREDUCTION_EGFR"]
res = run_predictor(y = y_ckd_static,command = run_survival_cv_by_batch,args = list(x_ckd_static,y_ckd_static,sample(c(1:10),nrow(x_ckd_static),replace=T),
			list_of_predictors = ckd_static_predictors,prediction_algo="coxph",time = "T_30PERCENTREDUCTION_EGFR",event = "EVENT_30PERCENTREDUCTION_EGFR"),
			keep_concordance=T)
ckd_static_rocs = res[1,]
ckd_static_concordance = res[2,]

#run dynamic CV model
cv_dynamic_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","SMOKE_3CAT","EGFR","SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")
x_cv_dynamic = x_6m_prior[,cv_dynamic_predictors] 
y_cv_dynamic = x_6m_prior[,"EVENT_PRIMARY"]
res = run_predictor(y = y_cv_dynamic,command = run_cv_by_batch,args = list(x_cv_dynamic,y_cv_dynamic,sample(c(1:10),nrow(x_cv_dynamic),replace=T),
			classification_function=simple_sampling_based_learner,prediction_args = list(type="response"),method = "glm"),
			keep_concordance=F)
cv_dynamic_t_6_rocs = res[1,]

x_cv_dynamic = x_12m_prior[,cv_dynamic_predictors] 
y_cv_dynamic = x_12m_prior[,"EVENT_PRIMARY"]
res = run_predictor(y = y_cv_dynamic,command = run_cv_by_batch,args = list(x_cv_dynamic,y_cv_dynamic,sample(c(1:10),nrow(x_cv_dynamic),replace=T),
			classification_function=simple_sampling_based_learner,prediction_args = list(type="response"),method = "glm"),
			keep_concordance=F)
cv_dynamic_t_12_rocs = res[1,]

cv_mean_rocs = c(mean(cv_dynamic_t_6_rocs),mean(cv_dynamic_t_12_rocs),mean(cv_static_rocs),mean(cv_uni_rocs))
cv_sd_rocs = c(sd(cv_dynamic_t_6_rocs),sd(cv_dynamic_t_12_rocs),sd(cv_static_rocs),sd(cv_uni_rocs))

#plot cv models comparison
model = c("Dynamic \n Multivariate \n (t=6)","Dynamic \n Multivariate \n (t=12)","Static \n Multivariate","Univariate")
df = data.frame(model,cv_mean_rocs,cv_sd_rocs)
df[,1] = factor(df$"model",levels = df$"model")
ggplot(df,aes(model,cv_mean_rocs,fill = model)) + geom_col(position = 'dodge')  +
     coord_cartesian(ylim=c(0.5,0.8)) + ylab("AUC ROC for Primary CV Outcome") + xlab("Prediction Model") + 
     ggtitle("Comparison of all predictive models in CV event prediction") + theme(legend.position="none",text = element_text(size=22),axis.text.x =element_text(size=20)) +
     geom_errorbar(aes(ymin=cv_mean_rocs-cv_sd_rocs, ymax=cv_mean_rocs+cv_sd_rocs),width=.1,size = 1,position=position_dodge(.9))

#plot comparison for all static models
performance = c(mean(cv_static_concordance),mean(cv_static_rocs),mean(cv_uni_concordance),mean(cv_uni_rocs),mean(ckd_static_concordance)
			,mean(ckd_static_rocs),mean(ckd_uni_concordance),mean(ckd_uni_rocs),mean(aki_static_concordance),
			mean(aki_static_rocs),mean(aki_uni_concordance),mean(aki_uni_rocs))
sd = c(sd(cv_static_concordance),sd(cv_static_rocs),sd(cv_uni_concordance),sd(cv_uni_rocs),sd(ckd_static_concordance)
			,sd(ckd_static_rocs),sd(ckd_uni_concordance),sd(ckd_uni_rocs),sd(aki_static_concordance),
			sd(aki_static_rocs),sd(aki_uni_concordance),sd(aki_uni_rocs))
performance_measure = c("concordance","AUC ROC","concordance","AUC ROC","concordance","AUC ROC","concordance","AUC ROC","concordance","AUC ROC","concordance","AUC ROC") 
names_rocs = c("Primary Outcome","Primary Outcome","Primary Outcome","Primary Outcome","Chronic Kidney Disease"
			,"Chronic Kidney Disease","Chronic Kidney Disease","Chronic Kidney Disease","Acute Kidney Failure","Acute Kidney Failure","Acute Kidney Failure","Acute Kidney Failure")
model = c("Multivariate","Multivariate","Univariate","Univariate","Multivariate","Multivariate","Univariate","Univariate",
		"Multivariate","Multivariate","Univariate","Univariate")
df = data.frame(names_rocs,model,performance_measure,performance,sd)
df[,1] = factor(df$"names_rocs",levels = c("Primary Outcome","Chronic Kidney Disease","Acute Kidney Failure"))
df[,2] = factor(df$"model",levels = c("Multivariate","Univariate"))
ggplot(df,aes(x = model,y = performance, fill = df$performance_measure)) + geom_col(position = 'dodge') +
coord_cartesian(ylim=c(0.5,0.8)) + ylab("Performance Score") + xlab("Outcome Predicted") + 
facet_wrap(~names_rocs, strip.position = "bottom", scales = "free_x") + scale_fill_discrete(name = "Performance measurement") +
geom_errorbar(aes(ymin=performance-sd, ymax=performance+sd),width=.2,size = 1,position=position_dodge(.9))+
ggtitle(" Comparison of predictive models in CV and kidney event prediction \n Static Multivariate vs. Univariate") + theme(panel.margin = unit(0, "lines"), strip.background = element_blank(),text = element_text(size=18))

# run treatment recommendation simulation
# theta is the treatment cutoff, i.e., prediction value > 0.5: intensive treatment  
# Lower theta- assuming higher sevirity for CV even. default theta = 0.5
new_assignment = treatment_reccomendation(x,0.5)
if("new_assignment" %in% colnames(x))
{
	x = x[,-which(colnames(x) == "new_assignment")]
}
x = cbind(x,new_assignment)
time = "AKI_ERS_DAYS"
event = "AKI_ERS_EVNT"
batch_surv = Surv(x[,time],x[,event])
batch_model = coxph(as.formula("batch_surv ~ new_assignment"),data=as.data.frame(x))
hr = c(exp(batch_model$coefficients))
summary(batch_model)[7][[1]][5]		
lower_ci = c(summary(batch_model)$"conf.int"[3])
upper_ci = c(summary(batch_model)$"conf.int"[4])

time = "T_PRIMARY"
event = "EVENT_PRIMARY"
batch_surv = Surv(x[,time],x[,event])
batch_model = coxph(as.formula("batch_surv ~ new_assignment"),data=as.data.frame(x))
hr = c(hr,exp(batch_model$coefficients))
lower_ci = c(lower_ci,summary(batch_model)$"conf.int"[3])
upper_ci = c(upper_ci,summary(batch_model)$"conf.int"[4])

df_data <- data.frame(Groups=c("Renal failure","CV event"),
                      HR=hr,
                      HR_lower=lower_ci,
                      HR_upper=upper_ci)
df_data[,1] = factor(df_data$"Group",levels = df_data$"Group")

ggplot(df_data, aes(x=Groups, y=HR, ymin=HR_lower, ymax=HR_upper)) + 
  	geom_linerange(size=1, colour="black",linetype=2) +
  	geom_hline(aes(x=0, yintercept=1), lty=1) + xlab("Reccomended Intensive vs. Standard") + ylab("HR (95% CI)") + 
  	geom_point(size=3, shape=16, fill="black", colour = "black", stroke = 1) +
  	geom_errorbar(aes(ymin=HR_lower, ymax=HR_upper),width=.1,size = 1,position=position_dodge(.6),linetype=1) + 
  	scale_y_continuous(limits = c(0, 4),breaks = seq(0, 4, by = 1)) +
  	coord_flip() +
  	theme_minimal() + ggtitle("Hazard ratios for CV event and Renal failure between the recommended intensive and standard groups")+
	theme(legend.position="none",axis.title.y =element_text(size=16),axis.title.x =element_text(size=14),axis.text.x =element_text(size=14),axis.text.y =element_text(size=14))

