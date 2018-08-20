needed_libs = c('PRROC','impute','survival','sas7bdat','ggplot2','glmnet','rms')
for(j in 1:length(needed_libs)){
	try({library(needed_libs[j],character.only=T)})
}
needed_scripts = c("./extract_dynamics_ACCORD.R","./extract_dynamics_SPRINT.R","./generate_sprint_raw_data.R","./generate_accord_raw_data.R",
			"./fix_framingham_score.R","./predictors.R","./run_predictor.R","./treatment_reccomendation.R")
sapply(needed_scripts,function(x) source(x))
#preprocess
options(warn=-1)
raw_data = generate_raw_data()
FORMER_SMOKER = (raw_data[,"SMOKE_3CAT"] == 2)*1
CURR_SMOKER = (raw_data[,"SMOKE_3CAT"] == 3)*1
feature_set = c("INTENSIVE","RISK10YRS","INCLUSIONFRS","SBP","DBP","N_AGENTS","NOAGENTS","ASPIRIN","EGFR","SCREAT","SUB_CKD","RACE_BLACK",
			"AGE","FEMALE","SUB_CVD","SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
			"CHR","GLUR","HDL","TRR","UMALCR","BMI","STATIN")
events_set = c("EVENT_PRIMARY","T_PRIMARY","EVENT_30PERCENTREDUCTION_EGFR","T_30PERCENTREDUCTION_EGFR",
			"AKI_ERS_EVNT","AKI_ERS_DAYS")
x = raw_data[,c(feature_set,events_set)]
x = cbind(x,FORMER_SMOKER,CURR_SMOKER)
x[,feature_set] = impute.knn(as.matrix(x[,feature_set]))$data
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

#The whole analysis might take up to ~90 minuts
#run univariate CV model
cv_uni_predictors = c("INTENSIVE")
y_cv_uni = x[,"EVENT_PRIMARY"]
res = run_predictor(y = y_cv_uni,command = run_survival_cv_by_batch,args = list(x=x,y=y_cv_uni,b=sample(c(1:10),nrow(x),replace=T),
			list_of_predictors = cv_uni_predictors,prediction_algo="coxph",time = "T_PRIMARY",event = "EVENT_PRIMARY"))
cv_uni_rocs = res[1,]


#run static CV model
cv_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","FORMER_SMOKER","CURR_SMOKER","EGFR","SCREAT","EGFR_SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","HDL","CHR")
y_cv_static = cbind(as.integer(x[,"T_PRIMARY"]+1),x[,"EVENT_PRIMARY"])
colnames(y_cv_static) = c("time","status")
rownames(y_cv_static) = rownames(x)
res = run_predictor(y = y_cv_static,command = run_cv_by_batch,args = list(x=as.matrix(x[,cv_static_predictors]),y=y_cv_static,b=sample(c(1:10),nrow(x),replace=T),
			classification_function=run_coxnet,method="cox"))
cv_static_rocs = res[1,]

#run univariate AKI model
aki_uni_predictors = c("INTENSIVE")
y_aki_uni= x[,"AKI_ERS_EVNT"]
res = run_predictor(y = y_aki_uni,command = run_survival_cv_by_batch,args = list(x=x,y=y_aki_uni,b=sample(c(1:10),nrow(x),replace=T),
			list_of_predictors = aki_uni_predictors,prediction_algo="coxph",time = "AKI_ERS_DAYS",event = "AKI_ERS_EVNT"))
aki_uni_rocs = res[1,]

#run static AKI model
aki_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","FORMER_SMOKER","CURR_SMOKER","EGFR","SCREAT","EGFR_SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","HDL","CHR")
y_aki_static = cbind(as.integer(x[,"AKI_ERS_DAYS"]+1),x[,"AKI_ERS_EVNT"])
colnames(y_aki_static) = c("time","status")
rownames(y_aki_static) = rownames(x)
res = run_predictor(y = y_aki_static,command = run_cv_by_batch,args = list(x=as.matrix(x[,aki_static_predictors]),y=y_aki_static,b=sample(c(1:10),nrow(x),replace=T),
			classification_function=run_coxnet,method="cox"))
aki_static_rocs = res[1,]

#run static CKD model
ckd_uni_predictors = c("INTENSIVE")
x_ckd_uni = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),]
y_ckd_uni = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),"EVENT_30PERCENTREDUCTION_EGFR"]
res = run_predictor(y = y_ckd_uni,command = run_survival_cv_by_batch,args = list(x=x_ckd_uni,y=y_ckd_uni,b=sample(c(1:10),nrow(x_ckd_uni),replace=T),
			list_of_predictors = ckd_uni_predictors,prediction_algo="coxph",time = "T_30PERCENTREDUCTION_EGFR",event = "EVENT_30PERCENTREDUCTION_EGFR"))
ckd_uni_rocs = res[1,]


#run static CKD model
ckd_static_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","FORMER_SMOKER","CURR_SMOKER","EGFR","SCREAT","EGFR_SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","HDL","CHR")
x_ckd_static = x[which(x[,"SUB_CKD"]==0 & !is.na(x[,"EVENT_30PERCENTREDUCTION_EGFR"])),]
y_ckd_static = cbind(as.integer(x_ckd_static[,"T_30PERCENTREDUCTION_EGFR"]+1),x_ckd_static[,"EVENT_30PERCENTREDUCTION_EGFR"])
colnames(y_ckd_static) = c("time","status")
rownames(y_ckd_static) = rownames(x_ckd_static)
res = run_predictor(y = y_ckd_static,command = run_cv_by_batch,args = list(x=as.matrix(x_ckd_static[,ckd_static_predictors]),y=y_ckd_static,b=sample(c(1:10),nrow(x_ckd_static),replace=T),
			classification_function=run_coxnet,method="cox"))
ckd_static_rocs = res[1,]

#run dynamic CV model
cv_dynamic_predictors = c("INTENSIVE","RISK10YRS","N_AGENTS","FORMER_SMOKER","CURR_SMOKER","EGFR","SCREAT","EGFR_SCREAT","SUB_CKD","RACE_BLACK","AGE","FEMALE","SUB_CVD"
				,"SUB_CLINICALCVD","SUB_SUBCLINICALCVD","SUB_SENIOR","ETHNIC_WHITE_RACE","ETHNIC_HISPANIC","ETHNIC_OTHER",
				"UMALCR","BMI","STATIN","PP","HDL_CHR_RATIO","HDL","CHR","pp_m","pp_mean","pp_R2","pp_f_pvalue","pp_max_min")
x_cv_dynamic = x_6m_prior[,cv_dynamic_predictors] 
y_cv_dynamic = x_6m_prior[,"EVENT_PRIMARY"]
names(y_cv_dynamic) = rownames(x_cv_dynamic)
res = run_predictor(y = y_cv_dynamic,command = run_cv_by_batch,args = list(x = as.matrix(x_6m_prior[,cv_dynamic_predictors]),y = y_cv_dynamic,b = sample(c(1:10),nrow(x_6m_prior),replace=T),
			classification_function=simple_sampling_based_learner,method="binomial"))
cv_dynamic_t_6_rocs = res[1,]

x_cv_dynamic = x_12m_prior[,cv_dynamic_predictors] 
y_cv_dynamic = x_12m_prior[,"EVENT_PRIMARY"]
names(y_cv_dynamic) = rownames(x_cv_dynamic)
res = run_predictor(y = y_cv_dynamic,command = run_cv_by_batch,args = list(x = as.matrix(x_12m_prior[,cv_dynamic_predictors]),y = y_cv_dynamic,b = sample(c(1:10),nrow(x_12m_prior),replace=T),
			classification_function=simple_sampling_based_learner,method="binomial",prediction_args=list(type="response")))
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
performance = c(mean(cv_static_rocs),mean(cv_uni_rocs),mean(ckd_static_rocs),mean(ckd_uni_rocs),mean(aki_static_rocs),mean(aki_uni_rocs))
sd = c(sd(cv_static_rocs),sd(cv_uni_rocs),sd(ckd_static_rocs),sd(ckd_uni_rocs),sd(aki_static_rocs),sd(aki_uni_rocs))
performance_measure = c("AUC ROC","AUC ROC","AUC ROC","AUC ROC","AUC ROC","AUC ROC") 
names_rocs = c("Primary Outcome","Primary Outcome","Chronic Kidney Disease","Chronic Kidney Disease","Acute Kidney Failure","Acute Kidney Failure")
model = c("Multivariate","Univariate","Multivariate","Univariate","Multivariate","Univariate")
df = data.frame(names_rocs,model,performance_measure,performance,sd)
df[,1] = factor(df$"names_rocs",levels = c("Primary Outcome","Chronic Kidney Disease","Acute Kidney Failure"))
df[,2] = factor(df$"model",levels = c("Multivariate","Univariate"))
ggplot(df,aes(x = model,y = performance, fill = df$performance_measure)) + geom_col(position = 'dodge') +
coord_cartesian(ylim=c(0.5,0.8)) + ylab("Performance Score") + xlab("Outcome Predicted") + 
facet_wrap(~names_rocs, strip.position = "bottom", scales = "free_x") + scale_fill_discrete(name = "Performance measurement") +
geom_errorbar(aes(ymin=performance-sd, ymax=performance+sd),width=.2,size = 1,position=position_dodge(.9))+
ggtitle(" Comparison of predictive models in CV and kidney event prediction \n Static Multivariate vs. Univariate") + theme(panel.margin = unit(0, "lines"), strip.background = element_blank(),text = element_text(size=18))

# run treatment recommendation simulation
# this might take a while (~30min)
num_of_iterations = 101
assignment_matrix = matrix(ncol = nrow(x))
colnames(assignment_matrix) = rownames(x)
for(j in 1:num_of_iterations)
{
	new_assignment = treatment_reccomendation(x,"EVENT_PRIMARY","T_PRIMARY","AKI_ERS_EVNT","AKI_ERS_DAYS")
	assignment_matrix = rbind(assignment_matrix,new_assignment)
}
assignment_matrix = assignment_matrix[-1,]
new_assignment = apply(assignment_matrix,2,function(x) length(which(x==1)) > floor(num_of_iterations/2))*1
if("new_assignment" %in% colnames(x))
{
	x = x[,-which(colnames(x) == "new_assignment")]
}
x = cbind(x,new_assignment)
time = "AKI_ERS_DAYS"
event = "AKI_ERS_EVNT"
batch_surv = Surv(x[,time],x[,event])
batch_model = coxph(as.formula("batch_surv ~ new_assignment"),data=as.data.frame(x))
aki_pvalue = summary(batch_model)$"sctest"[3]
hr = c(exp(batch_model$coefficients))
lower_ci = c(summary(batch_model)$"conf.int"[3])
upper_ci = c(summary(batch_model)$"conf.int"[4])
time = "T_PRIMARY"
event = "EVENT_PRIMARY"
batch_surv = Surv(x[,time],x[,event])
batch_model = coxph(as.formula("batch_surv ~ new_assignment"),data=as.data.frame(x))
cv_pvalue = summary(batch_model)$"sctest"[3]
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
  	scale_y_continuous(limits = c(0, 10),breaks = seq(0, 10, by = 1)) +
  	coord_flip() +
  	theme_minimal() + ggtitle("Hazard ratios for CV event and Renal failure between the recommended intensive and standard groups")+
	theme(legend.position="none",axis.title.y =element_text(size=16),axis.title.x =element_text(size=14),axis.text.x =element_text(size=14),axis.text.y =element_text(size=14))
print(paste("P-value for greater chance of CV event in the reccomended intensive group is:",formatC(signif(cv_pvalue,digits=3))))
print(paste("P-value for lower chance of AKI event in the reccomended intensive group is:",formatC(signif(aki_pvalue,digits=3))))