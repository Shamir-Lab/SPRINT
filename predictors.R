run_cv_by_batch<-function(x,y,b,classification_function=simple_sampling_based_learner,
	prediction_args=list(type="response"),class_of_interest = "1",method = "svm", ...)
{
	if(class(y)=="matrix") events = y[,"status"] else  events = y 
	weights = rep((table(events)[1]/length(events)),length(events))
	weights[events==1] = table(events)[2]/length(events)
	names(weights) = names(events)
	if(method != "cox" & method != "binomial")
	{ 
		y = as.factor(y)
	} else {
		cv = cv.glmnet(x,y,family=method,alpha=1,weights=weights)
	}
	preds = c()
	batches = unique(b)
	for(batch in batches){
		test_samples = which(b==batch)
		if(method != "cox")
		{
			tr_x = as.matrix(x[-test_samples,]);tr_y = y[-test_samples]
			te_x = as.matrix(x[test_samples,]);te_y = y[test_samples]
		} else { 
			tr_x = as.matrix(x[-test_samples,]);tr_y = y[-test_samples,]
			te_x = as.matrix(x[test_samples,]);te_y = y[test_samples,]
		}
		batch_model = classification_function(as.matrix(tr_x),tr_y,lambda = cv$lambda,weights=weights)
		if(method=="glm") { te_x = as.data.frame(te_x)}
		args = c(list(batch_model,te_x,cv$lambda.min),prediction_args)
		batch_prediction_obj = do.call(predict,args=args)
		dim(batch_prediction_obj)
		te_x[setdiff(rownames(te_x),rownames(batch_prediction_obj)),]
		if(class(batch_prediction_obj) == "factor")
		{
			batch_prediction_obj = getPredProbabilities(batch_prediction_obj)[,class_of_interest]
		}
		names(batch_prediction_obj) = rownames(te_x)
		preds = c(preds,batch_prediction_obj)
	}
	preds = preds[rownames(x)]
	return(preds)
}


getPredProbabilities <- function(pred_obj){
	if(class(pred_obj) == "numeric"){ return(pred_obj) }
	probs = attr(pred_obj,"prob")
	if (!is.null(probs)){return(probs)}
	probs = attr(pred_obj,"probabilities")
	if (!is.null(probs)){return(probs)}
	if (is.element("matrix",class(pred_obj))){return (pred_obj)}
	if (length(pred_obj)>=2) {probs = pred_obj[[2]]}
	if (!is.null(probs)){return(probs)}
	return (NULL)
}
run_survival_cv_by_batch<-function(x,y,b,list_of_predictors="INTENSIVE",prediction_algo="coxph",time,event, ...){
	if(is.matrix(x)){x = as.data.frame(x)}
	predictors = paste('',list_of_predictors,collapse="",sep="+")		
	predictors = substr(predictors,2,nchar(predictors))
	preds = c()
	batches = unique(b)
	averageConcordance = 0
	for(batch in batches){
		test_samples = which(b==batch)
		tr_x = x[-test_samples,];tr_y=y[-test_samples]
		te_x = x[test_samples,];te_y = y[test_samples]
		if(identical(prediction_algo,"ranger")){
			batch_predictions_obj=run_ranger(tr_x,te_x,list_of_predictors,time,event)
		}
		else{
			batch_predictions_obj=run_coxph(tr_x,te_x,predictors,time,event,prediction_algo)
			averageConcordance = averageConcordance + (survConcordance(Surv(te_x[,time],te_x[,event]) ~ batch_predictions_obj)$concordance*nrow(te_x))/nrow(x)
		}
		preds = c(preds,batch_predictions_obj)
	}
	preds = preds[rownames(x)]
	return(c(preds))
}


run_coxph<-function(tr_x,te_x,predictors,time,event,prediction_algo=coxph){
	batch_surv = Surv(tr_x[,time],tr_x[,event])
	batch_model = coxph(as.formula(paste("batch_surv ~ ",predictors)),data=tr_x)
	if(identical(prediction_algo,"survfit")){
		batch_fit = survfit(batch_model,newdata=te_x)
		batch_fit = t(batch_fit$surv)
		batch_predictions_obj = 1-batch_fit[,ncol(batch_fit)]
	}else{
		batch_fit = predict(batch_model,newdata=as.data.frame(te_x),type="risk")
		batch_predictions_obj = batch_fit
	}
	return(batch_predictions_obj)
}

simple_sampling_based_learner<-function(x,y,d=NULL,func = run_glmnet,reps=10,max_num_samples = 300,lambda=NULL,weights=NULL,...){
	if(is.null(d)){
		d = 1:nrow(x)
		names(d) = rownames(x)
	}
	positives = rownames(x)[y=="1"]
	datasets = unique(d[positives])
	negatives = rownames(x)[y=="0" & is.element(d[rownames(x)],set=datasets)]
	bgcs = setdiff(rownames(x),positives);bgcs = setdiff(bgcs,negatives)
	classifiers = list()
	samp_size = length(positives)
	for (i in 1:reps){
		newbgcs = NULL
		if(length(bgcs)>0){
			newbgcs = sample(bgcs)[1:(samp_size)]
		}
		pos_sample = positives
		neg_sample = NULL
		if(length(negatives)>0){
			neg_sample = sample(negatives)[1:min(samp_size,length(negatives))]
		}
		curr_samples = c(pos_sample,neg_sample,newbgcs)
		curr_inds = is.element(rownames(x),curr_samples)
		newx = x[curr_inds,]
		if(identical(func,run_coxnet))
		{
			newy = y[curr_inds,]
			classifiers[[i]] = func(newx,newy,lambda,weights,...)
		} else {
			newy = y[curr_inds]
			classifiers[[i]] = func(newx,newy,lambda,weights,...)
		}
	}
	obj = classifiers
	class(obj) = "simple_sampling_based_learner"
	return(obj)
}

run_coxnet<-function(x,y,lambda,weights)
{
	if(!is.null(weights))
	{
		model = glmnet(x,y,family="cox",alpha=1,lambda = lambda,weights = weights[rownames(y)])
	} else {
		model = glmnet(x,y,family="cox",alpha=1,lambda = lambda)
	}
	return(model)
}


run_glmnet <-function(x,y,lambda,weights)
{
	model = glmnet(x, y, alpha=1, family="binomial",lambda=lambda,weights = weights[names(y)])
	return(model)
}

run_glm <-function(x,y,...)
{
	model = glm(y ~.,family=binomial(link='logit'),data=as.data.frame(x))
	return(model)
}

predict.simple_sampling_based_learner<-function(obj,x,...){
	preds = c()
	counts = 0
	for (i in 1:length(obj)){
		if (length(preds)==0){
			try({
				preds = getPredProbabilities(predict(obj[[i]],x,...))
				counts = counts + 1
			})
		}
		else{
			try({
				preds = preds + getPredProbabilities(predict(obj[[i]],x,...))
				counts = counts + 1
			})
		}
	}
	preds = preds / counts
	return (preds)
}
