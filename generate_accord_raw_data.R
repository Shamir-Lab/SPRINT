# Preprocess of ACCORD data. All files must be in the currect working directory
generate_accord_raw_data<-function()
{
raw_accord_data = read.sas7bdat("./accord_key.sas7bdat")
raw_accord_data = raw_accord_data[which(raw_accord_data[,"arm"] %in% c(1,2,3,4)),]
raw_accord_data[raw_accord_data[,"arm"] %in% c(2,4),"arm"] = 0
raw_accord_data[raw_accord_data[,"arm"] %in% c(1,3),"arm"] = 1
ethnic_matrix = matrix(rep(0,4*nrow(raw_accord_data)),nrow=nrow(raw_accord_data),ncol=4)
colnames(ethnic_matrix) = c("ETHNIC_WHITE","ETHNIC_BLACK","ETHNIC_HISPANIC","ETHNIC_OTHER")
ethnic_matrix[which(raw_accord_data[,"raceclass"]=="White"),1] = 1
ethnic_matrix[which(raw_accord_data[,"raceclass"]=="Black"),2] = 1
ethnic_matrix[which(raw_accord_data[,"raceclass"]=="Hispanic"),3] = 1
ethnic_matrix[which(raw_accord_data[,"raceclass"]=="Other"),4] = 1
raw_accord_data = cbind(raw_accord_data,ethnic_matrix)
lipids = read.sas7bdat("./lipids.sas7bdat")
lipids = lipids[lipids[,"Visit"]=="BLR",]
lipids = lipids[which(lipids[,1] %in% raw_accord_data[,1]),]
raw_accord_data = merge(raw_accord_data,lipids,by="MaskID",all=T)
otherlabs = read.sas7bdat("./otherlabs.sas7bdat")
otherlabs = otherlabs[otherlabs[,"Visit"]=="BLR",]
otherlabs = otherlabs[which(otherlabs[,1] %in% raw_accord_data[,1]),]
raw_accord_data = merge(raw_accord_data,otherlabs,by="MaskID",all=T)
bp = read.sas7bdat("./bloodpressure.sas7bdat")
bp = bp[bp[,"Visit"]=="BLR",]
bp = bp[which(bp[,1] %in% raw_accord_data[,1]),]
raw_accord_data = merge(raw_accord_data,bp,by="MaskID",all=T)
concomitantmeds = read.sas7bdat("./concomitantmeds.sas7bdat")
concomitantmeds = concomitantmeds[concomitantmeds[,"Visit"]=="BLR",]
concomitantmeds = concomitantmeds[which(concomitantmeds[,1] %in% raw_accord_data[,1]),]
n_agents = apply(concomitantmeds,1,function(x) sum(as.numeric(x[c(3:16)])))
concomitantmeds = cbind(concomitantmeds,n_agents)
raw_accord_data = merge(raw_accord_data,concomitantmeds[,c(1,42,33,57)],by="MaskID",all=T)
cvdoutcomes = read.sas7bdat("./cvdoutcomes.sas7bdat")
cvdoutcomes = cvdoutcomes[which(cvdoutcomes[,1] %in% raw_accord_data[,1]),]
cvdoutcomes[,2] = 1*((1-cvdoutcomes[,10]) | (1-cvdoutcomes[,13]) | (1-cvdoutcomes[,17]) | (1-cvdoutcomes[,24]) | (1-cvdoutcomes[,31]))  
cvdoutcomes[,4] = apply(cvdoutcomes,1,function(x) min(round(as.numeric(x[11])*365),round(as.numeric(x[15])*365),
				round(as.numeric(x[19])*365),round(as.numeric(x[25])*365),round(as.numeric(x[33])*365)))
raw_accord_data = merge(raw_accord_data,cvdoutcomes[,c(1,2,4)],by="MaskID",all=T)
f07_baselinehistoryphysicalexam = read.sas7bdat("./f07_baselinehistoryphysicalexam.sas7bdat")
f07_baselinehistoryphysicalexam = f07_baselinehistoryphysicalexam[which(f07_baselinehistoryphysicalexam[,1] %in% raw_accord_data[,1]),]
curr_smoker = f07_baselinehistoryphysicalexam[,"cigarett"] -1
former_smoker = f07_baselinehistoryphysicalexam[,"smokelif"] -1
bmi = f07_baselinehistoryphysicalexam[,"wt_kg"]/((f07_baselinehistoryphysicalexam[,"ht_cm"]/100)^2)
f07_baselinehistoryphysicalexam = cbind(f07_baselinehistoryphysicalexam,curr_smoker,former_smoker,bmi)
raw_accord_data = merge(raw_accord_data,f07_baselinehistoryphysicalexam[,c(1,62:64)],by="MaskID",all=T)
kidney = read.sas7bdat("./microvascularoutcomes.sas7bdat")
kidney = kidney[which(kidney[,1] %in% raw_accord_data[,1]),]
kidney[,"Neph3Days"] = round(kidney[,"Neph3Days"]*365)
kidney[,"Neph4Days"] = round(kidney[,"Neph4Days"]*365) 
raw_accord_data = merge(raw_accord_data,kidney[,c(1,14,15,17,18)],by="MaskID",all=T)
rownames(raw_accord_data) = raw_accord_data[,1]
raw_accord_data = raw_accord_data[,-c(1,5,6,8,9,14,17,18,20:24,27,28,30,33)]
return(raw_accord_data)
}