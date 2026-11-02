library(survival)

## 1) SETTINGS
CFG <- list(
  input_file="data/green_space_cohort.rds", output_dir="results/green_space",
  id="patient_birth_ID", time="followup_days",
  outcomes=c(ASD="ASD_event", ADHD="ADHD_event", LD="LD_event", DSLD="DSLD_event",
             DCD="DCD_event", ID="ID_event", BD="BD_event"),
  dx_age=c(ASD="ASD_diagnosis_age", ADHD="ADHD_diagnosis_age", LD="LD_diagnosis_age",
           DSLD="DSLD_diagnosis_age", DCD="DCD_diagnosis_age", ID="ID_diagnosis_age",
           BD="BD_diagnosis_age"),
  ndvi=c(pre="ndvi_preconception", pregnancy="ndvi_pregnancy", post="ndvi_postnatal"),
  tree=c(pre="tree_preconception", pregnancy="tree_pregnancy", post="tree_postnatal"),
  iqr_ndvi=.12, censor="censor_event", censor_time="censor_time_days",
  strata=c("birth_yr","FIPS"), cluster="ZIP_Index",
  cat=c("Age_cat","race_cat","season_deliv","Alcohol_abuse_yn","Smoking_yn",
        "Substance_abuse_yn","Poor_nutrition_yn"),
  num=c("temperature","relative_humidity","pop_density","medhouseholdincome",
        "pct_no_highschool","mean_bmi","nearest_hospital_km"),
  subgroups=c("race_cat","infant_sex_cat","urbanicity")
)
dir.create(CFG$output_dir,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(CFG$input_file)
d$.time <- as.numeric(d[[CFG$time]])
d$.stratum <- interaction(d[CFG$strata],drop=TRUE)
d$.zip <- d[[CFG$cluster]]
needed <- unique(c(CFG$id,CFG$time,unname(CFG$outcomes),unname(CFG$dx_age),
 unname(CFG$ndvi),CFG$strata,CFG$cluster,CFG$cat,CFG$num,CFG$censor,CFG$censor_time))
miss <- setdiff(needed,names(d))
if(length(miss)) stop("Add or map required fields in CFG: ",paste(miss,collapse=", "))
if(anyNA(d[[CFG$id]])||anyDuplicated(d[[CFG$id]])) stop("Require one row per child.")
if(any(d$.time<=0,na.rm=TRUE)) stop("Follow-up must be positive days.")

## 2) VALIDATE NDVI AND COHORT
for(v in unname(CFG$ndvi)) if(!is.numeric(d[[v]])) stop(v," must be numeric.")
write.csv(data.frame(outcome=names(CFG$outcomes),N=nrow(d),
 events=vapply(CFG$outcomes,function(v)sum(d[[v]]==1,na.rm=TRUE),numeric(1))),
 file.path(CFG$output_dir,"cohort_counts.csv"),row.names=FALSE)
write.csv(data.frame(window=names(CFG$ndvi),field=unname(CFG$ndvi),
 mean=vapply(CFG$ndvi,function(v)mean(d[[v]],na.rm=TRUE),numeric(1)),
 IQR=vapply(CFG$ndvi,function(v)IQR(d[[v]],na.rm=TRUE),numeric(1))),
 file.path(CFG$output_dir,"ndvi_summary.csv"),row.names=FALSE)

## 3) MODEL FUNCTION: STRATIFIED COX, ZIP ROBUST SE
fit_model <- function(x,outcome,exposures,label,weight=NULL,event=NULL){
 y <- CFG$outcomes[[outcome]]
 vars <- c(".time",y,exposures,CFG$cat,CFG$num,".stratum",".zip")
 a <- x[complete.cases(x[,vars,drop=FALSE]),,drop=FALSE]
 if(!nrow(a)) stop("No complete cases: ",label)
 a$.event <- if(is.null(event)) a[[y]] else event[match(a[[CFG$id]],names(event))]
 rhs <- c(exposures,CFG$cat,CFG$num,"strata(.stratum)")
 f <- as.formula(paste("Surv(.time,.event)~",paste(rhs,collapse="+")))
 if(!is.null(weight)){a$.w<-weight[match(a[[CFG$id]],names(weight))]
   f<-update(f,.~.+cluster(.zip)); m<-coxph(f,a,weights=.w,robust=TRUE,ties="efron")
 } else m<-coxph(f,a,cluster=.zip,robust=TRUE,ties="efron")
 b<-coef(m)[exposures]; se<-sqrt(diag(vcov(m)))[exposures]
 delta<-vapply(exposures,function(v) if(v %in% unname(CFG$ndvi)) CFG$iqr_ndvi else IQR(a[[v]],na.rm=TRUE),numeric(1))
 data.frame(outcome=outcome,model=label,N=nrow(a),events=sum(a$.event),
  exposure=exposures,contrast=delta,HR=exp(b*delta),
  LCL=exp((b-1.96*se)*delta),UCL=exp((b+1.96*se)*delta))
}
res<-list(); status<-list()

## 4) MAIN NDVI MODELS, ONE WINDOW AT A TIME AND MUTUALLY ADJUSTED
for(o in names(CFG$outcomes)){
 for(w in names(CFG$ndvi)) res[[length(res)+1]]<-tryCatch(
  fit_model(d,o,CFG$ndvi[[w]],paste0("NDVI_",w)),error=function(e){message(conditionMessage(e));NULL})
 res[[length(res)+1]]<-tryCatch(fit_model(d,o,unname(CFG$ndvi),"NDVI_mutually_adjusted"),
  error=function(e){message(conditionMessage(e));NULL})
}

## 5) MULTICOLLINEARITY OF PRECONCEPTION, PREGNANCY, POSTNATAL NDVI
z<-d[complete.cases(d[unname(CFG$ndvi)]),unname(CFG$ndvi),drop=FALSE]
vif<-sapply(seq_along(CFG$ndvi),function(i)1/(1-summary(lm(z[[i]]~.,data=z[-i]))$r.squared))
write.csv(data.frame(window=names(CFG$ndvi),VIF=vif),file.path(CFG$output_dir,"NDVI_VIF.csv"),row.names=FALSE)

## 6) EFFECT MODIFICATION AND STRATUM-SPECIFIC MODELS
for(o in names(CFG$outcomes)) for(g in CFG$subgroups) if(g %in% names(d)){
 for(level in unique(as.character(d[[g]][!is.na(d[[g]])]))){
  s<-d[!is.na(d[[g]])&as.character(d[[g]])==level,,drop=FALSE]
  for(w in names(CFG$ndvi)) res[[length(res)+1]]<-tryCatch(
   fit_model(s,o,CFG$ndvi[[w]],paste("subgroup",g,level,w,sep="_")),error=function(e)NULL)
 }
}

## 7) TREE CANOPY ALTERNATIVE EXPOSURE
tree_fields<-unname(CFG$tree)
if(all(tree_fields %in% names(d))) for(o in names(CFG$outcomes)) for(w in names(CFG$tree))
 res[[length(res)+1]]<-tryCatch(fit_model(d,o,CFG$tree[[w]],paste0("tree_",w)),error=function(e)NULL)

## 8) DIAGNOSIS AGE 3+ AND 8+ RESTRICTIONS
for(o in names(CFG$outcomes)){
 agecol<-CFG$dx_age[[o]]
 if(agecol %in% names(d)) for(a0 in c(3,8)){
   e<-setNames(as.integer(d[[CFG$outcomes[[o]]]]==1 & !is.na(d[[agecol]]) &
                            d[[agecol]]>=a0),as.character(d[[CFG$id]]))
   res[[length(res)+1]]<-tryCatch(fit_model(d,o,unname(CFG$ndvi),
     paste0("diagnosis_age_",a0,"plus"),event=e),error=function(e)NULL)
 }
}

## 9) IPCW SENSITIVITY; CENSOR EVENT MUST MEAN LOSS TO FOLLOW-UP
ipcw<-function(x){
 ff<-as.formula(paste("Surv(",CFG$censor_time,",",CFG$censor,")~",
  paste(c(CFG$cat,CFG$num,"strata(.stratum)"),collapse="+")))
 fn<-as.formula(paste("Surv(",CFG$censor_time,",",CFG$censor,")~strata(.stratum)"))
 den<-coxph(ff,x); num<-coxph(fn,x)
 wd<-exp(-predict(den,newdata=x,type="expected"))
 wn<-exp(-predict(num,newdata=x,type="expected"))
 w<-wn/wd
 w<-pmin(w,quantile(w,.99,na.rm=TRUE))
 setNames(w,as.character(x[[CFG$id]]))
}
if(all(c(CFG$censor,CFG$censor_time)%in%names(d))){
 w<-tryCatch(ipcw(d),error=function(e){message("IPCW skipped: ",conditionMessage(e));NULL})
 if(!is.null(w)) for(o in names(CFG$outcomes)){
  for(win in names(CFG$ndvi)) res[[length(res)+1]]<-tryCatch(
    fit_model(d,o,CFG$ndvi[[win]],paste0("IPCW_",win),weight=w),error=function(e)NULL)
  res[[length(res)+1]]<-tryCatch(fit_model(d,o,unname(CFG$ndvi),
    "IPCW_mutually_adjusted",weight=w),error=function(e)NULL)
 }
}

## 10) SAVE RESULTS AND RUN SETTINGS
res<-Filter(Negate(is.null),res)
if(length(res)) write.csv(do.call(rbind,res),file.path(CFG$output_dir,"all_HRs.csv"),row.names=FALSE)
dput(CFG,file=file.path(CFG$output_dir,"settings_used.R"))
capture.output(sessionInfo(),file=file.path(CFG$output_dir,"sessionInfo.txt"))
