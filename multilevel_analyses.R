# Rosa Cheesman 2022
# EA-PGI->achievement + school n'hood municipality random effects
# outcome = composite maths, english and reading
# note this is for 2018 merged codes

# |==================================================|
lapply(paste('package:',names(sessionInfo()$otherPkgs),sep=""),detach,character.only=TRUE,unload=TRUE)
library(lmtest)
library(lme4)
library(stargazer)
library(psych)
library(foreign)
library(haven)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)
library(tidyverse)
library(MASS)

# READ IN DATA
# EA3_ch = child EA3 PGS; EA3_par_s=parent EA3
# poeng= national test results across subjects; time = grade test taken (3 timepoints) NB data long format with a new row for each timepoint.
# komm= municipality ID
# delom = delomraade / area ID
# gkrets = neighbourhood ID
# school = school iD
# m/p PCs = prinicpal components of genetic ancestry for mothers and fathers
# |==================================================|
b<-dat_2018[,c("IID","EA3_ch_s","EA3_par_s" ,
               "time","poeng", "komm_2018","delom_2018","gkrets","school" ,
               "mPC1","mPC2","mPC3","mPC4","mPC5","pPC1","pPC2","pPC3","pPC4","pPC5",
               "par_earned_income_2014_2017","par_EduYears11_2018")]
data<-na.omit(b)

# test random intercept structure
# |==================================================|
# STEP 1: no school/geo data
base <- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+ (1|IID),REML=FALSE,data=data)
# STEP 2: add school
sch_int <- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+ (1|IID)+(1|school),REML=FALSE,data=data)
# iid is crossed because children move schools.
# STEP 3: add neighbourhood--# nb school partically crossed with neigh.
nei_int<- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+ (1|IID)+(1|gkrets)+(1|school),REML=FALSE,data=data)
# STEP 4: add delom_2018rade (nested above neigh in the data)
del_int<- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+ (1|IID)+(1|delom_2018)+(1|gkrets)+(1|school),REML=FALSE,data=data)
# STEP 5: add municipality w(of residence)
mun_int<- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+ (1|IID)+(1|komm_2018)+(1|delom_2018)+(1|gkrets)+(1|school),REML=FALSE,data=data)

# check fit
anova(base,sch_int,mun_int,del_int,nei_int)
# test random slopes
# |==================================================|
sch_slp<- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+(1|IID)+(1|komm_2018)+(1|delom_2018)+(1|gkrets)+(1+EA3_ch_s|school),REML=FALSE,data=data)
nei_slp<- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+(1|IID)+(1|komm_2018)+(1|delom_2018)+(1+EA3_ch_s|gkrets)+(1+EA3_ch_s|school),REML=FALSE,data=data)
anova(mun_int,sch_slp,nei_slp)

summary(sch_slp)

# get iccs for best model
# |==================================================|
VarCorr(sch_slp) %>% as_data_frame() %>%
  mutate(icc=vcov/sum(vcov)) %>%
  select(grp, icc)

# now control for school-level covariates in best fitting slope model
# |==================================================|
b<-dat_2018[,c("IID","EA3_ch_s","EA3_par_s" ,
               "time","poeng", "komm_2018","delom_2018","gkrets","school" ,
               "mPC1","mPC2","mPC3","mPC4","mPC5","pPC1","pPC2","pPC3","pPC4","pPC5",
               "sch_avg_pedu10","sch_avg_prank6_9","sch_avg_nonwest",
               "orgnr5ginipeduc10","orgnr5giniprank6_9",
               "par_earned_income_2014_2017","par_EduYears11_2018")]
data<-na.omit(b)
baseEA3_ch_s    <- lmer(poeng ~ EA3_ch_s+EA3_par_s+par_earned_income_2014_2017+par_EduYears11_2018+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+time+ (1|IID)+(1|komm_2018)+(1|delom_2018)+(1|gkrets)+(1+EA3_ch_s|school),REML=FALSE,data=data)
sch_covsEA3_ch_s <- lmer(poeng ~ EA3_ch_s+EA3_par_s+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+par_earned_income_2014_2017+par_EduYears11_2018+sch_avg_pedu10+sch_avg_prank6_9+sch_avg_nonwest+orgnr5ginipeduc10+orgnr5giniprank6_9+time+ (1|IID)+(1|komm_2018)+(1|delom_2018)+(1|gkrets)+(1+EA3_ch_s|school),REML=FALSE,data=data)
sch_covs_slpEA3_ch_s<- lmer(poeng ~ EA3_ch_s+EA3_par_s+mPC1+mPC2+mPC3+mPC4+mPC5+pPC1+pPC2+pPC3+pPC4+pPC5+par_earned_income_2014_2017+par_EduYears11_2018+sch_avg_pedu10+sch_avg_prank6_9+sch_avg_nonwest+orgnr5ginipeduc10+orgnr5giniprank6_9+EA3_ch_s*sch_avg_pedu10+EA3_ch_s*sch_avg_prank6_9+EA3_ch_s*sch_avg_nonwest+EA3_ch_s*orgnr5ginipeduc10+EA3_ch_s*orgnr5giniprank6_9+time+(1|IID)+(1|komm_2018)+(1|delom_2018)+(1|gkrets)+(1+EA3_ch_s|school),REML=FALSE,data=data)

anova(baseEA3_ch_s,sch_covsEA3_ch_s,sch_covs_slpEA3_ch_s)

summary(sch_covsEA3_ch_s)


# plot each school line.
# |==================================================|

# intercept values from model (for average PGS)
randoms<-ranef(sch_slp, condVar = TRUE)
m_int_slp<-as.data.frame(randoms$school)

colnames(m_int_slp)<-c("Intercept","Slope")
head(m_int_slp)

df <- tibble::rownames_to_column(m_int_slp, "school")
df2<-merge(data, df, by = 'school')

fe<-fixef(sch_slp)
df2$pred<-df2$Intercept+df2$Slope * df2$EA3_ch_s+fe[1]+fe[2]*df2$EA3_ch_s 

df2<-df2 %>%
  mutate(weakest_slope = as.numeric(Slope >= quantile(Slope, 0.975)), 
         strongest_slope = as.numeric(Slope <= quantile(Slope, 0.025)))

df2$select_slps<-ifelse(df2$strongest_slope==1, 'top',
                        ifelse(df2$weakest_slope==1, 'bottom', 'normal'))

# PLOT DATA FOR ALL SCHOOLS
theme_set(theme_bw()) 
ggplot(data = df2, aes(x = EA3_ch_s, y=pred,group=school,color=select_slps))+
  geom_line(alpha=.4)+
  xlab("Within-family EA-PGS")+
  xlim(-4,4)+
  ylim(-2,2)+
  ylab("Predicted Achievement")+
  theme(plot.title = element_text(color="black", size=16, face="bold",hjust = 0.5),
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        strip.text = element_text(size = 16))+
  geom_vline(xintercept=0)+
  scale_color_manual(values=c("blue", "#999999", "red"))+
  theme(legend.position = "none")


# plot variation in school effect by PGS
# |==================================================|
# use the par SES controlled versions for each predictor


# Returns the intercept variance for given value of x
# x = covariate value
# v_u0 = intercept variance at x = 0
# v_u1 = slope variance
# c_u0u1 = intercept-slope covariance at x = 0
v_u0_x = function(x, v_u0, v_u1, c_u0u1) {
  v_u0 + 2 * c_u0u1 * x + v_u1 * x^2 
}

# Returns the intercept-slope correlation for given value of x
# x = covariate value
# v_u0 = intercept variance at x = 0
# v_u1 = slope variance
# c_u0u1 = intercept-slope covariance at x = 0
r_u0u1_x = function(x, v_u0, v_u1, c_u0u1) {
  (c_u0u1 + x * v_u1) / sqrt(v_u0_x(x, v_u0, v_u1, c_u0u1) * v_u1)
}

# Returns the % variance explained by intercept for given value of x
# x = covariate value
# v_u0 = intercept variance at x = 0
# v_u1 = slope variance
# c_u0u1 = intercept-slope covariance at x = 0
# v_e = residual variance
# rc added: v_iid is inidividual variance at x=0
frac_u0_x = function(x, v_u0, v_u1, c_u0u1, v_iid,v_gk,v_dl,v_km,v_e) {
  var_u0 = v_u0_x(x, v_u0, v_u1, c_u0u1)
  var_u0 / (var_u0 + v_iid+v_dl+v_gk+v_km+v_e)
}
fit0<-sch_slp

VC_est = as.data.frame(VarCorr(fit0))

c11 = VC_est$vcov[3]
c22 = VC_est$vcov[4]
c21 = VC_est$vcov[5]
eps = VC_est$vcov[8]
iid = VC_est$vcov[1]
gkt = VC_est$vcov[2]
del = VC_est$vcov[6]
kom = VC_est$vcov[7]

# % Intercept variance
curve(frac_u0_x(x, c11, c22, c21, iid,gkt,del,kom,eps), 
      from = -4, to = 3.7,
      xlab = "Within-family EA-PGS",
      ylab = "% Variance in achievement due to school",
      col='red')




