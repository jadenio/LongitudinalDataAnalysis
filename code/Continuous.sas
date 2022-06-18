
/*
Analysis of Renal Transplant Outcome Data, Keeping Outcome Continuous
*/


data kidney_data;
set renal.renal;
patient_id=id;
run;
*rename HC0=0 HC06=0.5 HC1=1 HC2=2 HC3=3 HC4=4 HC5=5 HC6=6 HC7=7 HC8=8 HC9=9 HC10=10;
*rename HC0='0' HC06='0.5' HC1='1' HC2='2' HC3='3' HC4='4' HC5='5' HC6='6' HC7='7' HC8='8' HC9='9' HC10='10';

proc transpose data=kidney_data out=kidney_data_t prefix=haema name=char_time;
	var HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by patient_id age cardio male reject;
run;

data kidney_data_t; 
set kidney_data_t;
if char_time='HC0' then delete;
if char_time='HC06' then delete;
if char_time='HC1' then time=1;
if char_time='HC2' then time=2;
if char_time='HC3' then time=3;
if char_time='HC4' then time=4;
if char_time='HC5' then time=5;
if char_time='HC6' then time=6;
if char_time='HC7' then time=7;
if char_time='HC8' then time=8;
if char_time='HC9' then time=9;
if char_time='HC10' then time=10;
timecls=time;
run;

* Mean structures with sdevs;
/*This isn't quite right because the data is not balanced. 
There are missing measurements for some subjects at the last few timepoints
*/
ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="mean_sex" imagefmt=png;
title "Mean Response by sex";
proc sgplot data=kidney_data_t;
   vline time / response=haema1 group=male stat=mean limitstat=stderr;
   yaxis label='Mean +/- SEM';
run;

ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="mean_cardio" imagefmt=png;
title "Mean Response by cardio";
proc sgplot data=kidney_data_t;
   vline time / response=haema1 group=cardio stat=mean limitstat=stderr;
   yaxis label='Mean +/- SEM';
run;
ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="mean_age" imagefmt=png;
title "Mean Response by age";
proc sgplot data=kidney_data_t;
   vline time / response=haema1 group=age stat=mean limitstat=stderr;
   yaxis label='Mean +/- SEM';
run;
ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="mean_reject" imagefmt=png;
title "Mean Response by reject";
proc sgplot data=kidney_data_t;
   vline time / response=haema1 group=reject stat=mean limitstat=stderr;
   yaxis label='Mean +/- SEM';
run;
	
* standard linear model. from this we can get the residuals;
* using just time to get an idea of the between/within covariance;
proc mixed data = kidney_data_t;
class time patient_id;
model haema1 = time / noint solution residual outp=predresid;
run;
* age reject cardio male;

data var_func;
set predresid;
est_var=resid**2;
run;

proc means data=var_func mean;
class time;
var est_var;
output out=avg_var mean=;
run;

ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="var_func" imagefmt=png;
proc sgplot data=avg_var;
   scatter x=time y=est_var;
   series x=time y=est_var;
run;

* correlation scatter plots;
* need to use residuals but not sure how to untranspose data once residuals have been calculated;
ods graphics / groupmax=1200;
ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="corr_time" imagefmt=png;
proc sgscatter data=kidney_data;
  title "Correlation Matrix over Time";
  matrix HC0 HC2 HC4 HC6 HC8 HC10
         / group=patient_id;
run;
* note that the variogram doesn't like missing data;


* individual profiles;
proc sort data=kidney_data_t out=kidney_data_t_sort;
by patient_id time;
run;

data kidney_subset;
set kidney_data_t_sort;
if patient_id > 24 then delete;
run;

* needed to extend default number of groups to plot all lines, too many;
ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="individual_curves" imagefmt=png;
ods graphics / GROUPMAX=1200 obsmax=500;
proc sgplot data=kidney_subset;
   scatter x=time y=haema1 / group=patient_id;
   series x=time y=haema1 / group=patient_id;
run;


* Individual profiles of ols residuals;
proc sort data=predresid out=predresid_sort;
by patient_id time;
run;

data resid_subset;
set predresid_sort;
if patient_id > 24 then delete;
run;

ods graphics / GROUPMAX=1200 maxobs=500;
ods listing gpath='/home/u49868933/longitudinal/figures';
ods graphics / imagename="individual_resid" imagefmt=png;
proc sgplot data=resid_subset;
   scatter x=time y=resid / group=patient_id;
   series x=time y=resid / group=patient_id;
run;
* The above indicates random intercepts are important, but not necessarily random slopes;

/*
3. Fitting a multivariate model, and finding the most parsimonious mean structure
*/

proc mixed data = kidney_data_t method=reml;
by patient_id;
class timecls;
model haema1 = time / solution noint;
repeated timecls / type=un subject=patient_id r rcorr;
run;
* how do I get the solution B into a table?;


/*
Using two-stage analysis to get an initial impression about trends and effects
of covariates
*/
* how tf to do this?;
* starts slide 94;
* dis doesn't work...;
proc sort data=kidney_data_t out=kidney_data_t_sort;
by patient_id;
run;
proc mixed data = kidney_data_t method=reml;
class timecls cardio male reject;
by patient_id;
model haema1 = time	/ solution;
repeated timecls / r rcorr;
run;

*age reject cardio male;
*	time*age time*reject time*cardio time*male ;

/*
5. Formulate a plausible random-effects model
*/
* Look at estimated covariance matrix with unstructured covariance matrix, slide 85;
proc mixed data = kidney_data_t method=reml;
class timecls patient_id;
model haema1 = time age reject cardio male time*age time*reject time*cardio time*male / noint solution;
repeated timecls / type=un subject=patient_id r rcorr;
run;

proc mixed data = kidney_data_t method=reml;
class patient_id cardio male reject;
model haema1 = time age reject cardio male time*age time*reject time*cardio time*male / noint solution;
random intercept time / type=un subject=patient_id g gcorr;
run;

* We can decide between random and repeated models with the same mean structure by doing a LR test;
* repeated: -2 Res Log Likelihood	41578.1, 54 df from null test;
* random: 40531.6, 3 df from null
* 41578.1 - 40531.6 ~ chi(54-3) = p value of approx 0. the random model is preferred (larger negative log likelihood).