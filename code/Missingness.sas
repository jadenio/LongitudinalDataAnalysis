*Longitudinal Data Analysis: Inspection of the Missing Data Mechanism;


data kidney_data;
set renal.renal;
patient_id=id;
run;


* Look at pattern of missing data before transpose;
ods select MissPattern;
proc mi data=kidney_data nimpute=0;
var HC0--HC10;
run;

* Sensitivity analysis;

* from https://www.lexjansen.com/wuss/2016/76_Final_Paper_PDF.pdf;
* For the missing values, say well these values should be;
* some transformation different than the observed;
* Say those with missing values had some unmeasure reaction to the;
* surgery that was not captured;
* If the estimates do end up being different, than we have evidence against MAR;
* by default will use mcmc;

/* Here we include multiple scenarios. One is that the overall haema1 level is lower.
for those who droput. 
Another is that those who dropout experienced symptoms of rejection or cardio issues.
Another scenario is that those who dropout are older. 
A final one could be that those who dropout are sampled from the outcome lower than baseline 
(binary_haema=0)
We also can incorporate NCMV.
*/

* Use fcs reg because non-monotone;
ods exclude all;
proc mi data=kidney_data seed=13 nimpute=10 out=mnar_mi_scale;
 fcs reg;
 mnar adjust(HC0 / scale=0.7);
 mnar adjust(HC06 / scale=0.7);
 mnar adjust(HC1 / scale=0.7);
 mnar adjust(HC2 / scale=0.7);
 mnar adjust(HC3 / scale=0.7);
 mnar adjust(HC4 / scale=0.7);
 mnar adjust(HC5 / scale=0.7);
 mnar adjust(HC6 / scale=0.7);
 mnar adjust(HC7 / scale=0.7);
 mnar adjust(HC8 / scale=0.7);
 mnar adjust(HC9 / scale=0.7);
 mnar adjust(HC10 / scale=0.7);
 var reject male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

* Then run transformation;
data mnar_z_scale;
set mnar_mi_scale;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=mnar_z_scale out=mnar_t_scale prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data mnar_months_scale; 
set mnar_t_scale;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data mnar_dichot_scale;
set mnar_months_scale;
if haema1<=0 then binary_haema = 0 ;
else binary_haema=1;
run;

proc gee data=mnar_dichot_scale;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* delete annoying params;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

* for output from GEE;
ods select VarianceInfo ParameterEstimates;
title "Missing pattern: scale=0.7";
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
 modeleffects male cardio reject months age;
run;


* rejection symptoms droput model;
* control-based pattern imputation;
* Don't include reject in var;
* have to include ods exclude everytime after I run ods select;
ods exclude all;
proc mi data=kidney_data seed=13 nimpute=10 out=mnar_mi_reject;
 class reject;
 fcs reg;
 mnar model(HC0 / modelobs=(reject='1'));
 mnar model(HC06 / modelobs=(reject='1'));
 mnar model(HC1 / modelobs=(reject='1'));
 mnar model(HC2 / modelobs=(reject='1'));
 mnar model(HC3 / modelobs=(reject='1'));
 mnar model(HC4 / modelobs=(reject='1'));
 mnar model(HC5 / modelobs=(reject='1'));
 mnar model(HC6 / modelobs=(reject='1'));
 mnar model(HC7 / modelobs=(reject='1'));
 mnar model(HC8 / modelobs=(reject='1'));
 mnar model(HC9 / modelobs=(reject='1'));
 mnar model(HC10 / modelobs=(reject='1'));
  var male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;


* Then run transformation;
data mnar_z_reject;
set mnar_mi_reject;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=mnar_z_reject out=mnar_t_reject prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data mnar_months_reject; 
set mnar_t_reject;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data mnar_dichot_reject;
set mnar_months_reject;
if haema1<=0 then binary_haema = 0 ;
else binary_haema=1;
run;

proc gee data=mnar_dichot_reject;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* delete annoying params;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

* for output from GEE;
ods select VarianceInfo ParameterEstimates;
title "Missing pattern: sample where reject=1";
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
 modeleffects male cardio reject months age;
run;


* cardio symptoms droput model;
ods exclude all;
proc mi data=kidney_data seed=13 nimpute=10 out=mnar_mi_cardio;
 class cardio;
 fcs reg;
 mnar model(HC0 / modelobs=(cardio='1'));
 mnar model(HC06 / modelobs=(cardio='1'));
 mnar model(HC1 / modelobs=(cardio='1'));
 mnar model(HC2 / modelobs=(cardio='1'));
 mnar model(HC3 / modelobs=(cardio='1'));
 mnar model(HC4 / modelobs=(cardio='1'));
 mnar model(HC5 / modelobs=(cardio='1'));
 mnar model(HC6 / modelobs=(cardio='1'));
 mnar model(HC7 / modelobs=(cardio='1'));
 mnar model(HC8 / modelobs=(cardio='1'));
 mnar model(HC9 / modelobs=(cardio='1'));
 mnar model(HC10 / modelobs=(cardio='1'));
 var reject male age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;


* Then run transformation;
data mnar_z_cardio;
set mnar_mi_cardio;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=mnar_z_cardio out=mnar_t_cardio prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data mnar_months_cardio; 
set mnar_t_cardio;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data mnar_dichot_cardio;
set mnar_months_cardio;
if haema1<=0 then binary_haema = 0 ;
else binary_haema=1;
run;

proc gee data=mnar_dichot_cardio;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* delete annoying params;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

* for output from GEE;
ods select VarianceInfo ParameterEstimates;
title "Missing pattern: sample where cardio=1";
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
 modeleffects male cardio reject months age;
run;

* Dropout model where hc06=0 in the dichotomized data version;
* Wide transformation and dichot;
data mnar_z_1;
set kidney_data;
HC06_z = log(HC06/HC0);
run;

data;
set mnar_z_1;
if HC06_z = . then output;
run;

* As previously observed, there are no missing values for hc06 which is why I can do this;
data wide_dichot;
set mnar_z_1;
if HC06_z<=. then HC06_di = .;
else if HC06_z<=0 then HC06_di = 0;
else HC06_di=1;
run;

ods exclude all;
proc mi data=wide_dichot seed=13 nimpute=10 out=mnar_mi_dichot;
 class HC06_di;
 fcs reg;
 mnar model( HC0 / modelobs= (HC06_di='0'));
 mnar model(HC06 / modelobs=(HC06_di='0'));
 mnar model(HC1 / modelobs=(HC06_di='0'));
 mnar model(HC2 / modelobs=(HC06_di='0'));
 mnar model(HC3 / modelobs=(HC06_di='0'));
 mnar model(HC4 / modelobs=(HC06_di='0'));
 mnar model(HC5 / modelobs=(HC06_di='0'));
 mnar model(HC6 / modelobs=(HC06_di='0'));
 mnar model(HC7 / modelobs=(HC06_di='0'));
 mnar model(HC8 / modelobs=(HC06_di='0'));
 mnar model(HC9 / modelobs=(HC06_di='0'));
 mnar model(HC10 / modelobs=(HC06_di='0'));
 var reject cardio male age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

* Have to redo the transposition and transformation after the imputation;
data mnar_z_dichot;
set mnar_mi_dichot;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0 HC06_di hc06_z;
run;

* Then transpose wide to long for modeling;
proc transpose data=mnar_z_dichot out=mnar_t_dichot prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data mnar_months_dichot; 
set mnar_t_dichot;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data mnar_dichot_dichot;
set mnar_months_dichot;
if haema1<=0 then binary_haema = 0 ;
else binary_haema=1;
run;

proc gee data=mnar_dichot_dichot;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* delete annoying params;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

* for output from GEE;
ods select VarianceInfo ParameterEstimates;
title "Missing pattern: sample from binary HC06=0";
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
 modeleffects male cardio reject months age;
run;



* Scale haema only for reject=1 patients;
ods exclude all;
proc mi data=kidney_data seed=13 nimpute=10 out=mnar_mi_scale_reject;
 class reject;
 fcs reg;
 mnar adjust(HC0 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC06 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC1 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC2 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC3 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC4 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC5 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC6 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC7 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC8 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC9 / scale=0.5 adjustobs=(reject='1'));
 mnar adjust(HC10 / scale=0.5 adjustobs=(reject='1'));
 var male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

data mnar_z_scale_reject;
set mnar_mi_scale_reject;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=mnar_z_scale_reject out=mnar_t_scale_reject prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data mnar_months_scale_reject; 
set mnar_t_scale_reject;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data mnar_dichot_scale_reject;
set mnar_months_scale_reject;
if haema1<=0 then binary_haema = 0 ;
else binary_haema=1;
run;

proc gee data=mnar_dichot_scale_reject;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* delete annoying params;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

* for output from GEE;
ods select VarianceInfo ParameterEstimates;
title "Missing pattern: scale=.05, reject=1";
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
 modeleffects male cardio reject months age;
run;

* MI with NCMV;
ods exclude all;
proc mi data=kidney_data seed=13 simple out=kidney_mono nimpute=10;
mcmc impute=monotone;
var reject male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

proc mi data=kidney_mono seed=13 simple out=mnar_mi_ncmv nimpute=1;
monotone reg;
mnar model( HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10 / modelobs=ncmv);
var reject male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

data mnar_z_ncmv;
set mnar_mi_ncmv;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=mnar_z_ncmv out=mnar_t_ncmv prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data mnar_months_ncmv; 
set mnar_t_ncmv;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data mnar_dichot_ncmv;
set mnar_months_ncmv;
if haema1<=0 then binary_haema = 0 ;
else binary_haema=1;
run;

proc gee data=mnar_dichot_ncmv;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* delete annoying params;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

* for output from GEE;
ods select VarianceInfo ParameterEstimates;
title "Missing pattern: NCMV restriction";
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
 modeleffects male cardio reject months age;
run;


* For future reference, example of tipping point analysis;
* https://documentation.sas.com/doc/en/statcdc/14.2/statug/statug_mianalyze_examples13.htm;

*Code for Weighted GEE, assuming MAR;
* First need to generate monotone data;
ods exclude all;
proc mi data=kidney_data seed=13 simple nimpute=10 round=0.1 out=kidney_monotone;
mcmc impute = monotone;
var reject male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

data monotone_z;
set kidney_monotone;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=monotone_z out=monotone_t prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data monotone_months; 
set monotone_t;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data monotone_dichot;
set monotone_months;
if haema1=. then binary_haema=.;
else if haema1<=0 then binary_haema = 0;
else binary_haema=1;
run;

proc sort data=monotone_dichot out=monotone_dichot;
by patient_id;
run;

* Create lagged variable;
data monotone_lag;
set monotone_dichot;
by patient_id;
prev_haema = lag(binary_haema);
if months=6 then prev_haema = binary_haema;
run;

* Need to resort;
proc sort data=monotone_lag out=monotone_lag;
by _imputation_;
run;

* slide 648;
proc gee data=monotone_lag plots=histogram;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
missmodel prev_haema age reject cardio male / type=obslevel;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* get rid of annoying parms;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

ods select VarianceInfo ParameterEstimates;
title 'Weighted GEE on monotone data with multiple imputation';
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
modeleffects male cardio reject months age;
run;

* Standard GEE + MI, with fcs reg;
ods exclude all;
proc mi data=kidney_data seed=13 simple nimpute=10 round=0.1 out=kidney_gee;
fcs reg;
var reject male cardio age HC0 HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;

data gee_z;
set kidney_gee;
HC06 = log(HC06/HC0);
HC1 = log(HC1/HC0);
HC2 = log(HC2/HC0);
HC3 = log(HC3/HC0);
HC4 = log(HC4/HC0);
HC5 = log(HC5/HC0);
HC6 = log(HC6/HC0);
HC7 = log(HC7/HC0);
HC8 = log(HC8/HC0);
HC9 = log(HC9/HC0);
HC10 = log(HC10/HC0);
drop HC0;
run;

* Then transpose wide to long for modeling;
proc transpose data=gee_z out=gee_t prefix=haema name=char_time;
	var HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
	by _imputation_ patient_id age cardio male reject;
run;


data gee_months; 
set gee_t;
if char_time='HC06' then do time=.5; months=6; end;
if char_time='HC1' then do time=1; months=12; end;
if char_time='HC2' then do time=2; months=24; end;
if char_time='HC3' then do time=3; months=36; end;
if char_time='HC4' then do time=4; months=48; end;
if char_time='HC5' then do time=5; months=60; end;
if char_time='HC6' then do time=6; months=72; end;
if char_time='HC7' then do time=7; months=84; end;
if char_time='HC8' then do time=8; months=96; end;
if char_time='HC9' then do time=9; months=108; end;
if char_time='HC10' then do time=10; months=120; end;
timecls=time;
run;

* Dichotomization with positive negative slope;
data gee_dichot;
set gee_months;
if haema1=. then binary_haema=.;
else if haema1<=0 then binary_haema = 0;
else binary_haema=1;
run;

proc sort data=gee_dichot out=gee_dichot;
by patient_id;
run;

* Create lagged variable;
data gee_lag;
set gee_dichot;
by patient_id;
prev_haema = lag(binary_haema);
if months=6 then prev_haema = binary_haema;
run;

* Need to resort;
proc sort data=gee_lag out=gee_lag;
by _imputation_;
run;

* slide 648;
proc gee data=gee_lag;
by _imputation_;
class cardio reject male timecls patient_id;
model binary_haema = months age reject cardio male / dist=bin;
repeated subject=patient_id / within=timecls corr=un ecovb modelse;
ods output GEEEmpPEst=gmparms parminfo=gmpinfo modelinfo=modelinfo GEERCov=gmcovb;
run;

* get rid of annoying parms;
data gmpinfo;
set gmpinfo;
if parameter='Prm5' then delete;
if parameter='Prm7' then delete;
if parameter='Prm9' then delete;
run;

data gmparms;
set gmparms;
if level1 = 1 then delete;
run;

ods select VarianceInfo ParameterEstimates;
title 'Standard GEE on + FCS multiple imputation';
proc mianalyze parms=gmparms parminfo=gmpinfo covb=gmcovb;
modeleffects male cardio reject months age;
run;




