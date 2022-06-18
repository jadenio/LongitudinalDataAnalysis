*****************************************************************************************
* DICHOTOMIZATION
*****************************************************************************************;

data kidney_data;
set renal.renal;
patient_id=id;
run;

* TRANSFORMATION * ; 

data renal.renal_dummy;
set renal.renal;
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


* DICHOTOMIZATION, depending on NEGATIVE OR POSITIVE GROWTH * ; 

data renal.renal_dummy_posneg;
set renal.renal_dummy;

if missing(HC06) then H06=.;
else if HC06<=0 then H06 = 0 ;
else H06=1;

if missing(HC1) then H1=.;
else if HC1<=0 then H1 = 0 ;
else H1=1;

if missing(HC2) then H2=.;
else if HC2<=0 then H2 = 0 ;
else H2=1;

if missing(HC3) then H3=.;
else if HC3<=0 then H3 = 0 ;
else H3=1;

if missing(HC4) then H4=.;
else if HC4<=0 then H4 = 0 ;
else H4=1;

if missing(HC5) then H5=.;
else if HC5<=0 then H5 = 0 ;
else H5=1;

if missing(HC6) then H6=.;
else if HC6<=0 then H6 = 0 ;
else H6=1;

if missing(HC7) then H7=.;
else if HC7<=0 then H7 = 0 ;
else H7=1;

if missing(HC8) then H8=.;
else if HC8<=0 then H8 = 0 ;
else H8=1;

if missing(HC9) then H9=.;
else if HC9<=0 then H9 = 0 ;
else H9=1;

if missing(HC10) then H10=.;
else if HC10<=0 then H10 = 0 ;
else H10=1;

drop HC06 HC1 HC2 HC3 HC4 HC5 HC6 HC7 HC8 HC9 HC10;
run;


* Create a narrow table;

data renal.renal_dummy_posneg (drop=H06 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10);
set renal.renal_dummy_posneg;
format haematocrit comma2.;
length Yearsafter $ 10;

Yearsafter="HC06";
haematocrit = H06;
output;

Yearsafter="HC1";
haematocrit = H1;
output;

Yearsafter="HC2";
haematocrit = H2;
output;

Yearsafter="HC3";
haematocrit = H3;
output;

Yearsafter="HC4";
haematocrit = H4;
output;

Yearsafter="HC5";
haematocrit = H5;
output;

Yearsafter="HC6";
haematocrit = H6;
output;

Yearsafter="HC7";
haematocrit = H7;
output;

Yearsafter="HC8";
haematocrit = H8;
output;

Yearsafter="HC9";
haematocrit = H9;
output;

Yearsafter="HC10";
haematocrit = H10;
output;

run;

data renal.renal_dummy_posneg (drop=Yearsafter);
set renal.renal_dummy_posneg;

if Yearsafter="HC06" then months=6;
else if Yearsafter="HC1" then months=12;
else if Yearsafter="HC2" then months=24;
else if Yearsafter="HC3" then months=36;
else if Yearsafter="HC4" then months=48;
else if Yearsafter="HC5" then months=60;
else if Yearsafter="HC6" then months=72;
else if Yearsafter="HC7" then months=84;
else if Yearsafter="HC8" then months=96;
else if Yearsafter="HC9" then months=108;
else months=120;

months_c = months;
run;


ods graphics on;
proc freq data=renal.renal_dummy_posneg;
tables haematocrit*months / nocum nopercent norow nofreq plots=freqplot(orient=horizontal scale=percent);
run;

*****************************************************************************************
* TRANSITIONAL MODELS
*****************************************************************************************;

* Create /* %dropout macro */ ;

%macro dropout(data=,id=,time=,response=,out=);
%if %bquote(&data)= %then %let data=&syslast;
proc freq data=&data noprint;
tables &id /out=freqid;
tables &time / out=freqtime;
run;
proc iml;
reset noprint;
use freqid;
read all var {&id};
nsub = nrow(&id);
use freqtime;
read all var {&time};
ntime = nrow(&time);
time = &time;
use &data;
read all var {&id &time &response};
n = nrow(&response);
dropout = j(n,1,0);
ind = 1;
do while (ind <= nsub);
    j=1;
    if (&response[(ind-1)*ntime+j]=.) then print "First Measurement is Missing";
    if (&response[(ind-1)*ntime+j]^=.) then
        do;
            j = ntime;
            do until (j=1);
                if (&response[(ind-1)*ntime+j]=.) then
                    do;
                        dropout[(ind-1)*ntime+j]=1;
                            j = j-1;
                        end;
                    else j = 1;
                end;
                    end;
        ind = ind+1;
end;
prev = j(n,1,1);
prev[2:n] = &response[1:n-1];
i=1;
do while (i<=n);
if &time[i]=time[1] then prev[i]=.;
i = i+1;
end;
create help var {&id &time &response dropout prev};
append;
quit;
data &out;
merge &data help;
run;
%mend;

* Create AR(1);

%dropout(data=renal.renal_dummy_posneg, id=id, time=months, response=haematocrit, out=test2);

data test2a;
set test2;
prev1=prev;
drop prev;
run;


* In order to account for unequal spacing:

* Option 1 --> include AR(2) 
--> 1) Create AR(2), 2) Fit the Model ;

%dropout(data=test2a, id=id, time=months, response=prev1, out=test3);

data test3a;
set test3;
prev2=prev;
drop prev;
run;

* Model-fitting;

proc genmod data=test3a descending;
class male(ref=first) cardio(ref=first) reject(ref=first) months id;

model haematocrit = male age cardio reject months_c prev1 prev2/ dist=binomial;

run;

* Option 2 --> include alpha1 for all the measurements and an additional alpha1a for the first measurement;;

data test2b;
set test2a;
prev1a=.;
if months_c<=12 then prev1a=prev1;
else prev1a=0;
run; 

proc genmod data=test2b descending;
class male(ref=first) cardio(ref=first) reject(ref=first) months id;

model haematocrit = male age cardio reject months_c prev1 prev1a/ dist=binomial;

run;

** Option 3 --> include two different models (depending on whether the previous measurement is equal to 0 or 1);;

proc genmod data=test3a descending;
class male(ref=first) cardio(ref=first) reject(ref=first) months id;

model haematocrit = male*prev1 age*prev1 cardio*prev1 reject*prev1 months_c*prev1 prev1/ noint dist=binomial;

run;


*
*****************************************************************************************
* MARGINAL MODELS
*****************************************************************************************;

**************************
* GEE (Liang and Zegler)
**************************

* With exchangeable working correlation structure;

proc genmod data=renal.renal_dummy_posneg descending;
class male(ref=first) cardio(ref=first) reject(ref=first) months id;

model haematocrit = male age cardio reject months_c/ dist=binomial;

repeated subject=id / withinsubject=months type=exch covb corrw modelse;
run;

* With unstructured working correlation structure;

proc genmod data=renal.renal_dummy_posneg descending;
class male(ref=first) cardio(ref=first) reject(ref=first) months id;

model haematocrit = male age cardio reject months_c/ dist=binomial;

repeated subject=id / withinsubject=months type=unstr covb corrw modelse;
run;


****************************************************
* GEE based on linearization: Nelder and Wedderburn
****************************************************

* With unstructured working correlation structure;

proc glimmix data=renal.renal_dummy_posneg method=RSPL empirical;
class male(ref=first) reject(ref=first) cardio(ref=first) months id;

model haematocrit (event="1") = male age cardio reject months_c/ dist=binary solution;
random _residual_ / subject=id type=un;
run;

* Naive SE;
proc glimmix data=renal.renal_dummy_posneg method=RSPL;
class male(ref=first) reject(ref=first) cardio(ref=first) months id;

model haematocrit (event="1") = male age cardio reject months_c/ dist=binary solution;
random _residual_ / subject=id type=un;
run;

* With compound symmetry working correlation structure;


proc glimmix data=renal.renal_dummy_posneg method=RSPL empirical;
class male(ref=first) reject(ref=first) cardio(ref=first) months id;

model haematocrit (event="1") = male age cardio reject months_c/ dist=binary solution;
random _residual_ / subject=id type=cs;
run;

* To get naive SE;
proc glimmix data=renal.renal_dummy_posneg method=RSPL;
class male(ref=first) reject(ref=first) cardio(ref=first) months id;

model haematocrit (event="1") = male age cardio reject months_c/ dist=binary solution;
random _residual_ / subject=id type=cs;
run;



***************************************************************************
* Alternating logistic regression --> inference about correlation structure
* We don't do GEE2, because its computation is complicated;
***************************************************************************

* With exchangeable working correlation structure;

proc genmod data=renal.renal_dummy_posneg descending;
class male(ref=first) reject(ref=first) cardio(ref=first) months id;

model haematocrit = male age reject cardio months_c/ dist=binomial;

repeated subject=id / withinsubject=months logor=exch covb corrw modelse;
run;

* With unstructured working correlation structure;

proc genmod data=renal.renal_dummy_posneg descending;
class male(ref=first) reject(ref=first) cardio(ref=first) months id;

model haematocrit = male age reject cardio months_c/ dist=binomial;

repeated subject=id / withinsubject=months logor=fullclust covb corrw modelse;
run;



