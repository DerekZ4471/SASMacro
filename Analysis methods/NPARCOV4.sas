*************************************************************************************************************
*************************************************************************************************************
*************************************************************************************************************

NParCov4
Non-Parametric Randomization-Based Analysis of Covariance

Richard C. Zink, Ph.D.
JMP Life Sciences
SAS Institute, Inc.

and 

Gary G. Koch, Ph.D.
Biometric Consulting Laboratory
Department of Biostatistics
University of North Carolina at Chapel Hill

Copyright 2001 - 2017

Inputs:

outcomes      List of outcome variables. Numeric.  At least 1 required. Need to be (0,1) outcomes for 
              LOGISTIC, PODDS, LOGRANK, WILCOXON

covars        List of covariates.  Numeric

exposures     Numeric exposure variables for (0,1) outcome for WILCOXON, LOGRANK and INCDENS. 
              Need same number of variables as OUTCOMES.

strata        Variable that defines strata.  Numeric.  NONE (default) for no 
              stratification (i.e. a single stratum)

trtgrps       Defines the treatments. Numeric. Differences defined 
              are higher number (trt2) - lower number (trt1).  2 treatments required.

hypoth        NULL (default), ALT.  NULL groups treatments for estimating variance.  
              ALT estimates variance within treatment

transform     NONE, LOGISTIC, PODDS, LOGRANK, LOGRATIO, WILCOXON, INCDENS.  Default is NONE. 

combine       NONE, LAST, FIRST, PRETRANSFORM.  Only usable when strata is not NONE Default is NONE.
 
c             Used in calculation of strata weights. 0 <= c <= 1. 0 assigns equal 
              weights to strata, 1 is for Mantel-Haenszel weights. Default is 1.       

alpha         For CI (0.05 default)

seed          Seed for exact p-values and confidence intervals.  (Default 0)
 
nreps         Number of random data sets to generate.  (Default 1000)

exact         NO (Default), YES. Must be NO for TRANSFORM = PODDS.  Only available
              for a single outcome. Under hypoth = ALT performs sampling with replacement within treatment 
              (bootstrap). Under hypoth = NULL performs sampling without replacement, essentially shuffling 
              treatment codes among observations (permutation).

symsize       Define space requirements for IML (Default 20000)

dsnin         Input Data (required)

dsnout        Output data (required)

details       YES (default), NO.Print analysis details to the log.

*************************************************************************************************************
*************************************************************************************************************
************************************************************************************************************;

%macro NPARCOV4( outcomes = , 
                   covars = , 
                   strata = %str(NONE), 
                  trtgrps = , 
                   hypoth = %str(NULL), 
                transform = %str(NONE),
                  combine = %str(NONE),
                        c = 1, 
                    dsnin = ,
                   dsnout = ,
                    alpha = 0.05,
                exposures =,
                     seed = 0,
                    nreps = 1000,
                  symsize = 20000,
                    exact = %str(NO),
                  details = %str(YES)
               );

   %global numstrat numtrts numobstot numresp numcov trt1 trt2 workoutcomes missind;

****************;
****************;
****************;

   %let outcomes = %upcase(&outcomes);
   %let strata = %upcase(&strata);
   %let exposures = %upcase(&exposures);
   %let covars = %upcase(&covars);

   %if ^%length(&dsnin) or ^%length(&dsnout) %then %do;
      %put ERROR: DSNIN and DSNOUT must be specified;
      %goto stopmac;
   %end;

   %if (%upcase(&hypoth) ne NULL) and (%upcase(&hypoth) ne ALT) %then %do;
      %put ERROR: HYPOTH must be specified as NULL or ALT;
      %goto stopmac;
   %end;

   %if (%upcase(&combine) ne NONE) and (%upcase(&strata) = NONE) %then %do;
      %put ERROR: STRATA = NONE so COMBINE should = NONE;
      %goto stopmac;
   %end;

   %if ^%index(NONE LOGISTIC PODDS LOGRATIO LOGRANK WILCOXON INCDENS, %upcase(&transform)) %then %do;
      %put ERROR: TRANSFORM must be one of NONE LOGISTIC PODDS LOGRATIO LOGRANK WILCOXON INCDENS;
      %goto stopmac;
   %end;
    
   %if ^%index(NONE LAST FIRST PRETRANSFORM, %upcase(&combine)) %then %do;
      %put ERROR: COMBINE must be one of NONE LAST FIRST PRETRANSFORM;
      %goto stopmac;
   %end;

   %if ^%length(&outcomes) %then %do;
      %put ERROR: OUTCOMES must have at least one variable specified;
      %goto stopmac;
   %end;

   %if ^%length(&trtgrps) %then %do;
      %put ERROR: TRTGRPS must be specified;
      %goto stopmac;
   %end;
   
   %if %sysevalf(&c < 0) or %sysevalf(&c > 1) %then %do;
      %put ERROR: C must be between 0 and 1, inclusive;
      %goto stopmac;
   %end;

   %if %sysevalf(&alpha < 0) or %sysevalf(&alpha > 1) %then %do;
      %put ERROR: ALPHA must be between 0 and 1;
      %goto stopmac;
   %end;

   %if ^%index(LOGRANK WILCOXON INCDENS, %upcase(&transform)) %then %do;
      %if %length(&exposures) ne 0 %then %do;
         %put ERROR: EXPOSURES should only be specified for WILCOXON LOGRANK INCDENS;
         %goto stopmac;
      %end;
   %end;

   %if %index(LOGRANK WILCOXON INCDENS, %upcase(&transform)) %then %do;
      %if ^%length(&exposures) %then %do;
         %put ERROR: EXPOSURES should be specified for WILCOXON LOGRANK INCDENS;
         %goto stopmac;
      %end;
   %end;

   %if ^%index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) and %upcase(&combine) = PRETRANSFORM %then %do;
      %put ERROR: COMBINE = PRETRANSFORM only available for LOGRATIO LOGISTIC PODDS INCDENS;
      %goto stopmac;
   %end;

   %if %upcase(&strata) ne NONE and %upcase(&combine) = NONE %then %do;
      %put ERROR: STRATA specified but COMBINE = NONE;
      %goto stopmac;
   %end;

   %if ^%index(NO YES, %upcase(&exact)) %then %do;
      %put ERROR: EXACT must be one of NO YES;
      %goto stopmac;
   %end;
   
   %if %upcase(&transform) = PODDS and %upcase(&exact) = YES %then %do;
      %put ERROR: PODDS must have EXACT = NO;
      %goto stopmac;
   %end;   

****************;
****************;
****************;
   
   ods noresults;

   data _null_;
      call symput("numresp", left(countw(compbl("&outcomes"), " ")));
      call symput("numcov", left(countw(compbl("&covars"), " ")));
      %if %index(LOGRANK WILCOXON INCDENS, %upcase(&transform)) %then %do;
         call symput("numtime", left(countw(compbl("&exposures"), " "))); 
      %end;
      call symput("numstratwrd", left(countw(compbl("&strata"), " ")));
   run;

   %if %upcase(&exact) = YES and &numresp > 1 %then %do;
      %put ERROR: OUTCOMES can only have one variable for EXACT = YES;
      %goto stopmac;
   %end;

   %if %sysevalf(&numstratwrd ne 1) %then %do;
      %put ERROR: STRATA can only have a single variable specified;
      %goto stopmac;
   %end; 

   %if %index(LOGRANK WILCOXON INCDENS, %upcase(&transform)) %then %do;
      %if &numtime ne &numresp %then %do;
         %put ERROR: EXPOSURES and RESPONSE must have same number of variables specified for WILCOXON LOGRANK INCDENS;
         %goto stopmac;
      %end;
   %end;

   %if %index(LOGRANK WILCOXON, %upcase(&transform)) %then %do;
      %let workoutcomes = ;
      %do i = 1 %to &numresp; 
         %let tempoutcomes = &transform._%scan(&outcomes, &i);
         %let workoutcomes = &workoutcomes &tempoutcomes; 
      %end;
   %end;
   %else %if %index(INCDENS, %upcase(&transform)) %then %do;
      %let workoutcomes = ;
      %do i = 1 %to &numresp; 
         %let tempoutcomes = INC_%scan(&outcomes, &i);
         %let workoutcomes = &workoutcomes &tempoutcomes; 
      %end;
   %end;
   %else %do;
      %let workoutcomes = &outcomes;
   %end;

   data _work_temp;
      set &dsnin end = eof;
      %if %upcase(&strata) = NONE or %upcase(&combine) = NONE %then %do;
         _NONE = 1;
      %end;
      if eof then call symput("numobstot", left(_n_));
   run;

   %if %index(LOGRANK WILCOXON PODDS LOGISTIC, %upcase(&transform)) %then %do;      
      %do i = 1 %to &numresp; 
         proc freq noprint data = _work_temp;
            tables %scan(&outcomes,&i) / out = binarychk;
         run;

         proc means n sum noprint data = binarychk;
            var %scan(&outcomes,&i);
            output out = binarychk n = nbinchk sum = sumbinchk;
         run;

         data _null_;
            set binarychk;
            call symput("nbinchk", left(nbinchk)); 
            call symput("sumbinchk", left(sumbinchk)); 
         run;

         %if not(&nbinchk = 2 and &sumbinchk = 1) %then %do;
            %put ERROR: %upcase(%scan(&outcomes,&i)) not a (0,1) OUTCOME;
            %goto stopmac;
         %end;
      %end;
       
      proc datasets nolist;
         delete binarychk;
      quit;   
   %end;

   %if %upcase(&strata) = NONE %then %do;
      %MISSDAT(varlist = &outcomes &covars &trtgrps &exposures, data = _work_temp);
   %end;
   %else %do;
      %MISSDAT(varlist = &outcomes &covars &strata &trtgrps &exposures, data = _work_temp);
   %end;
   %if &missind %then %do;
      %put ERROR: Missing values are present in at least one of OUTCOMES COVARS STRATA TRTGRPS EXPOSURES;
      %goto stopmac;
   %end;   

   %if %index(INCDENS, %upcase(&transform)) %then %do;
      data _work_temp; 
         set _work_temp(rename = (%do i = 1 %to &numresp; %scan(&outcomes, &i) = %scan(&workoutcomes, &i) %end;)); 
      run;
   %end;

   proc freq nlevels data = _work_temp;
      ods exclude all;
      ods output NLevels = _work_ntrts_chk;
      tables &trtgrps;
   run;
   ods select all;

   data _null_;
      set _work_ntrts_chk;
      call symput("numtrtlevels", left(nlevels)); 
   run; 

   %if &numtrtlevels ne 2 %then %do;
      %put ERROR: TRTGRPS must have TWO levels;
      %goto stopmac;
   %end;
 
   %if %upcase(&strata) = NONE or %upcase(&combine) = NONE %then %do;
      %let strata = _NONE;
   %end;

   proc freq noprint data = _work_temp;
      tables &trtgrps / out = _work_ntrts;
   run;

   data _null_;
      set _work_ntrts;
      if _n_ = 1 then call symput("trt1", trim(left(&trtgrps)));
      else if _n_ = 2 then call symput("trt2", trim(left(&trtgrps)));
   run;

   %if %upcase(&transform) = LOGRANK or %upcase(&transform) = WILCOXON %then %do;
      %SURVSCORE(eventlist = &outcomes, timelist = &exposures, score = %upcase(&transform));
   %end;

   proc sort data = _work_temp;
      by &strata &trtgrps;
   run;
   
   proc sql noprint;
      select distinct &strata 
      into :stratalist separated by ','
      from _work_temp;
   quit; 

   proc freq nlevels data = _work_temp;
      ods exclude all;
      ods output NLevels = _work_nstrt;
      tables &strata;
   run;   
   ods select all;

   data _null_;
      set _work_nstrt;
      call symput("numstrat", left(nlevels)); 
   run; 

   proc freq noprint data = _work_temp;
       tables &trtgrps * &strata / out = _work_trtstrat_chk(drop = percent);
   run;

   proc means min n noprint data = _work_trtstrat_chk;
      var count;
      where count > 0;
      output out = _work_trtstrat_chk n = trtstratcount min = trtstratmin;
   run;

   data _null_;
      set _work_trtstrat_chk;
      call symput("treatin", left(trtstratcount)); 
      call symput("treatinmin", left(trtstratmin)); 
   run; 

   %if &treatin ^= %sysevalf(2 * &numstrat) %then %do;
      %put ERROR: Both TRTGRPS should be present in all strata;
      %goto stopmac;
   %end;

   %if &treatinmin = 1 and %upcase(&hypoth) = ALT %then %do;
      %put ERROR: Need at least two observations for each strata*treatment combination for variance estimation under HYPOTH = ALT;
      %goto stopmac;
   %end;

   proc datasets nolist;
      delete _work_trtstrat_chk;
   quit;
   
   %if %upcase(&details) = YES %then %do;
      %put NOTE NOTE NOTE:  Treatment1 = &trt1 and Treatment2 = &trt2;
      %put NOTE NOTE NOTE:  Total Number of OBSERVATIONS for BOTH TREATMENTS: &numobstot;
      %put NOTE NOTE NOTE:  Number of RESPONSE VARIABLES: &numresp;
      %put NOTE NOTE NOTE:  Number of COVARIATES: &numcov;
      %put NOTE NOTE NOTE:  Number of STRATA: &numstrat;
   %end;

   %if %upcase(&exact) = YES and %upcase(&hypoth) = ALT %then %do;   
      proc sort data = _work_temp;
         by &trtgrps;
      run;
   %end;

   %if %upcase(&exact) = YES %then %do;
      proc multtest %if %upcase(&hypoth) = NULL %then %do; permutation %end; %else %do; bootstrap %end; 
         nocenter noprint seed = &seed nsample = &nreps outsamp = _work_outsamp(drop = _obs_) data = _work_temp;
         class &trtgrps;
         strata &strata;
         %if %upcase(&hypoth) = ALT %then %do;
            ***** select samples within treatment under alternative *****;
            by &trtgrps;
         %end;
          %if %upcase(&transform) = INCDENS %then %do; 
            test mean (&workoutcomes &exposures &covars);
          %end;
          %else %do;
            test mean (&workoutcomes &covars);
          %end;
      run;

      data _work_outsamp(drop = _class_ _stratum_);
         set _work_outsamp %if %upcase(&hypoth) = ALT %then %do; (drop = &trtgrps) %end;;
         &trtgrps = _class_ + 0;
         &strata = _stratum_ + 0;
      run;

      proc sort data = _work_outsamp;
         by _sample_ &strata &trtgrps;
      run;
   %end;
   
   %if %upcase(&exact) = YES and %upcase(&hypoth) = ALT %then %do;
      proc sort data = _work_temp;
         by &strata &trtgrps;
      run;

      data _work_jack(drop = i);
         do i = 1 to &numobstot;
            _sample_ = i + &nreps;
            do j = 1 to &numobstot;
               if j ne i then do;
                  set _work_temp point = j;
                  output;
                  end;
               end;
            end;
         stop;
      run;
   %end;   

   data _work_temp;
      set _work_temp(in = ina)
      %if %upcase(&exact) = YES %then %do; 
         _work_outsamp(in = inb) 
         %if (%upcase(&exact) = YES and %upcase(&hypoth) = ALT) %then %do;
            _work_jack(in = inc)
         %end; 
      %end;       
      ;
      if ina then _sample_ = 0;
   run;

   ***** Means *****;

   proc sort data = _work_temp;
      by _sample_ &strata &trtgrps;
   run;

   proc means mean noprint data = _work_temp;
      %if %upcase(&transform) = INCDENS %then %do; 
         var &workoutcomes &exposures &covars;
      %end;
      %else %do;
         var &workoutcomes &covars;
      %end;
      by _sample_ &strata &trtgrps;
      %if %upcase(&transform) = INCDENS %then %do; 
         output out = _work_means mean = &workoutcomes &exposures &covars;
      %end;
      %else %do;
         output out = _work_means mean = &workoutcomes &covars;
      %end;
   run;
   
   ***** Covariances *****;

   proc corr cov outp = _work_covar(where = (_type_ = 'COV')) noprint data = _work_temp;
      %if %upcase(&transform) = INCDENS %then %do; 
         var &workoutcomes &exposures &covars;
      %end;
      %else %do;
         var &workoutcomes &covars;
      %end;
      %if %upcase(&hypoth) = ALT %then %do;
         by _sample_ &strata &trtgrps;
      %end;
      %else %if %upcase(&hypoth) = NULL %then %do;
         by _sample_ &strata;         
      %end;
   run;
      
   %if %upcase(&hypoth) = NULL %then %do;
      data _work_covar;
         set _work_covar(in = ina) _work_covar(in = inb);
         if ina then &trtgrps = "&trt1";
         else if inb then &trtgrps = "&trt2";
      run;

      proc sort data = _work_covar;
         by _sample_ &strata &trtgrps;
      run;
   %end;

****************;
****************;
****************;

*** Output Data ***;

%if &numcov >= 1 %then %do;
   data _&dsnout._covtest;
      length type $ 10 trt1 trt2 strata covariates outcomes exposures $ 50 hypothesis $ 4;
      format h df 8.0 Q 8.4 pvalue pvalue8.4;
      %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) = 0 %then %do;
         drop exposures;
      %end;
   run;
%end;

data _&dsnout._deptest;
   length type $ 10 transform trt1 trt2 strata covariates outcomes exposures $ 50 hypothesis $ 4;
   format h df 8.0 beta sebeta Q_j ratio 8.4 pvalue pvalue8.4;
   %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) = 0 %then %do;
      drop exposures;
   %end;
   %if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) = 0 %then %do; 
      drop ratio;
   %end;
run;

%if %upcase(&hypoth) = ALT %then %do;
   data _&dsnout._ci;
      length type $ 10 trt1 trt2 strata covariates outcomes exposures $ 50 hypothesis $ 4;
      format h 8.0 beta sebeta lower upper 8.4 alpha 8.2;
      %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) = 0 %then %do;
         drop exposures;
      %end;
   run;
%end;

%if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) > 0 and %upcase(&hypoth) = ALT %then %do;
   data _&dsnout._ratioci;
      length type $ 10 trt1 trt2 strata covariates outcomes exposures $ 50 hypothesis $ 4;
      format h 8.0 ratio ratio_lower ratio_upper 8.4 alpha 8.2;
      %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) = 0 %then %do;
         drop exposures;
      %end;
      run;
%end;
   
%if %index(PODDS, %upcase(&transform)) > 0 %then %do;
   data _&dsnout._homogen;
      length type $ 10 trt1 trt2 strata covariates outcomes $ 50;
      format df 8.0 Q_c 8.4 pvalue pvalue8.4;
   run;
%end;

%if %upcase(&exact) = YES %then %do;
   data _&dsnout._exact;
      length type $ 10 trt1 trt2 strata covariates outcomes $ 50 hypothesis $ 4;
      format twosided one_lower one_upper pvalue8.4 h nreps seed 8.0 bias accel bca_lower bca_upper pct_lower pct_upper alpha_low alpha_hi 8.6 seed 8.0;
      %if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) and &hypoth = ALT %then %do;
         format bca_ratio_lower bca_ratio_upper pct_ratio_lower pct_ratio_upper 8.6;
      %end;
      %if &hypoth = NULL %then %do;
          drop bias accel bca_upper bca_lower alpha_hi alpha_low pct_lower pct_upper;
      %end;
   run;
%end;

****************;
****************;
****************;

   proc iml symsize = &symsize;
      start NPARCOV4(h, r, t, c, hypoth, f, V, n, varnames, combine, transform); 
         wh = CREATEWH(h, c, n);
         if t>0 then X = I(r)//J(t,r,0);
         else if t=0 then X = I(r);

         if transform = 'LOGISTIC' & combine = 'PRETRANSFORM' then run TRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames, 0, 0, 'LOGISTIC'); 
         else if transform = 'PODDS' & combine = 'PRETRANSFORM' then run PROTRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames);
         else if transform = 'LOGRATIO' & combine = 'PRETRANSFORM' then run TRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames, 0, 0, 'LOGRATIO');
         else if transform = 'INCDENS' & combine = 'PRETRANSFORM' then run TRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames, 0, 0, 'INCDENS');
         else do;
            if (transform = 'LOGISTIC' | transform = 'PODDS') & combine ^= 'PRETRANSFORM' then run TRANSFORM(h, r, t, f, V, hypoth, n, 'LOGISTIC');
            else if (transform = 'LOGRATIO') & combine ^= 'PRETRANSFORM' then run TRANSFORM(h, r, t, f, V, hypoth, n, 'LOGRATIO');
            else if (transform = 'INCDENS') & combine ^= 'PRETRANSFORM' then run TRANSFORM(h, r, t, f, V, hypoth, n, 'INCDENS');

            Vdh = VARSTRAT(h, r, t, V, n);
            dh = CREATEDH(h, r, t, f);

            if transform = 'LOGRANK' & (combine = 'LAST' | combine = 'NONE') then run LAST(n, h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 0); 
            else if transform = 'LOGRANK' & combine = 'FIRST' then run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 0); 
            else if transform = 'WILCOXON' & (combine = 'LAST' | combine = 'NONE') then run LAST(n, h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 0); 
            else if transform = 'WILCOXON' & combine = 'FIRST' then run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 0);
            else if transform = 'LOGISTIC' & (combine = 'LAST' | combine = 'NONE') then run LAST(n, h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 1);
            else if transform = 'LOGISTIC' & combine = 'FIRST' then run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 1);   
            else if transform = 'LOGRATIO' & (combine = 'LAST' | combine = 'NONE') then run LAST(n, h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 1);
            else if transform = 'LOGRATIO' & combine = 'FIRST' then run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 1);
            else if transform = 'PODDS' & combine = 'FIRST' then run PROFIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames);
            else if transform = 'PODDS' & (combine = 'LAST' | combine = 'NONE') then run PROLAST(h, r, t, dh, wh, Vdh, X, hypoth, varnames);
            else if transform = 'NONE' & (combine = 'LAST' | combine = 'NONE') then run LAST(n, h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 0);
            else if transform = 'NONE' & combine = 'FIRST' then run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 0); 
            else if transform = 'INCDENS' & (combine = 'LAST' | combine = 'NONE') then run LAST(n, h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 1);
            else if transform = 'INCDENS' & combine = 'FIRST' then run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 0, 0, 1);
         end;
      finish NPARCOV4;

****************;
****************;
****************;

      start VARSTRAT(h, r, t, V, n);
         VAR = J(h*(r+t), r+t, 0);
         count = 1;
         do i = 1 to h;
            VAR[1+(i-1)*(r+t):i*(r+t),] = V[1+(count-1)*(r+t):count*(r+t),] / n[count] + V[1+(count)*(r+t):(count+1)*(r+t),] / n[count+1]; 
            count = count + 2;
         end;
         return(VAR);
      finish VARSTRAT;

****************;
****************;
****************;

      start CREATEDH(h, r, t, f); 
         D = J(h*(r+t), 1, 0);
         count = 1;
         do i = 1 to h;
            D[1+(i-1)*(r+t):i*(r+t)] =  - f[1+(count-1)*(r+t):count*(r+t)] + f[1+(count)*(r+t):(count+1)*(r+t)]; 
            count = count + 2;
         end;
         return(D);
      finish CREATEDH;

****************;
****************;
****************;

      start CREATEWH(h, c, n);
         wh = J(h, 1, 0);
         count = 1;
         do i = 1 to h;
            wh[i] = ( (n[count] * n[count+1])/(n[count] + n[count+1]) ) ## c;   
            count = count + 2;
         end;
         return(wh);
      finish CREATEWH;

****************;
****************;
****************;

      start LAST(nval, h, r, t, dh, wh, Vdh, X, hypoth, varnames, homoyes, proptitl, oddsrats);
         betaall = J(h*NCOL(X), 1, 0);
         vbetaall = J(h*NCOL(X), NCOL(X), 0);
         if t > 0 then Q = 0;

         n_h = J(h, 1, 0);
         count = 1;
         do i = 1 to h;
            VEE = Vdh[1+(i-1)*(r+t):i*(r+t),];
            dee = dh[1+(i-1)*(r+t):i*(r+t)];
            betaall[1+(i-1)*NCOL(X):i*NCOL(X)] = ESTIMATE(X, VEE, dee); 
            vbetaall[1+(i-1)*NCOL(X):i*NCOL(X),] = VAREST(X, VEE); 

            n_h[i] = nval[count] + nval[count+1];   
            count = count + 2;
     
            if t > 0 then do;
               run COVTESTS(dee, VEE, X, betaall[1+(i-1)*NCOL(X):i*NCOL(X)], t, varnames, Qpart, 0, proptitl, h);
               Q = Q + Qpart;
            end;
         end;

         if h > 1 then do;
            strata_h = {&stratalist};
            beta_h = betaall;         
            create _&dsnout._stratabeta var{strata_h,beta_h,n_h};            
               append;
            close _&dsnout._stratabeta;
            free beta_h strata_h;
         end;
                           
         beta = WGTSUM(h, betaall, wh, 1);
         varbeta = WGTSUM(h, vbetaall, wh, 2);

         create _&dsnout._covbeta from varbeta[colname = varnames];         
            append from varbeta;
         close _&dsnout._covbeta;

         if t > 0 & homoyes = 0 then do;
            df = (NROW(X) - NCOL(X))*h; 
            pvalue = 1 - PROBCHI(Q, df);
            edit _&dsnout._covtest;
            type = "COVTEST";
            trt1 = "&trt1";
            trt2 = "&trt2";
            %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) %then %do;
               exposures = "&exposures";
            %end;
            outcomes = "&workoutcomes";
            covariates = "&covars";
            strata = "&strata";
            hypothesis = "&hypoth";
            append;
            close _&dsnout._covtest;
         end;

         if homoyes = 0 then run DEPTESTS(beta, varbeta, varnames, ncol(X));
         if hypoth = 'ALT' & homoyes = 0 then run CI(beta, varbeta, varnames, oddsrats, ncol(X));
         if homoyes = 1 then run HOMOGEN(beta, varbeta);
      finish LAST;

****************;
****************;
****************;

      start ESTIMATE(X, V, d);
         return(INV(X`*INV(V)*X)*X`*INV(V)*d);
      finish ESTIMATE;

****************;
****************;
****************;

      start VAREST(X, V);
         return(INV(X`*INV(V)*X));
      finish VAREST;

****************;
****************;
****************;

      start COVTESTS(d, Vd, X, beta, t, varnames, Q, printyes, proptitl, h);
         Q = (d - X*beta)` * INV(Vd) * (d - X*beta);
         df = NROW(X)-NCOL(X);
         pvalue = 1 - PROBCHI(Q, df);
  
         if printyes = 1 then do; 
            edit _&dsnout._covtest;
            type = "COVTEST";
            trt1 = "&trt1";
            trt2 = "&trt2";
            covariates = "&covars";
            %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) %then %do;
               exposures = "&exposures";
            %end;
            outcomes = "&workoutcomes";
            strata = "&strata";
            hypothesis = "&hypoth";
            append;
            close _&dsnout._covtest;
         end;
      finish COVTESTS;

****************;
****************;
****************;

      start DEPTESTS(beta, varbeta, varnames, num);
         sebeta = sqrt(VECDIAG(varbeta));
         Q_j = beta##2/VECDIAG(varbeta);
         df = J(NROW(beta), 1, 1);
         pvalue = 1 - PROBCHI(Q_j, df);
         outcomes = varnames[1:num]`;
         %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) %then %do;  
            exposures = {&exposures}`;  
         %end;
         %if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) %then %do; 
            ratio = exp(beta); 
         %end;
         edit _&dsnout._deptest;
         transform = J(NROW(beta), 1, "&transform");
         type = J(NROW(beta), 1, "DEPTEST");
         trt1 = J(NROW(beta), 1, "&trt1");
         trt2 = J(NROW(beta), 1, "&trt2");
         covariates = J(NROW(beta), 1, "&covars");
         h = J(NROW(beta), 1, &numstrat);
         strata = J(NROW(beta), 1, "&strata");
         label = J(NROW(beta), 1, "Test Between Groups for Outcomes");
         hypothesis = J(NROW(beta), 1, "&hypoth");
         append;
         close _&dsnout._deptest;    
      finish DEPTESTS;

****************;
****************;
****************;

      start CI(beta, varbeta, varnames, oddsrats, num);
         sebeta = sqrt(VECDIAG(varbeta));
         lower = beta - probit(%sysevalf(1-&alpha/2))*sebeta;
         upper = beta + probit(%sysevalf(1-&alpha/2))*sebeta;

         edit _&dsnout._ci;
         outcomes = varnames[1:num]`;
         %if %index(INCDENS LOGRANK WILCOXON, %upcase(&transform)) %then %do;  
            exposures = {&exposures}`;  
         %end;
         type = J(NROW(beta), 1, "CI");
         trt1 = J(NROW(beta), 1, "&trt1");
         trt2 = J(NROW(beta), 1, "&trt2");
         covariates = J(NROW(beta), 1, "&covars");
         h = J(NROW(beta), 1, &numstrat);
         alpha = J(NROW(beta), 1, &alpha);
         strata = J(NROW(beta), 1, "&strata");
         hypothesis = J(NROW(beta), 1, "&hypoth");
         append;
         close _&dsnout._ci;

         if oddsrats=1 then do;
            ratio = EXP(beta);
            ratio_lower = EXP(lower);
            ratio_upper = EXP(upper);

            edit _&dsnout._ratioci;
            type = J(NROW(beta), 1, "RATIOCI");      
            trt1 = J(NROW(beta), 1, "&trt1");
            trt2 = J(NROW(beta), 1, "&trt2");
            covariates = J(NROW(beta), 1, "&covars");
            h = J(NROW(beta), 1, &numstrat);
            alpha = J(NROW(beta), 1, &alpha);
            strata = J(NROW(beta), 1, "&strata");
            hypothesis = J(NROW(beta), 1, "&hypoth");
            append;
            close _&dsnout._ratioci;
         end;
      finish CI;

****************;
****************;
****************;

      start WGTSUM(h, mats, wh, expo);
         total = J(NROW(mats)/h, NCOL(mats), 0);
         do i = 1 to h;
            total = total + (wh[i] ## expo) * mats[1+(i-1)*(NROW(mats)/h):i*(NROW(mats)/h),];
         end;
         total = total / (sum(wh) ## expo);
         return(total); 
      finish WGTSUM;

****************;
****************;
****************;

      start WGTSUM2(h, r, t, mats, wh, expo);
         total = J(2*(r+t), NCOL(mats), 0);
         count = 1;
         do i = 1 to h;
            total[1:(r+t),] = total[1:(r+t),] + (wh[i] ## expo) * mats[1+(count-1)*(r+t):count*(r+t),];
            total[(r+t+1):2*(r+t),] = total[(r+t+1):2*(r+t),] + (wh[i] ## expo) * mats[1+count*(r+t):(count+1)*(r+t),];
            count = count + 2;
         end;
         total = total / (sum(wh) ## expo);
         return(total); 
      finish WGTSUM2;

****************;
****************;
****************;

      start FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, homoyes, proptitl, oddsrats);
         d = WGTSUM(h, dh, wh, 1); 
         Vd = WGTSUM(h, Vdh, wh, 2);         
         beta = ESTIMATE(X, Vd, d); 
         varbeta = VAREST(X, Vd); 
         create _&dsnout._covbeta from varbeta[colname = varnames];            
            append from varbeta;
         close _&dsnout._covbeta;
         if t > 0 & homoyes = 0 then run COVTESTS(d, Vd, X, beta, t, varnames, Q, 1, proptitl, h);
         if homoyes = 0 then run DEPTESTS(beta, varbeta, varnames, ncol(X));
         if hypoth = 'ALT' & homoyes = 0 then run CI(beta, varbeta, varnames, oddsrats, ncol(X));
         if homoyes = 1 then run HOMOGEN(beta, varbeta);
      finish FIRST;

****************;
****************;
****************;

      start PROFIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames);
         run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 1, 1, 1);
         X = J(r+t, 1, 0);
         X[1:r] = 1;
         newname = J(t+1, 1, 'TREATMENT');
         if t > 0 then newname[2:t+1] = varnames[r+1:r+t];
         run FIRST(h, r, t, dh, wh, Vdh, X, hypoth, newname, 0, 1, 1);
      finish PROFIRST;

****************;
****************;
****************;

      start PROLAST(h, r, t, dh, wh, Vdh, X, hypoth, varnames);
         run LAST(h, r, t, dh, wh, Vdh, X, hypoth, varnames, 1, 0, 1);
         X = J(r+t, 1, 0);
         X[1:r] = 1;
         newname = J(t+1, 1, 'TREATMENT');
         if t > 0 then newname[2:t+1] = varnames[r+1:r+t];
         run LAST(h, r, t, dh, wh, Vdh, X, hypoth, newname, 0, 1, 1);
      finish PROLAST;

****************;
****************;
****************;

      start TRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames, homoyes, proptitl, transtype); 
         %if %index(LOGISTIC LOGRATIO, %upcase(&transform)) > 0 %then %do;
            VEE = J(2*h*(r+t), r+t, 0);

            do i = 1 to 2*h;
               VEE[1+(i-1)*(r+t):i*(r+t),] = V[1+(i-1)*(r+t):i*(r+t),] / n[i];
            end;

            fstar = WGTSUM2(h, r, t, f, wh, 1);
            Vstar = WGTSUM2(h, r, t, VEE, wh, 2);            
         %end;
         %else %if %index(INCDENS, %upcase(&transform)) %then %do;
            VEE = J(2*h*(2*r+t), 2*r+t, 0);

            do i = 1 to 2*h;
               VEE[1+(i-1)*(2*r+t):i*(2*r+t),] = V[1+(i-1)*(2*r+t):i*(2*r+t),] / n[i];
            end;

            fstar = WGTSUM2(h, 2*r, t, f, wh, 1);
            Vstar = WGTSUM2(h, 2*r, t, VEE, wh, 2);
         %end;
         run TRANSFORM(1, r, t, fstar, Vstar, hypoth, n, transtype);
         Vdh = VARSTRAT(1, r, t, Vstar, J(2,1,1));
         dh = CREATEDH(1, r, t, fstar);
         
         run FIRST(1, r, t, dh, 1, Vdh, X, hypoth, varnames, homoyes, proptitl, 1);
      finish TRANSFIRST;

****************;
****************;
****************;

      start PROTRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames);
         run TRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, varnames, 1, 0, 'LOGISTIC');
         X = J(r+t, 1, 0);
         X[1:r] = 1;
         newname = J(t+1, 1, 'TREATMENT');
         if t > 0 then newname[2:t+1] = varnames[r+1:r+t];
         run TRANSFIRST(h, r, t, n, f, wh, V, X, hypoth, newname, 0, 1, 'LOGISTIC');
      finish PROTRANSFIRST;

****************;
****************;
****************;

      start HOMOGEN(beta, varbeta);
         C = I(NROW(beta)-1)||J(NROW(beta)-1, 1, -1);
         Q_c = beta` * C` * INV(C * varbeta * C`) * C * beta;
         df = NROW(beta) - 1;
         pvalue = 1 - PROBCHI(Q_c, df);
         outcomes = "TREATMENT";

         edit _&dsnout._homogen;
            type = "HOMOGEN";
            trt1 = "&trt1";
            trt2 = "&trt2";
            covariates = "&covars";
            strata = "&strata";
            hypoth = "&hypoth";
            append;
         close _&dsnout._homogen;
      finish HOMOGEN;

****************;
****************;
****************;

      start TRANSFORM(h, r, t, f, V, hypoth, n, transtype);
         if hypoth = 'ALT' then run TRANSVARALT(h, r, t, f, V, transtype);
         else if hypoth = 'NULL' then run TRANSVARNULL(h, r, t, f, V, n, transtype);
         run TRANSF(h, r, t, f, transtype);
      finish TRANSFORM;

****************;
****************;
****************;

      start TRANSVARALT(h, r, t, f, V, transtype);
         if transtype = "LOGISTIC" | transtype = "LOGRATIO" then do;
            d = J(r+t, 1, 0);
            if t > 0 then d[r+1:r+t] = 1;
            do i = 1 to 2*h;
               temp = f[1+(i-1)*(r+t):(i*r)+(i-1)*t];
               if transtype = "LOGISTIC" then d[1:r] = 1/(f[1+(i-1)*(r+t):(i*r)+(i-1)*t]#(1-f[1+(i-1)*(r+t):(i*r)+(i-1)*t]));
               else if transtype = "LOGRATIO" then d[1:r] = 1/f[1+(i-1)*(r+t):(i*r)+(i-1)*t];
               V[1+(i-1)*(r+t):i*(r+t),] = (V[1+(i-1)*(r+t):i*(r+t),] # d) # d`;
            end;
         end;
         else if transtype = "INCDENS" then do;
            Vtemp = J(2*h*(r+t), (r+t), 0); 
            do i = 1 to 2*h;
               dy = diag(1/f[1+(i-1)*(2*r+t):(i*r)+(i-1)*(r+t)]);
               dt = -diag(1/f[(r+1)+(i-1)*(2*r+t):(2*i*r)+(i-1)*t]);
               if t = 0 then M = dy||dt;
               else do;
                  top = dy||dt||J(r,t,0); 
                  bottom = J(t,2*r,0)||I(t);
                  M = top // bottom;
                  free top bottom;
               end;
               Vtemp[1+(i-1)*(r+t):i*(r+t),] = M * V[1+(i-1)*(2*r+t):i*(2*r+t),] * M`;
            end;
            V = Vtemp;
            free Vtemp M;
         end;
      finish TRANSVARALT;

****************;
****************;
****************;

      start TRANSVARNULL(h, r, t, f, V, n, transtype);
         if transtype = "LOGISTIC" | transtype = "LOGRATIO" then do;
            d = J(r+t, 1, 0);
            if t > 0 then d[r+1:r+t] = 1;
            count = 1;
            do i = 1 to h;
               ybar = (n[count] * f[1+(count-1)*(r+t):(count*r)+(count-1)*t] + n[count+1] * f[1+(count)*(r+t):((count+1)*r)+(count)*t]) /
                      (n[count] + n[count+1]);
               if transtype = "LOGISTIC" then d[1:r] = 1/(ybar # (1-ybar));
               else if transtype = "LOGRATIO" then d[1:r] = 1/ybar;
               V[1+(count-1)*(r+t):count*(r+t),] = (V[1+(count-1)*(r+t):count*(r+t),] # d) # d`;
               V[1+(count)*(r+t):(count+1)*(r+t),] = (V[1+(count)*(r+t):(count+1)*(r+t),] # d) # d`;
               count = count + 2;
            end;
         end;
         else if transtype = "INCDENS" then do;
            Vtemp = J(2*h*(r+t), (r+t), 0); 
            count = 1;
            do i = 1 to h;
               ytbar = (n[count] * f[1+(count-1)*(2*r+t):(count*2*r)+(count-1)*t] + n[count+1] * f[1+(count)*(2*r+t):((count+1)*2*r)+(count)*t]) /
                       (n[count] + n[count+1]);
               dy = diag(1/ytbar[1:r]);
               dt = -diag(1/ytbar[r+1:2*r]);
               if t = 0 then M = dy||dt;
               else do;
                  top = dy||dt||J(r,t,0); 
                  bottom = J(t,2*r,0)||I(t);
                  M = top // bottom;
                  free top bottom;
               end;
               Vtemp[1+(count-1)*(r+t):count*(r+t),] = M * V[1+(count-1)*(2*r+t):count*(2*r+t),] * M`;
               Vtemp[1+(count)*(r+t):(count+1)*(r+t),] = M * V[1+(count)*(2*r+t):(count+1)*(2*r+t),] * M`;
               count = count + 2;
            end;
            V = Vtemp;
            free Vtemp M;
         end;
      finish TRANSVARNULL;

****************;
****************;
****************;

      start TRANSF(h, r, t, f, transtype);
         do i = 1 to 2*h;
            if transtype = "LOGISTIC" then f[1+(i-1)*(r+t):(i*r)+(i-1)*t] = LOG(f[1+(i-1)*(r+t):(i*r)+(i-1)*t] / (1 - f[1+(i-1)*(r+t):(i*r)+(i-1)*t]));
            else if transtype = "LOGRATIO" then f[1+(i-1)*(r+t):(i*r)+(i-1)*t] = LOG(f[1+(i-1)*(r+t):(i*r)+(i-1)*t]);
            else if transtype = "INCDENS" then f[1+(i-1)*(2*r+t):(i*2*r)+(i-1)*t] = LOG(f[1+(i-1)*(2*r+t):(i*2*r)+(i-1)*t]);
         end;
         if transtype = "INCDENS" then do;
            ftemp = J(2*h*(r+t), 1, 0);
            if t = 0 then M = I(r)||-I(r);
            else do;
               top = I(r)||-I(r)||J(r,t,0); 
               bottom = J(t,2*r,0)||I(t);
               M = top // bottom;
               free top bottom;
            end;
            do i = 1 to 2*h;
               ftemp[1+(i-1)*(r+t):i*(r+t)] = M * f[1+(i-1)*(2*r+t):i*(2*r+t)];
            end;
            f = ftemp;
            free ftemp M;
         end;
      finish TRANSF;
      
****************;
****************;
****************;

      start RANDPARCOV(h, r, t, c, hypoth, f, V, n, varnames, transform, seed, nreps, combine); 
         if t>0 then X = I(r)//J(t,r,0);
         else if t=0 then X = I(r);

         if hypoth = "ALT" then loopcount = nreps + &numobstot + 1;
         else loopcount = nreps + 1;

         twosided = 0;
         one_lower = 0;
         one_upper = 0;
         
         betasamp = J(loopcount, 1, 0);
         if t > 0 then covxact = 0;
         
         if combine = 'FIRST' then do;          
            do i = 1 to loopcount;     
               if (transform = 'LOGISTIC' | transform = 'LOGRATIO') then do;
                  fsub = f[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];
                  Vsub = V[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];               
                  run TRANSFORM(h, r, t, fsub, Vsub, hypoth, n[1+(i-1)*2*h:i*2*h], transform);
               end;
               else if transform = 'INCDENS' then do;
                  fsub = f[1+(i-1)*2*h*(2*r+t):i*2*h*(2*r+t),];
                  Vsub = V[1+(i-1)*2*h*(2*r+t):i*2*h*(2*r+t),];
                  run TRANSFORM(h, r, t, fsub, Vsub, hypoth, n[1+(i-1)*2*h:i*2*h], transform);               
               end;
               else do;
                  fsub = f[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];
                  Vsub = V[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];                              
               end;

               dh = CREATEDH(h, r, t, fsub);
               Vdh = VARSTRAT(h, r, t, Vsub, n[1+(i-1)*2*h:i*2*h]);               
               wh = CREATEWH(h, c, n[1+(i-1)*2*h:i*2*h]);                 
               Vd = WGTSUM(h, Vdh, wh, 2); 
               d = WGTSUM(h, dh, wh, 1);
               
               betasamp[i] = ESTIMATE(X, Vd, d);
               if i = 1 & t > 0 then covxact0 = (d - X*betasamp[1+(i-1)*r:i*r])`*INV(Vd)*(d - X*betasamp[1+(i-1)*r:i*r]);
               if (i >= 2) & (i <= (nreps + 1)) then do;
                  twosided = twosided + (abs(betasamp[i]) >= abs(betasamp[1]));
                  one_lower = one_lower + (betasamp[i] <= betasamp[1]);
                  one_upper = one_upper + (betasamp[i] >= betasamp[1]);
                  if t > 0 then covxact = covxact + ((d - X*betasamp[i])`*INV(Vd)*(d - X*betasamp[i]) >= covxact0);
               end; 
            end;
         end;
         else if (combine = 'LAST' | combine = 'NONE') then do;
            if hypoth = "ALT" then do;
               jackmean = J(h, 1, 0);
               jackest = J(&numobstot,1,0);
               l = 1;
               ntemp = sum(n[1:2]);
               ntot = ntemp;
            end;

            do i = 1 to loopcount;
               betaall = J(h, 1, 0);
               vbetaall = J(h, r, 0);
               if t > 0 then Qxact = 0;

               if (transform = 'LOGISTIC' | transform = 'LOGRATIO') then do;
                  fsub = f[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];
                  Vsub = V[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];
                  run TRANSFORM(h, r, t, fsub, Vsub, hypoth, n[1+(i-1)*2*h:i*2*h], transform);
               end;
               else if transform = 'INCDENS' then do;                   
                  fsub = f[1+(i-1)*2*h*(2*r+t):i*2*h*(2*r+t),];
                  Vsub = V[1+(i-1)*2*h*(2*r+t):i*2*h*(2*r+t),];
                  run TRANSFORM(h, r, t, fsub, Vsub, hypoth, n[1+(i-1)*2*h:i*2*h], transform);
               end;
               else do;
                  fsub = f[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];
                  Vsub = V[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];                              
               end;
               
               dh = CREATEDH(h, r, t, fsub);
               Vdh = VARSTRAT(h, r, t, Vsub, n[1+(i-1)*2*h:i*2*h]);
               
               wh = CREATEWH(h, c, n[1+(i-1)*2*h:i*2*h]);

               do j = 1 to h;
                  VEE = Vdh[1+(j-1)*(r+t):j*(r+t),];
                  dee = dh[1+(j-1)*(r+t):j*(r+t)];
                  betaall[1+(j-1)*r:j*r] = ESTIMATE(X, VEE, dee); 
                  vbetaall[1+(j-1)*r:j*r,] = VAREST(X, VEE);      

                  if t > 0 then Qxact = Qxact + (dee - X*betaall[1+(j-1)*r:j*r])`*INV(VEE)*(dee - X*betaall[1+(j-1)*r:j*r]);

                  if hypoth = "ALT" then do;
                     if (j = l) & (i >= (nreps + 2)) then do;
                        jackest[i-nreps-1] = ESTIMATE(X, VEE, dee);
                        jackmean[l] = jackmean[l] + jackest[i-nreps-1];
                     end;
                  end;
               end;

               if hypoth = "ALT" then do;
                  if i = (nreps + 1 + ntot) then do;
                     jackmean[l] = jackmean[l] / ntemp;
                     l = l + 1;
                     ntemp = sum(n[1+(l-1)*2:2*l]);
                     ntot = ntot + ntemp;
                  end; 
               end;

               betasamp[i] = WGTSUM(h, betaall, wh, 1);
               if i = 1 & t > 0 then covxact0 = Qxact;

               if (i >= 2) & (i <= (nreps + 1)) then do;
                  twosided = twosided + (abs(betasamp[i]) >= abs(betasamp[1]));
                  one_lower = one_lower + (betasamp[i] <= betasamp[1]);
                  one_upper = one_upper + (betasamp[i] >= betasamp[1]);
                  if t > 0 then covxact = covxact + (Qxact >= covxact0);
               end; 
            end;                           
         end;
         else if (transform = 'LOGISTIC' | transform = 'LOGRATIO'| transform = 'INCDENS') & combine = 'PRETRANSFORM' then do;         
            if (transform = 'LOGISTIC' | transform = 'LOGRATIO') then do;            
               VEE = J(loopcount*2*h*(r+t), r+t, 0);
            end;
            else if transform = 'INCDENS' then do;
               VEE = J(loopcount*2*h*(2*r+t), 2*r+t, 0);
            end;
            l = 1;
            do i = 1 to loopcount;
               do j = 1 to 2*h;
                  ntemp = n[1+(i-1)*2*h:i*2*h];
                  if (transform = 'LOGISTIC' | transform = 'LOGRATIO') then do;
                     VEE[1+(l+j-2)*(r+t):(l+j-1)*(r+t),] = V[1+(l+j-2)*(r+t):(l+j-1)*(r+t),] / ntemp[j];
                  end;
                  else if transform = 'INCDENS' then do;
                     VEE[1+(l+j-2)*(2*r+t):(l+j-1)*(2*r+t),] = V[1+(l+j-2)*(2*r+t):(l+j-1)*(2*r+t),] / ntemp[j];
                  end;
               end;     
               l = l + (2 * h);
            end;            
           
            do i = 1 to loopcount; 
               wh = CREATEWH(h, c, n[1+(i-1)*2*h:i*2*h]);                             
               if (transform = 'LOGISTIC' | transform = 'LOGRATIO') then do;       
                  fsub = f[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];
                  Vsub = VEE[1+(i-1)*2*h*(r+t):i*2*h*(r+t),];                 
                  
                  fstar = WGTSUM2(h, r, t, fsub, wh, 1);           
                  Vstar = WGTSUM2(h, r, t, Vsub, wh, 2);                  
               end;
               else if transform = 'INCDENS' then do;
                  fsub = f[1+(i-1)*2*h*(2*r+t):i*2*h*(2*r+t),];
                  Vsub = VEE[1+(i-1)*2*h*(2*r+t):i*2*h*(2*r+t),];
                  fstar = WGTSUM2(h, 2*r, t, fsub, wh, 1);           
                  Vstar = WGTSUM2(h, 2*r, t, Vsub, wh, 2);               
               end;             
                              
               run TRANSFORM(1, r, t, fstar, Vstar, hypoth, n, transform); 
               Vdh = VARSTRAT(1, r, t, Vstar, J(2,1,1));
               dh = CREATEDH(1, r, t, fstar);             
               d = WGTSUM(1, dh, 1, 1);
               Vd = WGTSUM(1, Vdh, 1, 2);                               
               betasamp[i] = ESTIMATE(X, Vd, d);               
               
               if i = 1 & t > 0 then covxact0 = (d - X*betasamp[i])`*INV(Vd)*(d - X*betasamp[i]);
               if (i >= 2) & (i <= (nreps + 1)) then do;
                  twosided = twosided + (abs(betasamp[i]) >= abs(betasamp[1]));
                  one_lower = one_lower + (betasamp[i] <= betasamp[1]);
                  one_upper = one_upper + (betasamp[i] >= betasamp[1]);
                  if t > 0 then covxact = covxact + ((d - X*betasamp[i])`*INV(Vd)*(d - X*betasamp[i]) >= covxact0);
               end; 
            end;
         end;

         twosided = twosided / nreps;
         one_lower = one_lower / nreps;
         one_upper = one_upper / nreps;
         if t > 0 then covxact = covxact / nreps;

         create _&dsnout._betasamp var{betasamp};
            append;
         close _&dsnout._betasamp;

         if t > 0 then do;
            create _&dsnout._covxact var{covxact};
            append;
            close _&dsnout._covxact;
         end;

         if hypoth = "ALT" & (combine = "FIRST" | combine = "PRETRANSFORM") then do;
            bias = probit(sum((betasamp[2:(nreps+1)] < betasamp[1]) ) / nrow(betasamp[2:(nreps+1)]));
            jackmean = sum(betasamp[(nreps+2):nreps+&numobstot+1])/&numobstot;  
            accel = sum((betasamp[(nreps+2):nreps+&numobstot+1]-jackmean)##3) / (6*(sum((betasamp[(nreps+2):nreps+&numobstot+1]-jackmean)##2))**(1.5));
            alpha_low = probnorm(bias + (bias + probit(&alpha/2))/(1-accel*(bias + probit(&alpha/2))));
            alpha_hi = probnorm(bias + (bias + probit(1-&alpha/2))/(1-accel*(bias + probit(1-&alpha/2))));
         end;
         else if hypoth = "ALT" & (combine = 'LAST' | combine = 'NONE') then do;
            jackdelta = J(&numobstot,1,0);
            l = 1;
            ntemp = sum(n[1:2]);
            ntot = ntemp;

            do i = (nreps + 2) to (nreps + &numobstot + 1);
               jackdelta[i-nreps-1] = jackest[i-nreps-1] - jackmean[1+(l-1)*r:l*r];

               if i = (nreps + 1 + ntot) then do;
                  l = l + 1;
                  ntemp = sum(n[1+(l-1)*2:2*l]);
                  ntot = ntot + ntemp;
               end; 
            end;

            l = 1;
            ntemp = sum(n[1:2]);
            ntot = ntemp;
            jacksum2 = J(r*h,1,0);
            jacksum3 = J(r*h,1,0);
            do i = (nreps + 2) to (nreps + &numobstot + 1);
               jacksum2[l] = jacksum2[l] + jackdelta[i-nreps-1]##2;
               jacksum3[l] = jacksum3[l] + jackdelta[i-nreps-1]##3;

               if i = (nreps + 1 + ntot) then do;
                  l = l + 1;
                  ntemp = sum(n[1+(l-1)*2:2*l]);
                  ntot = ntot + ntemp;
               end; 
            end;

            jacktotsum2 = J(r,1,0);
            jacktotsum3 = J(r,1,0);

            do l = 1 to h;
               ntemp = sum(n[1+(l-1)*2:2*l]);
               jacktotsum2 = jacktotsum2 + jacksum2[l] / ntemp##2;
               jacktotsum3 = jacktotsum3 + jacksum3[l] / ntemp##3;
            end;

            bias = probit(sum((betasamp[2:(nreps+1)] < betasamp[1]) ) / nrow(betasamp[2:(nreps+1)]));
            accel = jacktotsum3 / (6*jacktotsum2##(3/2));      
            alpha_low = probnorm(bias + (bias + probit(&alpha/2))/(1-accel*(bias + probit(&alpha/2))));
            alpha_hi = probnorm(bias + (bias + probit(1-&alpha/2))/(1-accel*(bias + probit(1-&alpha/2))));
         end; 

         type = "EXACT";
         trt1 = "&trt1";
         trt2 = "&trt2";
         outcomes = "&outcomes";
         covariates = "&covars";
         strata = "&strata";
         hypothesis = "&hypoth";
         seed = &seed;
         edit _&dsnout._exact;
            append;
         close _&dsnout._exact;
      finish RANDPARCOV;

****************;
****************;
****************;

      use _work_means where(_sample_ = 0);
         read all var{&workoutcomes %if %upcase(&transform) = INCDENS %then %do; &exposures %end; &covars} into f;
         read all var{_freq_} into n;
      close;
      use _work_covar where(_sample_ = 0);
         read all var{&workoutcomes %if %upcase(&transform) = INCDENS %then %do; &exposures %end; &covars} into V;
      close;

      f = COLVEC(f);
      varnames = {&workoutcomes &covars}`;
      strat = "%upcase(&strata)";
      combine = "%upcase(&combine)";
      transform = "%upcase(&transform)";
      hypoth = "%upcase(&hypoth)";

      run NPARCOV4(&numstrat, &numresp, &numcov, &c, hypoth, f, V, n, varnames, combine, transform);
            
      %if %upcase(&exact) = YES %then %do;
         free f n V;
         use _work_means;
            read all var{&workoutcomes %if %upcase(&transform) = INCDENS %then %do; &exposures %end; &covars} into f;
            read all var{_freq_} into n;
         close;
         use _work_covar;
            read all var{&workoutcomes %if %upcase(&transform) = INCDENS %then %do; &exposures %end; &covars} into V;
         close;

         f = COLVEC(f);
         run RANDPARCOV(&numstrat, &numresp, &numcov, &c, hypoth, f, V, n, varnames, transform, &seed, &nreps, combine);
      %end;
      quit;
      
      %if &numcov >= 1 %then %do;
         data _&dsnout._covtest(where = (type ne ""));
            set _&dsnout._covtest;
         run;
      %end;
      
      %if %upcase(&hypoth) = ALT %then %do;
         data _&dsnout._ci(where = (type ne ""));
            set _&dsnout._ci;
         run;
      %end;
      
      %if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) > 0 and %upcase(&hypoth) = ALT %then %do;
         data _&dsnout._ratioci(where = (type ne ""));
            set _&dsnout._ratioci;
         run;      
      %end;
      
      %if %index(PODDS, %upcase(&transform)) > 0 %then %do;
         data _&dsnout._homogen(where = (type ne ""));
            set _&dsnout._homogen;
         run;            
      %end;

      %if %upcase(&exact) = YES %then %do;
         data _&dsnout._exact(where = (type ne ""));
            set _&dsnout._exact;
         run;            
      %end;      
           
      %if %upcase(&exact) = YES %then %do;
         data _&dsnout._betasamp;
            length flag $ 12;
            set _&dsnout._betasamp;
            if _n_ = 1 then flag = "OBSERVED";
            %if %upcase(&hypoth) = ALT %then %do; 
               else if 2 <= _n_ <= %sysevalf(&nreps + 1) then flag = "BOOTSTRAP";
            %end;
            %else %if %upcase(&hypoth) = NULL %then %do; 
               else if 2 <= _n_ <= %sysevalf(&nreps + 1) then flag = "PERMUTATION";
            %end;            
            %if %upcase(&hypoth) = ALT %then %do; 
               else if _n_ > %sysevalf(&nreps + 1) then flag = "JACKKNIFE";
            %end;
         run;
         
         %if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) %then %do;
            data _&dsnout._betasamp;
               set _&dsnout._betasamp;
               EXPBETASAMP = exp(betasamp);
            run;         
         %end;
      %end;

      %if %upcase(&exact) = YES and &numcov >= 1 %then %do;
         data _&dsnout._covtest;
            merge _&dsnout._covtest _&dsnout._covxact;
         run;

         proc datasets nolist;
            delete _&dsnout._covxact;
         quit;
      %end;

      %if %upcase(&exact) = YES and %upcase(&hypoth) = ALT %then %do;
         data _&dsnout._exact;
            set _&dsnout._exact;
            call symput("alpha_low", alpha_low);
            call symput("alpha_hi", alpha_hi);
            drop one_lower one_upper twosided;
         run;

         %PCTCI(flag = BCA, low = %sysevalf(&alpha_low * 100), hi = %sysevalf(&alpha_hi * 100));

         %PCTCI(flag = PCT, low = %sysevalf(&alpha/2 * 100), hi = %sysevalf((1 - &alpha/2) * 100));

         data _&dsnout._exact;
            merge _&dsnout._exact 
                  _&dsnout._pctfinal
                  _&dsnout._bcafinal;
            %if %index(LOGISTIC PODDS LOGRATIO INCDENS, %upcase(&transform)) %then %do;
               bca_ratio_lower = exp(bca_lower);
               bca_ratio_upper = exp(bca_upper);
               pct_ratio_lower = exp(pct_lower);
               pct_ratio_upper = exp(pct_upper);
            %end;                  
         run; 

         data _&dsnout._covtest;
            set _&dsnout._covtest;
            drop covxact;
         run;
         
         proc datasets nolist;
            delete _&dsnout._bcafinal _&dsnout._pctfinal;
         quit;
      %end;      

      data _&dsnout._covbeta;
         length type $ 10;
         set _&dsnout._covbeta;
         type = "COVBETA";
      run;

      proc datasets nolist;
         delete _WORK_MEANS _WORK_COVAR _WORK_TEMP _WORK_NSTRT _WORK_NTRTS _WORK_NTRTS_CHK _WORK_OUTSAMP _WORK_JACK;
      quit;
      
      %if %upcase(&details) = YES %then %do;
         proc datasets;
         quit;
      %end;
      
      data _&dsnout._deptest(where = (type ne ""));
         set _&dsnout._deptest;
      run;
      
      ods results;
   %stopmac:;
%mend NPARCOV4;

****************;
****************;
****************;

%macro PCTCI(flag = BCA, low = %sysevalf(&alpha_low * 100), hi = %sysevalf(&alpha_hi * 100));
   proc univariate noprint data = _&dsnout._betasamp;
      var betasamp;
      output out = _&dsnout._&flag.final pctlpts = &low &hi pctlpre = &flag._; 
      where flag = "BOOTSTRAP";
   run;

   proc transpose data = _&dsnout._&flag.final out = _&dsnout._&flag.final;
   run;

   data _&dsnout._&flag.final;
      set _&dsnout._&flag.final;
      if _n_ = 1 then ci = "&flag._LOWER";
      else ci = "&flag._UPPER";
   run;

   proc transpose data = _&dsnout._&flag.final out = _&dsnout._&flag.final(drop = _name_);
      var col1;
      id ci;
   run;
%mend PCTCI;

****************;
****************;
****************;

%macro MISSDAT(varlist = , data = );
   data _work_miss;
      set &data;
      totmiss = nmiss(of &varlist);
   run;
   
   proc means sum noprint data = _work_miss;
      var totmiss;
      output out = _work_miss(drop = _TYPE_ _FREQ_) sum = sum;
   run;
    
   data _null_;
      set _work_miss;
      call symput('missind', trim(left(put(sum,8.))));
   run;

   proc datasets nolist;
      delete _work_miss;
   quit;
%mend MISSDAT;   

****************;
****************;
****************;

%macro SURVSCORE(eventlist = , timelist = , score = %str(LOGRANK));
   %do sloop = 1 %to &numresp;
      %let event = %scan(&eventlist, &sloop);
      %let time = %scan(&timelist, &sloop);

      proc freq noprint data = _work_temp;
         tables &time * &event / out = _work_surv(drop = percent);
      run;

      proc means n noprint data = _work_temp;
         var &time;
         output out = _work_tot(drop = _freq_ _type_) n = total;
      run;

      proc sort data = _work_surv out = _work_surv;
         by &time descending &event;
      run;

      %if &score = LOGRANK %then %do;
         data _work_surv;
            set _work_surv;
            if _n_ = 1 then set _work_tot;
            ratio = 0;   
            retain survscore cumsum atrisk cumratio;
            if _n_ = 1 then cumsum = count;
            else cumsum = cumsum + count;
            if &event = 1 then do;
               atrisk = total - cumsum + count;
               ratio = count / atrisk;
            end;
            if _n_ = 1 then cumratio = ratio;
            else cumratio = cumratio + ratio;
            if _n_ = 1 and &event = 0 then survscore = 0;
            else do;
               if &event = 1 then survscore = 1 - cumratio;
               else if &event = 0 then survscore = - cumratio;
            end;
            keep &time &event count total survscore;
         run;
      %end;
      %else %if &score = WILCOXON %then %do;      
         data _work_surv;
            set _work_surv;
            if _n_ = 1 then set _work_tot;
            ratio = 1;
            retain survscore cumsum atrisk cumratio;
            if _n_ = 1 then cumsum = count;
            else cumsum = cumsum + count;
            if &event = 1 then do;
               atrisk = total - cumsum + count;
               ratio = (atrisk - count) / atrisk;
            end;
            if _n_ = 1 then cumratio = ratio;
            else cumratio = cumratio * ratio;
            if _n_ = 1 and &event = 0 then survscore = 0;   
            else do;
               if &event = 1 then survscore = 2 * cumratio - 1;
               else if &event = 0 then survscore = cumratio - 1;
            end;
            keep &time &event count total survscore;
         run;
      %end;

      proc sort data = _work_temp;
         by &time descending &event;
      run;

      data _work_temp;
         merge _work_temp _work_surv(keep = &time &event survscore rename = (survscore = &score._&event));
         by &time descending &event;
      run;

      data _&dsnout._surv;
         set _work_temp;
         %if &numstrat = 1 %then %do;
            drop _NONE;
         %end;
      run;
   %end;
%mend SURVSCORE;

************************************************************************************************************;
************************************************************************************************************;
************************************************************************************************************;