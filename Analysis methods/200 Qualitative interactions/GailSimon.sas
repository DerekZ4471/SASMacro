/***************************************************

GailSimon macro


Dmitrienko, A., Koch, G. G., and SAS Institute (eds.) (2017), Analysis of clinical trials using SAS: a practical guide, Cary, NC: SAS Institute.

Macro written by Alan Cantor

***************************************************/

%macro GailSimon(dataset,est,stderr,testtype);
/*
Inputs:

DATASET  = Data set with test statistics and associated standard 
           errors for each stratum.

EST      = Name of the variable containing the test statistics.

STDERR   = Name of the variable containing the standard errors.

TESTTYPE = P, N, T to carry out the one-sided Gail-Simon test
           for positive or negative differences or the two-sided
           Gail-Simon test, respectively.
*/
data pvalue;
    set &dataset nobs=m;
    format stat 6.3 p 6.4;
    retain qminus 0 qplus 0;    
    qminus=qminus+(&est>0)*(&est/&stderr)**2;
    qplus=qplus+(&est<0)*(&est/&stderr)**2;
    if _n_=m then do;
        if upcase(&testtype)='P' then do; stat=qplus; df=m+1; end;
        if upcase(&testtype)='N' then do; stat=qminus; df=m+1; end;
        if upcase(&testtype)='T' then do; stat=min(qminus,qplus); df=m; end;
        p=0;
        do i=1 to df-1;
            p=p+pdf("binomial",i,0.5,df-1)*(1-probchi(stat,i));
        end;
    end;
    label stat='Test statistic' p='P-value';
    if _n_=m;
    keep stat p;
proc print data=pvalue noobs label;
    %if %upcase(&testtype)='P' %then %do; 
        title 'One-sided Gail-Simon test for positive differences'; %end;
    %if %upcase(&testtype)='N' %then %do; 
        title 'One-sided Gail-Simon test for negative differences'; %end;
    %if %upcase(&testtype)='T' %then %do; 
        title 'Two-sided Gail-Simon test'; %end;
    run; 
%mend;
