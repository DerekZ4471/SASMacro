/***************************************************

LinRank macro

***************************************************/
%macro LinRank(dataset=_last_, time=time, event=  ,
               censval=, groupvar=, method=logrank,rho=1,
               stratvar= _none_, stratmis=no,trend=order );

/*       Delete invalid observations and 
         print list of observations deleted. */  

data xx deleted; 
     set &dataset;
     obsnumb = _n_;
     if &time <0  or &groupvar=''   or &event = . then delete=1;
     if "&stratvar" ne "_none_" and  "&stratmis" = "no"
     and  &stratvar = '' then delete =1;
     _none_ = 1;
     if delete=1 then output deleted;
     else output xx;
proc print data=deleted; 
     title 'Deleted Observations';
     var  obsnumb &time &event &groupvar
     %if "&stratvar" ne "_none_" %then &stratvar;;

/*       Determine number of groups, their names, 
         and the weights to use for trend test.   */ 

proc sort data = xx; 
     by &groupvar;
data y; 
     set xx; 
     by &groupvar;
     if first.&groupvar then do;
        n+1;
        call symput('ngrps',left(n));
     end;
run;
data grpnames; 
     set y; 
     by &groupvar;
     keep &groupvar n;
     if first.&groupvar;
data groupwts; 
     set grpnames;
     keep
     %if "&trend" = "order" %then n;
     %else &groupvar;;

/* Find number of strata */

proc sort data=xx; 
     by &stratvar;
data xx; 
     set xx;
     by &stratvar;
     retain stratn 0 ;
     if first.&stratvar then do;
          stratn+1;
          call symput('stratcnt', left(stratn));
     end;
run;

/*     Start loop on strata        */

%do ii = 1 %to &stratcnt;

/*    Form stratum subset, find number of 
      groups, number in each group, and      
      group weights in stratum               */

     data x; 
          set xx; 
          if stratn = &ii;
          call symput('stratval', &stratvar);
     run;
     proc freq; 
          table &groupvar/ noprint out= counts;
     proc sort data=x; 
          by &groupvar;
     data x; 
          set x;
          by &groupvar;
          retain grpn 0 ;
          if first.&groupvar then do;
                         grpn+1;
                call symput('grpcount', left(grpn));
                call symput('grpname'||left(grpn), &groupvar);
               end;
     run;
     data grpnames; 
          set x; 
          by &groupvar;
          keep &groupvar grpn;
          if first.&groupvar;
     data grpwts; 
          set grpnames;
          keep
          %if "&trend" = "order" %then grpn;
          %else &groupvar;;

/*    Create table         */

     proc sort data=x; 
          by descending &time;
     data y; 
          set x; 
          keep r1-r&grpcount rtot;
          array r{*} r1-r&grpcount;
          retain r1-r&grpcount rtot 0;
          %let countsq = %eval(&grpcount*&grpcount);
          r{grpn}+1;
          rtot+1;
     data x; 
          merge x y;
     proc sort; 
          by &time;
     data x; 
          set x; 
          by &time;
          array d{*} d1-d&grpcount;
          retain d1-d&grpcount dtot;
          if first.&time then do i=1 to &grpcount;
               d{i}=0;
               dtot=0;
               end;
          if &event not in (&censval) then do;
               d{grpn}+1;
               dtot+1;
               end;
          if last.&time then output;
     data x; 
          set x; 
          if dtot>0; 
          retain km km_  1; 
          all=1;
          array e{*} e1-e&grpcount;
          array diff{*} diff1-diff&grpcount;
          array r{*} r1-r&grpcount;
          array d{*} d1-d&grpcount;
          array wdiff{*} wdiff1-wdiff&grpcount;
          array s{*} sum1-sum&grpcount;
          array cov{&grpcount, &grpcount} cov1-cov&countsq;
          array sumcov{&grpcount,&grpcount}       sumcov1-sumcov&countsq;
          if _n_ = 1 then km_ = 1;
          else km_ = km;
          km=km*(rtot-dtot)/rtot;
          do j=1 to &grpcount;
               e{j} = dtot*r{j}/rtot;
               diff{j} = d{j} - e{j};
               if "&method"="logrank" then w=1;
               if "&method"="gehan" then w=rtot;
               if "&method"="tarone" then w=sqrt(rtot);
               if "&method"="harrington" then w=km_**&rho;
               wdiff{j} = w*diff{j};
               s{j} + wdiff{j};
          do l=1 to &grpcount;
                    if dtot=1 then c=1; else                     
                       c=(rtot-dtot)/(rtot-1);
                    if j=l then                                                
                       cov{j,l}=w**2*(dtot*(rtot*r{j}-r{j}**2)*c)/rtot**2;
                    else                     
                       cov{j,l}=-w**2*(r{j}*r{l}*dtot*c)/rtot**2;
                    sumcov{j,l}+cov{j,l};
                    end;
           end;

   /*       Sum over times and reformat for printout        */

     proc means sum noprint; 
          var d1-d&grpcount e1-e&grpcount diff1-diff&grpcount
          wdiff1-wdiff&grpcount;
          output out = out sum=;
     data out; 
     set out;
          array e{*} e1-e&grpcount;
          array d{*} d1-d&grpcount;
          array difff{*} diff1-diff&grpcount;
          array wdif{*} wdiff1-wdiff&grpcount;
          do j = 1 to &grpcount;
               group = j;
               events = d{j};
               expected = e{j};
               diff = difff{j};
               wdiff = wdif{j};
               output;
               end;
          label wdiff = 'Weighted Diff';
          label events = 'Events';
          label expected = 'Expected';
          label diff = 'Diff';
     data xxx; 
          merge out grpnames counts;
/* ad
     proc print l noobs; 
          var &groupvar count percent events expected diff wdiff;
          sum  count events;
          title1 'Summary of Events vs Expected';
          %if "&stratvar" ne "_none_" %then title2 "&stratvar = &stratval";;
          title3 "Method = &method, Rho=&rho";
*/
     run;

/*      Accumulate vectors and matrices for pooled stats  */

     %if "&ii" = "1" %then %do;
          data pooled; 
               set xxx; 
          %end;
     %else %do; 
          data pooled; 
          set pooled xxx; 
          %end;
     data x; 
          set x;
     proc sort; 
          by all;
     data s (keep = sum1-sum&grpcount) cov (keep =               
                 col1-col&grpcount);
          set x; 
          by all; 
          if last.all;
          array s{*} sum1-sum&grpcount;
          array sumcov{&grpcount, &grpcount}                   
                      sumcov1-sumcov&countsq;
          array col{*} col1-col&grpcount;
          output s;
          do j=1 to &grpcount;
               do l=1 to &grpcount;
                    col{l}=sumcov{j,l};
                    end;
          output cov;
          end;
     data yy; 
          merge grpnames cov;

                              /* Give columns of covariance matrix group names */
     %do j = 1 %to &grpcount;
          label col&j = "&&grpname&j";
          %end;
/* ad
     proc print l noobs; 
          var &groupvar col1-col&grpcount;
          title1 'Covariance Matrix';
          %if "&stratvar" ne "_none_" %then title2 "&stratvar=
               &stratval";;
          title3 "Method = &method, Rho=&rho";
*/
     %if "&ii" = "1" %then %do; 
          data poolcov; 
               set yy; 
          %end ;
     %else %do; 
          data poolcov; 
          set poolcov yy; 
          %end;

/*     Use proc iml to do matrix calculations 
       for test statistic.                       */  

     proc iml;
          reset noprint;
          use s;
          read all into x;
          use cov;
          read all into v;
          use grpwts;
          read all var _all_ into grpwts;

/*      Omit first row and column      */

          xx=x[1:1,2:&grpcount];
          vv=v[2:&grpcount,2:&grpcount];
          stat= xx*inv(vv)*xx`;
          df = &grpcount - 1;
          p_val = 1-probchi(stat,df);
          results = stat||df||p_val;
          cols={ChiSquare df p_value};
          title1 ' ';
          %if "&stratvar" ne "_none_" %then title1 "&stratvar=
               &stratval";;
          title2 "Method = &method, Rho=&rho";
          print results[colname=cols];

/*     Test for trend.        */

          if %eval(&grpcount) > 2 then do;
               wts=grpwts[2:&grpcount, 1:1];
               xxx=xx*wts;
               vvv=wts`*vv*wts;
               stat = xxx*xxx/vvv;
               df=1;
               p_val= 1-probchi(stat,df);
               trend = stat||df||p_val;
               print trend[colname=cols];
               end;
          quit;
     %end; 

 /* end of loop on strata */

/*     Pooled results if stratified analyis   */

%if "&stratvar" ne "_none_" %then %do;
     proc freq data=xx;
          table &groupvar / noprint out=counts;
     proc sort data=pooled; 
          by &groupvar;
proc print data=pooled;
     proc means noprint sum data=pooled; 
          var count events expected diff wdiff;
          by group; 
          output out=pooled1 sum=;
     data; 
          merge pooled1 grpnames counts;
/* ad
     proc print l noobs; 
            var &groupvar count percent events expected diff wdiff;
          sum count events;
          title1 'Summary of Events vs Expected';
          title2 "Pooled Over All Values of &stratvar";
*/
proc print data=poolcov;
     proc sort data=poolcov; 
          by &groupvar;
     proc means noprint sum data=poolcov;
          var col1-col&ngrps;
          by &groupvar; 
          output out=pooled2 sum=;
/* ad
     proc print l noobs; 
          var &groupvar col1-col&ngrps;
          title1 'Covariance Matrix';
          title2 "Pooled Over All Values of &stratvar";
*/
     data pooled2; 
          set pooled2; 
          keep col1-col&ngrps;
     run;
     proc iml;
          reset noprint;
          use pooled1;
          read all var {wdiff} into x;
          use pooled2;
          read all into v;
print x v;
          xx=x[2:&ngrps,1:1];
          vv=v[2:&ngrps,2:&ngrps];
print xx vv;
          stat = xx`*inv(vv)*xx;
          df = &ngrps - 1;
          p_val=1-probchi(stat,df);
          cols={ChiSquare df p_value};
          title1 'Pooled Results';
          title2 "Method = &method, Rho=&rho";
          results = stat||df|| p_val;
          print results[colname=cols];

/*     Test for trend.         */

          if %eval(&ngrps) > 2 then do;
               use groupwts;
               read all var _all_ into weights;
               wts = weights[2:&ngrps, 1:1];
               xtrend = xx`*wts;
               vtrend = wts`*vv*wts;
               stattrnd = xtrend**2/vtrend;
               p_valtrd = 1-probchi(stattrnd,1);
               df=1;
               trend=stattrnd||df||p_valtrd;
               print trend[colname=cols];
               run;
               end;
     %end;
%mend;
