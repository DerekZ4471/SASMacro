/***************************************************

MinRisk macro

Dmitrienko, A., Koch, G. G., and SAS Institute (eds.) (2017), Analysis of clinical trials using SAS: a practical guide, Cary, NC: SAS Institute.

Macro written by Alan Cantor

***************************************************/

%macro MinRisk(dataset);
/*
Inputs:

DATASET  = Data set with event rates observed in each stratum.
*/
proc iml;
    use &dataset;
    read all var {event1 noevent1 event2 noevent2} into data;    
    m=nrow(data);
    p=j(m,2,0); n=j(m,2,0);
    n[,1]=data[,1]+data[,2];
    n[,2]=data[,3]+data[,4];
    total=sum(n);
    p[,1]=data[,1]/n[,1];
    p[,2]=data[,3]/n[,2];
    delta=p[,1]-p[,2];
    v=p[,1]#(1-p[,1])/n[,1]+p[,2]#(1-p[,2])/n[,2];
    pave=(p[,1]#n[,1]+p[,2]#n[,2])/(n[,1]+n[,2]);
    v0=pave#(1-pave)#(1/n[,1]+1/n[,2]);
    alpha=delta*sum(1/v)-sum(delta/v);    
    c=1+alpha*sum((n[,1]+n[,2])#delta/total);
    h=diag(v*sum(1/v))+alpha*delta`;
    wmr=inv(h)*c;
    dmr=wmr`*delta;  
    zmr1=abs(dmr)-3/(16*sum(n[,1]#n[,2]/(n[,1]+n[,2])));
    zmr2=sqrt(sum(wmr#wmr#v0));
    zmr=zmr1/zmr2;
    pmr=2*(1-probnorm(zmr));
    title={"Estimate", "Statistic", "P-value"};
    minrisk=dmr||zmr||pmr;
    print minrisk [colname=title format=best6.];
    win=(1/v)/sum(1/v);
    din=win`*delta;    
    zin=abs(din)/sqrt(sum(win#win#v0));
    pin=2*(1-probnorm(zin));
    invar=din||zin||pin;
    print invar [colname=title format=best6.]; 
    wss=(n[,1]#n[,2]/(n[,1]+n[,2]))/sum(n[,1]#n[,2]/(n[,1]+n[,2]));
    dss=wss`*delta;    
    zss=abs(dss)/sqrt(sum(wss#wss#v0));
    pss=2*(1-probnorm(zss));
    ssize=dss||zss||pss;
    print ssize [colname=title format=best6.]; 
    quit;
%mend;
