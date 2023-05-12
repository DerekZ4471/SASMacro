%macro BayesFutilityBin(data,par,delta,eta,alpha,prob);
    /*********************************************************
    Inputs:
    
    DATA    = Data set to be analyzed (includes number of patients and
              observed event counts in two treatment groups at each 
              interim look)
    
    PAR     = Name of an input data set with projected sample size in 
              each treatment group and parameters of prior distributions
    
    DELTA   = Clinically significant difference (required by Bayesian 
              predictive probability method and ignored by predictive 
              power method)
    
    ETA     = Confidence level of Bayesian predictive probability method
              (required by Bayesian predictive probability method and 
              ignored by predictive power method)
    
    ALPHA   = One-sided Type I error probability of the significance 
              test carried out at the end of the trial (required by 
              predictive power method and ignored by Bayesian predictive
              probability method
    
    PROB    = Name of an output data set containing predictive power and 
              predictive probability at each interim look
    
    *********************************************************/
    proc iml;
        start integral(p) global(ast,bst);
            i=p**(ast[1]-1)*(1-p)**(bst[1]-1)*
                probbeta(p-&delta,ast[2],bst[2]);
            return(i);
        finish;
        start beta(a,b);
            beta=exp(lgamma(a)+lgamma(b)-lgamma(a+b));
            return(beta);
        finish;
        use &data;
        read all var {n1 count1 n2 count2} into data;
        n=data[,1]||data[,3];    
        s=data[,2]||data[,4];
        m=nrow(n);
        use &par;
        read all var {nn1 nn2 alpha1 alpha2 beta1 beta2} into par;
        nn=t(par[1:2]);
        a=t(par[3:4]);
        b=t(par[5:6]);    
        t=j(1,2,0);
        output=j(m,4,0);
        range=j(1,2,0);
        range[1]=&delta; range[2]=1;
        do i=1 to m;
        output[i,1]=i;
        output[i,2]=(n[i,1]+n[i,2])/(nn[1]+nn[2]);
        do t1=0 to nn[1]-n[i,1];
        do t2=0 to nn[2]-n[i,2];
            t[1]=t1; t[2]=t2;
            ast=s[i,]+t+a; 
            bst=nn-s[i,]-t+b;                
            b1=beta(ast,bst);                    
            b2=exp(lgamma(nn-n[i,]+1)-lgamma(nn-n[i,]-t+1)-lgamma(t+1));
            b3=beta(s[i,]+a,n[i,]-s[i,]+b);
            pr=(b1#b2)/b3;
            mult=pr[1]*pr[2];
            p=(s[i,]+t)/nn;
            ave=(p[1]+p[2])/2;
            aven=(nn[1]+nn[2])/2;
            teststat=(p[1]-p[2])/sqrt(2*ave*(1-ave)/aven);        
            output[i,3]=output[i,3]+(teststat>probit(1-&alpha))*mult;
            call quad(value,"integral",range) eps=1e-8;
            output[i,4]=output[i,4]+(value/beta(ast[1],bst[1])>&eta)*mult;
        end;
        end;
        end;    
        varnames={"Analysis", "Fraction", "PredPower", "PredProb"};
        create &prob from output[colname=varnames];
        append from output;
        quit;
    data &prob;
        set &prob;
        format Fraction 4.2 PredPower PredProb 6.4;
        label Fraction='Fraction of total sample size'          
              PredPower='Predictive power'
              PredProb='Predictive probability';
    %mend BayesFutilityBin;
    