
%macro BayesFutilityCont(data,par,delta,eta,alpha,prob);
    /*********************************************************
    Inputs:
    
    DATA    = Data set to be analyzed (includes number of patients,
              estimated mean treatment effects and sample standard 
              deviations in two treatment groups at each interim look)
    
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
        use &data;
        read all var {n1 mean1 sd1 n2 mean2 sd2} into data;
        n=(data[,1]+data[,4])/2;
        s=sqrt(((data[,1]-1)#data[,3]#data[,3]+(data[,4]-1)#data[,6]#data[,6])/(data[,1]+data[,4]-2));	
        z=(data[,2]-data[,5])/(s#sqrt(1/data[,1]+1/data[,4]));    
        m=nrow(data);
        use &par;
        read all var {nn1 nn2 mu1 mu2 sigma} into par;
        nn=(par[,1]+par[,2])/2;
        sigma=par[,5];
        a1=(n-nn)/nn+(nn-n)/(nn#(1+(s/sigma)#(s/sigma)/n));
        b1=sqrt(n/(2*nn))#(nn-n)#(par[,3]-par[,4])/(1+n#(sigma/s)#(sigma/s));
        c1=1/(1+n#(sigma/s)#(sigma/s));
        output=j(m,4,0);
        output[,1]=t(1:m);
        output[,2]=n/nn;
        num1=sqrt(nn)#z#(1+a1)+b1/s-sqrt(n)*probit(1-&alpha);
        den1=sqrt((nn-n)#((nn-n)#(1-c1)+n)/nn);
        output[,3]=probnorm(num1/den1);      
        an=1/(1+n#(sigma/s)#(sigma/s));
        ann=1/(1+nn#(sigma/s)#(sigma/s));	
        b2=sqrt(n#nn/2)#(par[,3]-par[,4])/(1+n#(sigma/s)#(sigma/s));
        c2=1-1/sqrt(1+(s/sigma)#(s/sigma)/nn);
        num2=sqrt(nn)#z#(1-an)+b2/s-(&delta/s)#sqrt(n#nn/2)-sqrt(n)#(1-c2)*probit(&eta);
        den2=sqrt((1-ann)#(1-ann)#(nn-n)#((nn-n)#(1-an)+n)/nn);
        output[,4]=probnorm(num2/den2);      
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
    %mend BayesFutilityCont;
    