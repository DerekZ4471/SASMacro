
%macro SeqPower(data, analysis, fraction, out);
    /*********************************************************
    Inputs:
    
    DATA     = Data set to be analyzed (includes interim analysis number, 
               treatment indicator, number of patients and number of patients
               with the event of interest)
    
    ANALYSIS = Number of the current interim analysis
    
    FRACTION = Information fraction at the current interim analysis
    
    OUT      = Data set for storing conditional and predictive power
    
    *********************************************************/
    
        proc seqdesign altref=0.09;
            design method=errfuncobf
            nstages=2
            info=cum(&fraction 1)
            alt=upper
            stop=reject
            alpha=0.1
            beta=0.2;
        samplesize model(ceiladjdesign=include)
                         =twosamplefreq(nullprop=0.70 test=prop ref=avgprop);
        ods output adjustedboundary=obf_boundary;
        run;
    
        ods select none;
    
          proc genmod data=&data( where=( _stage_ = &analysis));
              model nsurv/nobs= treatment / dist=bin link=identity;
              ods output nobs=ntotal parameterestimates=parmest;
          run;
    
          data ntotall;
             set ntotal;
             nobs= n;
             if (label='Number of Trials');
             keep nobs;
          run;
    
          data parmsl;
             set parmest;
             if parameter='treatment';
             _scale_ ='MLE';
             keep _scale_ parameter estimate;
          run;
    
          data parms;
              merge ntotall parmsl;
              _stage_ = 1;
          run;
    
          proc seqtest boundary=obf_boundary
              parms(testvar=treatment infovar=nobs)=parms
              boundaryadj=errfuncobf
              boundarykey=alpha
              nstages=2 order=stagewise
              condpower(type=finalstage) 
              predpower;
          ods output condpower=cp predpower=pp;
          run;
    
        data cnull;
            set cp;
            if ref="Null";
            cnull=cpower;
            keep cnull;
        run;    
    
        data cmle;
            set cp;
            if ref="MLE";
            cmle=cpower;
            keep cmle;
        run;    
            
        data calt;
            set cp;
            if ref="Alternative";
            calt=cpower;
            keep calt;
        run;    
            
        data &out;
            merge cnull cmle calt pp;   
            analysis=&analysis;
        run;    
        
        ods select all;
    
    %mend SeqPower; 
    