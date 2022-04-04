
data{
 int<lower=1> nc; // number of control subjects
 int<lower=1> nt; // number of treated subjects
 int<lower=0> ncovU; // number of covariates
 int<lower=0> ncovY; // number of covariates  

 real YctlY[nc]; // control outcomes
 real YtrtY[nt]; // treatment outcomes
 
// real YctlU[nc]; // control outcomes <- this doesn't exist yet, right?
 real YtrtU[nt]; // treatment outcomes
 
 matrix[nc,ncovU] XctlU; // covariate matrix for controls
 matrix[nt,ncovU] XtrtU; //covariate matrix for treated

 matrix[nc,ncovY] XctlY; // covariate matrix for controls
 matrix[nt,ncovY] XtrtY; //covariate matrix for treated

 int<lower=0,upper=1> bottomOuter[nt]; // 1=firster 0=skipper <- WHat this the difference between this and YtrtU
}


parameters{
 real alphaTBO; // Y-intercept for treatment Bottom Outers
 real alphaTNBO; // Y-intercept for treatment Not  Bottom Outers
 real alphaCBO; // Y-intercept for control  Bottom Outers
 real alphaCNBO; // Y-intercept for control Not Bottom Outers

 real alphaU; // intercept for usage model
 vector[ncovU] betaU; // coefficients for covariates in usage model
 vector[ncovY] betaY; // coefficients for covariates in outcome model <- add dummy codes HERE, take out virtual (use model.matrix() function to get predictor matrix)
                                                                        // ^ I'm not sure when to use the model.matrix() function bc the mode is inside STAN

 // residual standard deviation:
 real<lower=0> sigTBO;
 real<lower=0> sigTNBO;
 real<lower=0> sigCBO;
 real<lower=0> sigCNBO;
}


transformed parameters{
 //declare new parameters:
 real bottomOuterATE; //Avg. Effect for skippers
 real notbottomOuterATE; //Avg. Effect for firsters
 real ATEdiff; //Difference btw Avg Effects
 vector[nc] piC; //Pr(notbottomOuterATE) for controls
 vector[nt] piT; //Pr(bottomOuter) for treateds

 //define new parameters:
 bottomOuterATE=alphaTBO-alphaCBO; // effect of treatement on bottom outers
 notbottomOuterATE=alphaTNBO-alphaCNBO; // effect of treatment on not bottom outers
 ATEdiff=bottomOuterATE-notbottomOuterATE; // differences in effects
 piC=inv_logit(alphaU+XctlU*betaU);
 piT=inv_logit(alphaU+XtrtU*betaU);// prop of bing a bottom outer
}


model{
 // vectors of intercepts and residual SDs for treated
 // students. useful for vectorizing:
 vector[nt] alphaT;
 vector[nt] sigT;

 // The ? functions as if-else A?B:C returns B if A=1 and C otherwise https://mc-stan.org/docs/2_29/reference-manual/conditional-statements.html
 for(i in 1:nt){
   alphaT[i] = bottomOuter[i]?alphaTBO:alphaTNBO;
   sigT[i] = bottomOuter[i]?sigTBO:sigTNBO;
 }


// producing error
 //bottomOuterATE~normal(0,.5);
 //notbottomOuterATE~normal(0,.5);


 bottomOuter~bernoulli(piT); // model for who is a bottom outer
 YtrtY~normal(alphaT+XtrtY*betaY,sigT); // model for outcomes in treatment group


 for(i in 1:nc)
  target += log_sum_exp(
   log(piC[i]) + normal_lpdf(YctlY[i] // I'm not sure what this should
          | alphaCBO+XctlY[i,]*betaY,sigCBO),
   log((1-piC[i])) + normal_lpdf(YctlY[i] |alphaCNBO+XctlY[i,]*betaY,sigCNBO));
} 

generated quantities{
 int<lower=0,upper=1> first_repC[nc]; //M_T for controls
 real Ytrt_rep[nt]; // outcomes for treateds
 real Yctl_rep[nc]; // outcomes for controls

 for(i in 1:nt){
  if(bottomOuter[i]==1)
   Ytrt_rep[i]=normal_rng(alphaTBO+XtrtY[i,]*betaY,sigTBO);
  else
   Ytrt_rep[i]=normal_rng(alphaTNBO+XtrtY[i,]*betaY,sigTNBO);
 }
 for(i in 1:nc){
  first_repC[i]=bernoulli_rng(piC[i]);
  if(first_repC[i]==1)
   Yctl_rep[i]=normal_rng(alphaCBO+XctlY[i,]*betaY,sigCBO);
  else
   Yctl_rep[i]=normal_rng(alphaCNBO+XctlY[i,]*betaY,sigCNBO);
 }
}

