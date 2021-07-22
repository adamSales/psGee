data {
     int<lower=1> nctl;
     int<lower=1> ntrt;
     int<lower=1> p;
     int<lower=1> nscl;
     matrix[ntrt,p] xt;
     matrix[nctl,p] xc;	

     int<lower=1,upper=nscl> sclt[ntrt];
     int<lower=1,upper=nscl> sclc[nctl];     

     vector[ntrt] Ytrt;
     vector[nctl] Yctl;

     int<lower=0,upper=1> St[ntrt];

}
parameters {
 
 vector[p] a2;

 real a00;
 real a01;
 real a10;
 real a11;

 real b0;
 vector[p] b1;

 real sclEffY[nscl];
 real sclEffS[nscl];
 
 real<lower=0> sigSclY;
 real<lower=0> sigSclS;
 real<lower=0> sigt;
 real<lower=0> sigc;
}

transformed parameters{
 real eff0=a10-a00;
 real eff1=a10+a11-a01;
 real effDiff=eff1-eff0;
}

model{

 vector[ntrt] xb1;
 vector[nctl] xb2;
 vector[ntrt] xb3;
 vector[nctl] prob;

 a00~std_normal();
 a10~std_normal();
 a01~std_normal();
 a11~std_normal();

 sclEffY~normal(0,sigSclY);
 sclEffS~normal(0,sigSclS);

 xb1=a10+a11*to_vector(St)+xt*a2;
 xb2=xc*a2;
 xb3=b0+xt*b1;

 for(i in 1:ntrt){
  xb1[i]+=sclEffY[sclt[i]];
  xb3[i]+=sclEffS[sclt[i]];
 }
 for(i in 1:nctl){
  xb2[i]+=sclEffY[sclc[i]];
 }
  
 St~bernoulli_logit(xb3);
 Ytrt~normal(xb1,sigt);

 prob=b0+xc*b1;

 for(i in 1:nctl){
       target+=log_sum_exp(
        log(prob[i])+normal_lpdf(Yctl[i]|a01+xb2[i],sigc),
        log(1-prob[i])+normal_lpdf(Yctl[i]|a00+xb2[i],sigc)
	);
 }
}

