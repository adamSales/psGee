data {
     int<lower=1> nctl;
     int<lower=1> ntrt;
     int<lower=1> p;
     matrix[ntrt,p] xt;
     matrix[nctl,p] xc;	

     int<lower=0,upper=1> Ytrt[ntrt];
     int<lower=0,upper=1> Yctl[nctl];	

     int<lower=0,upper=1> St[ntrt];

     real<lower=0> lowBound;
     real<upper=1> upBound;

}
parameters {
 
 vector[p] a2;

 real a00;
 real a01;
 real a10;
 real a11;

 real b0;
 vector[p] b1;

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
 vector[nctl] probS0;
 //vector[nctl] probY00;
 //vector[nctl] probY01;	
 //real probS0;
 real probY00;
 real probY01;

 a00~std_normal();
 a10~std_normal();
 a01~std_normal();
 a11~std_normal();
 a2~std_normal();
 b0~std_normal();
 b1~std_normal();

 xb1=a10+a11*to_vector(St)+xt*a2;
 xb2=xc*a2;
 xb3=b0+xt*b1;

  
 St~bernoulli_logit(xb3);

 for(i in 1:ntrt)
  Ytrt[i]~bernoulli(xb1[i]<lowBound?lowBound:xb1[i]>upBound?upBound:xb1[i]);

 probS0=inv_logit(b0+xc*b1);	

 for(i in 1:nctl){
  probY00=(a00+xb2[i])<lowBound?lowBound:(a00+xb2[i])>upBound?upBound:(a00+xb2[i]);
  probY01=(a01+xb2[i])<lowBound?lowBound:(a01+xb2[i])>upBound?upBound:(a01+xb2[i]);	

  target+=log_sum_exp(
   log(probS0[i])+bernoulli_lpmf(Yctl[i]|probY01),
   log(1-probS0[i])+bernoulli_lpmf(Yctl[i]|probY00)
  );
 }
}

