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
 vector[nctl] probY00;
 vector[nctl] probY01;	

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
 Ytrt~bernoulli(xb1<0?0:xb1>1?1:xb1);

 probS0=inv_logit(b0+xc*b1);
 probY00=(a00+xb2)<0?0:(a00+xb2)>1?1:(a00+xb2);
 probY01=(a01+xb2)<0?0:(a01+xb2)>1?1:(a01+xb2);	


 for(i in 1:nctl){
       target+=log_sum_exp(
        log(prob[i])+bernoulli_lpdf(Yctl[i]|probY01),
        log(1-prob[i])+normal_lpdf(Yctl[i]|probY00)
	);
 }
}

