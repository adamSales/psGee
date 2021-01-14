data {
     int<lower=1> nctl;
     int<lower=1> ntrt;
     vector[ntrt] x1t;
     vector[ntrt] x2t;
     vector[nctl] x1c;
     vector[nctl] x2c;
     vector[ntrt] Ytrt;
     vector[nctl] Yctl;

     int<lower=0,upper=1> St[ntrt];

}
parameters {
 
 real b1y;
 real b2y;

 real mu00;
 real mu01;
 real mu10;
 real mu11;

 real b0s;
 real b1s;
 real b2s;

 real<lower=0> sigt;
 real<lower=0> sigc;
}

model{

 mu00~std_normal();
 mu10~std_normal();
 mu01~std_normal();
 mu11~std_normal();

 St~bernoulli_logit(b0s+b1s*x1t+b2s*x2t);
 Ytrt~normal(mu10+b1y*x1t+b2y*x2t+(mu11-mu10)*to_vector(St),sigt);

 for(i in 1:nctl){
       real prob=inv_logit(b0s+b1s*x1c[i]+b2s*x2c[i]);
       target+=log_sum_exp(
        log(prob)+normal_lpdf(Yctl[i]|mu01+b1y*x1c[i]+b2y*x2c[i],sigc),
        log(1-prob)+normal_lpdf(Yctl[i]|mu00+b1y*x1c[i]+b2y*x2c[i],sigc)
	);
 }
}

