clc;
clear variables;
beep off;

t0=clock; 

strike=100; %Strike Price
n=6; %Time steps per year    
T=5; %Years to maturity

%=======================
%Heston model parameters
rho=-0.3;
sigma=1; 
kappa=2; 
theta=0.09;
r=0.05;
%=======================

int=1/n; %Time step size
S=zeros(T/int+1,1); %Price
V=zeros(T/int+1,1); %Volatility
S(1)=100; %Initial price
V(1)=0.09; %Initial volatility

paths=100000; %Number of paths
ES=zeros(paths,1); %S(T) stored for each path

%==================================
%Bernoulli scheme
%==================================
for i=1:paths
     k=2;

    mu1=0.675;
    u1=mu1+1/mu1;
    g1=rand(1,T/int+1);
    g1=(mu1+1/mu1)*(g1<(mu1^2/(1+mu1^2)));
    mean1=mu1;
    
    mu2=1;
    u2=mu2+1/mu2;
    g2=rand(1,T/int+1);
    g2=(mu2+1/mu2)*(g2<(mu2^2/(1+mu2^2)));
    mean2=mu2;

    g2=(rho*g1+sqrt(1-rho^2)*g2)-(mean1*rho+sqrt(1-rho^2)*mean2);

    for j=int:int:T
        
        V(k)= V(k-1)+(kappa/n)*(theta-V(k-1))+(1/sqrt(n))*sigma*sqrt(V(k-1))*(g1(k-1)-mean1);
        S(k)=S(k-1)*exp( (1/n)*(r-( 1/2 )*V(k-1)  ) +sqrt( V(k-1)/n )*g2(k-1) );
        k=k+1;
       
    end

  ES(i)=exp(-r*T)*max(S(k-1)-strike,0);
  
end

ES_mean=mean(ES) %Expected Discount Call Option Price
STD=std(ES)/sqrt(paths) %Standard deviation of ES_mean
Bias=ES_mean-34.9998 %Estimated Bias of ES_mean
Lower=ES_mean-1.96*STD % 95% lower CI on ES_mean
Upper=ES_mean+1.96*STD % 95% upper CI on ES_mean
ms = round(etime(clock,t0) * 1000) %Computation Time


