clc;
clear variables;
beep off;

t0=clock; 

strike=100; %Strike Price
n=5; %Time steps per year    
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

paths=1000000; %Number of paths
ES=zeros(paths,1); %S(T) stored for each path

%==================================
%Euler scheme
%==================================
for i=1:paths
    
    k=2;
    z1=normrnd(0,1,1,(T/int+1));
    z2=normrnd(0,1,1,(T/int+1));
    z2=(rho*z1+sqrt(1-rho^2)*z2);
    
    for j=int:int:T
        
        V(k)=V(k-1)+(kappa/n)*( theta - V(k-1)*(V(k-1)>0) )+sigma*(sqrt( V(k-1)*(V(k-1)>0)/n ) )*z1((k-1));
        S(k)=exp(log(S(k-1)) + (1/n)*(r- ( 1/2 )*V(k-1)*(V(k-1)>0) ) +sqrt( V(k-1)*(V(k-1)>0)/n )*z2((k-1)) );
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