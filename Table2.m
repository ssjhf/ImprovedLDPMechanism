clear
clc

epsilon=[0.01 0.05 0.25 0.5];
pi_A=0.1;
var=0.1;
x=exp(epsilon);


%modified Warner
N1=((x./(x-1)./(x-1)))/var;

%modified Simmons
N2=((x./(x-1)./(x-1)))/var;

%modified Christofides
p2=0.01;
N3=(((x+1).*(x+1)./(x-1)./(x-1)/(1-p2)-1)/4)/var;

%%Theory improved Christofides
p2=0.01;
N4=pi_A*(1-pi_A)*((x+1).*(x+1)./(x-1)./(x-1)/(1-p2)-1)/var+1;

N_result=ceil([N1;N2;N3;N4]);