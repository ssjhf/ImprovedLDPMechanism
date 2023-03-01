clear
clc

N=10000;
epsilon=[0.01 0.05 0.25 0.5];
p2=[0.01 0.05 0.1];
for i=1:4
    for j=1:3
        left(i,j)=0.5-0.5*sqrt(1-((1-p2(j))./(1+p2(j)/4*(exp(epsilon(i))+1/exp(epsilon(i))-2))));
        right(i,j)=0.5+0.5*sqrt(1-((1-p2(j))./(1+p2(j)/4*(exp(epsilon(i))+1/exp(epsilon(i))-2))));
        length(i,j)=sqrt(1-(N-1)/N*((1-p2(j))./(1+p2(j)/4*(exp(epsilon(i))+1/exp(epsilon(i))-2))));
    end
end
length=roundn(length,-3)';