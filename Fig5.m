clear
clc

load health_insurance.dat
health_insurance=health_insurance';
N = length(health_insurance);

%modified Warner
epsilon=0.05:0.01:0.5;
REPEAT_TIMES=2;
p_Warner=exp(epsilon)./(exp(epsilon)+1);
var_Warner_experiment=zeros(length(epsilon),1);
for j=1:length(epsilon)
    j
    E_pi=zeros(REPEAT_TIMES,1);
    for repeat=1:REPEAT_TIMES
        rand_num=rand(1,N);
        health_insurance_Warner=(rand_num<p_Warner(j)).*health_insurance+(rand_num>p_Warner(j)).*(1-health_insurance);
        lambda=mean(health_insurance_Warner);
        E_pi(repeat)=(lambda-(1-p_Warner(j)))/(2*p_Warner(j)-1);
    end
    var_Warner_experiment(j)=var(E_pi);
end
plot(epsilon,var_Warner_experiment,'-db','LineWidth',1,'MarkerSize',8)
hold on

%modified Simmons
% clear
% clc
% 
% load health_insurance.dat
% health_insurance=health_insurance';
% N = length(health_insurance);
% epsilon=0.05:0.01:0.5;
% pi_B = 0.5;
% REPEAT_TIMES=2;
% var_Simmons_experiment=zeros(length(epsilon),1);
% p_Simmons=(1-pi_B)*(exp(epsilon)-1)./((1-pi_B)*(exp(epsilon)-1)+1);
% 
% for j=1:length(epsilon)
%     j
%     E_pi=zeros(REPEAT_TIMES,1);
% 
%     for repeat=1:REPEAT_TIMES
%         rand_num1=rand(1,N);
%         rand_num2=rand(1,N);
%         health_insurance_Simmons=(rand_num1<p_Simmons(j)).*((rand_num1<p_Simmons(j))+(1-(rand_num1<p_Simmons(j))).*(rand_num2<pi_B))+(rand_num1>p_Simmons(j)).*(rand_num1>p_Simmons(j)).*(rand_num2<pi_B);
%         lambda=mean(health_insurance_Simmons);
%         E_pi(repeat)=(lambda-(1-p_Simmons(j)*pi_B))/p_Simmons(j);
%     end
%     var_Simmons_experiment(j)=var(E_pi);
% end
% plot(epsilon,var_Simmons_experiment,'-','LineWidth',1,'Color',[0 0.5 0],'MarkerSize',8)

%modified Christofides
clear
clc

load health_insurance.dat
health_insurance=health_insurance';
N = length(health_insurance);
L=3;
REPEAT_TIMES=2;
epsilon=0.05:0.01:0.5;
x=exp(epsilon);
P2=0.01;
p1=(1-P2)./(x+1);
p3=x.*p1;
p=[p1;1-p1-p3;p3]';
for j=1:length(epsilon)
    for l=1:L
        q(j,l+1)=sum(p(j,1:l));
    end
    EY(j)=sum([1:L].*p(j,:));
end

var_Christofides_experiment=zeros(46,1);

for j=1:46
    j
    lambda=zeros(REPEAT_TIMES,1);
    piA_estimation=zeros(REPEAT_TIMES,1);

    for r=1:REPEAT_TIMES
        card_draw=zeros(1,N);
        for i=1:N
            rand1=rand;
            if rand1<round(p1(j)*N)/N
                card_draw(i)=1;
            elseif (rand1>(round(p1(j)*N)/N))&&(rand1<(round(p1(j)*N)/N+P2))
                card_draw(i)=2;
            else
                card_draw(i)=3;
            end
        end
        X=(L+1-card_draw).*health_insurance+card_draw.*(1-health_insurance);
        lambda(r)=mean(X);
        piA_estimation(r)=(lambda(r)-EY(j))/(L+1-2*EY(j));
    end
    var_Christofides_experiment(j)=var(piA_estimation);
end
plot(epsilon,var_Christofides_experiment,'-or','lineWidth',1,'MarkerSize',8)

%the improved Christofides
clear
clc

load health_insurance.dat
health_insurance=health_insurance';
N = length(health_insurance);
epsilon=0.05:0.01:0.5;
L=3;
REPEAT_TIMES=2;
x=exp(epsilon);
P2=0.01;
p1=(1-P2)./(x+1);
p3=x.*p1;
p=[p1;1-p1-p3;p3]';
lambda=zeros(REPEAT_TIMES,1);
pi_A_estimation=zeros(REPEAT_TIMES,1);
for j=1:length(epsilon)
    for l=1:L
        q(j,l+1)=sum(p(j,1:l));
    end
    EY(j)=sum([1:L].*p(j,:));
end
card_draw=zeros(1,N);


var_the_improved_Christofides_experiment=zeros(46,1);
for j=1:46
    j
    lambda=zeros(REPEAT_TIMES,1);
    pi_A_estimation=zeros(REPEAT_TIMES,1);
    for r=1:REPEAT_TIMES
        card=zeros(1,N);
        for l=1:L
            card(1+round(q(j,l)*N):round(q(j,l+1)*N))=[ones(round(q(j,l+1)*N-round(q(j,l)*N)),1)]*l;
        end
        card_draw=card(randperm(N));
        X=(L+1-card_draw).*health_insurance+card_draw.*(1-health_insurance);
        lambda(r)=mean(X);
        pi_A_estimation(r)=(lambda(r)-EY(j))/(L+1-2*EY(j));
    end
    var_the_improved_Christofides_experiment(j)=var(pi_A_estimation);
end
plot(epsilon,var_the_improved_Christofides_experiment,'-sk','lineWidth',1,'MarkerSize',8)

xlabel('privacy budget (\epsilon)','FontSize',14)
ylabel('variance','FontSize',14)
title('comparison of four mechanisms')
set(gca,'FontSize',14);
xlim([0.05 0.5])
set(gca,'XTick',[0.05:0.05:0.5])
legend({'modified Warner/Simmons','modified Christofides(L=3)','the improved Christofides(L=3)'},'FontSize',14)