clear
clc

%Warner by numerical simulation
n=1;
N=9;
data=zeros(1,N);
data(1:n)=1;
data=data(randperm(N));
epsilon=0.05:0.01:0.5;
REPEAT_TIMES=10000;

p_Warner=exp(epsilon)./(exp(epsilon)+1);
var_Warner_experiment=zeros(length(epsilon),1);
for j=1:length(epsilon)
    j
    E_pi=zeros(REPEAT_TIMES,1);
    for repeat=1:REPEAT_TIMES
        rand_num=rand(1,N);
        dataWarner=(rand_num<p_Warner(j)).*data+(rand_num>p_Warner(j)).*(1-data);
        lambda=mean(dataWarner);
        E_pi(repeat)=(lambda-(1-p_Warner(j)))/(2*p_Warner(j)-1);
    end
    var_Warner_experiment(j)=var(E_pi);
end
plot(epsilon,var_Warner_experiment,'-db','LineWidth',1)
hold on

%modified Christofides by numerical simulation
clear
clc

n=1;
N=9;
data=zeros(1,N);
data(1:n)=1;
data=data(randperm(N));
P2=0.36;
epsilon=0.05:0.01:0.5;
L=3;
REPEAT_TIMES=10000;
x=exp(epsilon);

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
            elseif (rand1<round(p1(j)*N)/N+P2)&&(rand1>round(p1(j)*N)/N)
                card_draw(i)=2;
            else
                card_draw(i)=3;
            end
        end
        X=(L+1-card_draw).*data+card_draw.*(1-data);
        lambda(r)=mean(X);
        piA_estimation(r)=(lambda(r)-EY(j))/(L+1-2*EY(j));
    end
    var_Christofides_experiment(j)=var(piA_estimation);
end
plot(epsilon,var_Christofides_experiment,'-*r','lineWidth',1)

%the improved Christofides by numerical simulation
clear
clc
n=1;
pi_A=n/9;
N=9;
data=zeros(1,N);
data(1:n)=1;
data=data(randperm(N));
P2=0.36;
epsilon=0.05:0.01:0.5;
L=3;
REPEAT_TIMES=10000;
x=exp(epsilon);

p1=(1-P2)./(x+1);
p3=x.*p1;
p=[p1;1-p1-p3;p3]';
lambda=zeros(REPEAT_TIMES,1);
pi_A_estimate=zeros(REPEAT_TIMES,1);
for j=1:length(epsilon)
    for l=1:L
        q(j,l+1)=sum(p(j,1:l));
    end
    EY(j)=sum([1:L].*p(j,:));
end
card_draw=zeros(N,1);


var_the_improved_Christofides_experiment=zeros(46,1);
for j=1:46
    j
    lambda=zeros(REPEAT_TIMES,1);
    pi_A_estimate=zeros(REPEAT_TIMES,1);
    for r=1:REPEAT_TIMES
        card=zeros(1,N);
        for l=1:L
            card(1+round(q(j,l)*N):round(q(j,l+1)*N))=[ones(round(q(j,l+1)*N-round(q(j,l)*N)),1)]*l;
        end
        card_draw=card(randperm(N));
        X=(L+1-card_draw).*data+card_draw.*(1-data);
        lambda(r)=mean(X);
        pi_A_estimate(r)=(lambda(r)-EY(j))/(L+1-2*EY(j));
    end
    var_the_improved_Christofides_experiment(j)=var(pi_A_estimate);
end
plot(epsilon,var_the_improved_Christofides_experiment,'-sk','lineWidth',1)

xlabel('privacy budget (\epsilon)','FontSize',14)
ylabel('variance','FontSize',14)
title('comparison of four mechanisms')
set(gca,'FontSize',14);
xlim([0.05 0.5])
set(gca,'XTick',[0.05:0.05:0.5])
legend({'modified Warner/Simmons','modified Christofides(L=3)','the improved Christofides(L=3)'},'FontSize',14)