pi = 202; %recruitment rate
mu=0.0089; %natural death rate
kappa=0.156986; %rate of development of clinical symptoms
alpha=0.156986; %hospitalization rate of quarantine
phi=0.20619; %hospitalization rate of infectious
psi=0.001; %loss of infection acquired immunity
d1=0.04227; %death rate infectious
d2=0.027855; %death rate hospitalized
hnu=.6; %relative infectiousness e.g. 0=perfect intervention
qnu=.0;
gamma1=0.03521; %recovery rate of infectious
gamma2=0.027855; %recovery rate of hospitalized
mean_cont=0.9575;
transm_risk=0.20;
beta=mean_cont*transm_risk;
st=0.70; %succes rate of test
ss=0.782; %succes rate of sampling
rt=1; %result of test e.g. 1=Same day
tr=0.14552*rt; %testing rate e.g. T/N, 1=Everybody is tested
sigma=st*ss*tr; %quarantine rate

N=20000; %population


tspan = 0:1:300; % Time Span of the Simulation
xinit = [19992,0,8,0,0,0,N]; % Initial Conditons in order; 
                        % ('S(t)', 'E(t)', 'I(t)', 'Q(t)', 'H(t)', 'R(t)', 'N(t)')

f = @(t1,y1) [pi+psi*y1(6)-(beta*(y1(3)+hnu*y1(5))/y1(7))*y1(1)-mu*y1(1);
    (beta*(y1(3)+hnu*y1(5))/y1(7))*y1(1)-(kappa+sigma+mu)*y1(2);
    kappa*y1(2)-(gamma1+phi+mu+d1)*y1(3);
    sigma*y1(2)-(alpha+mu)*y1(4);
    alpha*y1(4)+phi*y1(3)-(gamma2+mu+d2)*y1(5);
    gamma1*y1(3)+gamma2*y1(5)-(psi+mu)*y1(6);
    pi-mu*y1(7)-(d1*y1(3)+d2*y1(5))];

beds=ones(300,1)* 66.4;

[t1,y1]=ode45(f, tspan, xinit);

plot(y1(:,:))
xlabel('Days')
ylabel('Number of Individuals')
legend('S(t)', 'E(t)', 'I(t)', 'Q(t)', 'H(t)', 'R(t)', 'N(t)')
%plot(y1(:,5))
hold on
%ylim([0 800])
%plot(beds)
        
