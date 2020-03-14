
N=1000; % Total Initial Population
r=0.3; % Effectiveness of Isolation
l=0.5; % Relative Transmissibility of Isolated Individuals
mu=4.98*10^-5; % Natural Death Rate
alpha=1/10; % 
gamma=1/14; % Removal Rate of Infectious Individuals
gammar=1/10; % Removal Rate of Isolated Individuals
beta=0.3335; % Mean Transmission Rate
delta=0.15; % Probability of a Recovered Individuals becomes Susceptible
d=1/14; % Removal Rate of Recovered Individuals
T=0; % Number of Tests
T1=80;
ss=0.95; % Succes Rate of Sampling
st=0.85; % Succes Rate of Test
c=ss*st; % Succes Rate
r1=(1-r)*l; 
rate=T*c; 
rate1=T1*c;
tD=25; % Start of Drone Operation
te=50; % End of Drone Operation

f = @(t1,y1) [mu*N-mu*y1(1)-beta*(y1(4)+r1*y1(3))/(N-r*y1(3))*y1(1)+delta*d*y1(5);...
    (beta*(y1(4)+r1*y1(3))/(N-r*y1(3))*y1(1))-(mu+alpha*(1+rate/N))*y1(2);...
    alpha*rate/N*(y1(2)+y1(4))-(gammar+mu)*y1(3);...
    alpha*y1(2)-(gamma+mu+alpha*rate/N)*y1(4);...
    gammar*y1(3)+gamma*y1(4)-mu*y1(5)-delta*d*y1(5)];

g = @(t2,y2) [mu*N-mu*y2(1)-beta*(y2(4)+r1*y2(3))/(N-r*y2(3))*y2(1)+delta*d*y2(5);...
    (beta*(y2(4)+r1*y2(3))/(N-r*y2(3))*y2(1))-(mu+alpha*...
    (1+(rate1/N)*(t2>tD)*(t2<te)))*y2(2);...
    alpha*rate1/N*(t2>tD)*(y2(2)+y2(4))-(gammar+mu)*y2(3);...
    alpha*y2(2)-(gamma+mu+alpha*rate1/N)*y2(4);...
    gammar*y2(3)+gamma*y2(4)-mu*y2(5)-delta*d*y2(5)];

tspan = 0:1:150; % Time Span of the Simulation
conds = [960,0,0,40,0]; % Initial Conditons in order; 
                        % ('S(t)', 'E(t)', 'J(t)', 'I(t)', 'R(t)')
[t1,y1]=ode45(f, tspan, conds); 
[t2,y2]=ode45(g, tspan, conds);
plot(y2(:,:))
%hold on
%plot(y1(:,4),'g')
legend('S(t)', 'E(t)', 'J(t)', 'I(t)', 'R(t)')
Sum=0;
Sum2=0;
for i=1:75
sn=y1(i,3)+y1(i,4);
Sum=Sum+sn;
end
for i=1:75
sn=y2(i,3)+y2(i,4);
Sum2=Sum2+sn;
end

Infection_Toll = Sum2- Sum



