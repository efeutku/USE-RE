
N=1500;
r=0.6;
l=0.5;
mu=4.98*10^-5;
alpha=1/4;
zeta=1/3;
gamma=1/6;
gammar=1/5;
mean_cont=0.9575;
transm_risk=0.6;
beta=mean_cont*transm_risk; %R0
delta=0.15;
kappa=1/7;
kappar=1/9;
tau=1/10;
h=0.29;
d=1/14;
T=0;
T1=850;
ss=0.92; % Succes Rate of Sampling
st=0.95; % Succes Rate of Test
c=ss*st;
r1=(1-r)*l;
rate=T*c; 
rate1=T1*c;
tD1=20; % Start of Drone Operation
te1=70; % End of Drone Operation
tD2=20; % Start of Drone Operation
te2=70; % End of Drone Operation


tspan = 0:1:800;
xinit = [1495,0,0,5,0,0];
options = odeset('NonNegative',4);

f = @(t1,y1) [mu*N-mu*y1(1)-beta*(y1(4)+r1*y1(3))/(N-r*y1(3))*y1(1)+delta*d*y1(5);...
        (beta*(y1(4)+r1*y1(3))/(N-r*y1(3))*y1(1))-(mu+alpha+(rate/N)*(t1>tD1)*(t1<te1))...
        *y1(2);...
        rate/N*(t1>tD1)*(t1<te1)*(y1(2))+zeta*rate/N*(t1>tD1)*(t1<te1)*y1(4)-(gammar+mu+kappar*h)*y1(3);...
        alpha*y1(2)-(kappa*h+gamma+mu+zeta*rate/N*(t1>tD1)*(t1<te1))*y1(4);...
        tau*y1(6)+gammar*y1(3)+gamma*y1(4)-mu*y1(5)-delta*d*y1(5);
        h*(kappar*y1(3)+kappa*y1(4))-(mu+tau)*y1(6)];
    [t1,y1]=ode45(f, tspan, xinit,options);
    
    
g = @(t2,y2) [mu*N-mu*y2(1)-beta*(y2(4)+r1*y2(3))/(N-r*y2(3))*y2(1)+delta*d*y2(5);...
        (beta*(y2(4)+r1*y2(3))/(N-r*y2(3))*y2(1))-(mu+alpha+(alpha*rate1/N)*(t2>tD2)*(t2<te2))...
        *y2(2);...
        alpha*rate1/N*(t2>tD2)*(t2<te2)*(y2(2))+zeta*rate1/N*(t2>tD2)*(t2<te2)*y2(4)-(gammar+mu+kappar*h)*y2(3);...
        alpha*y2(2)-(kappa*h+gamma+mu+zeta*rate1/N*(t2>tD2)*(t2<te2))*y2(4);...
        tau*y2(6)+gammar*y2(3)+gamma*y2(4)-mu*y2(5)-delta*d*y2(5);
        h*(kappar*y2(3)+kappa*y2(4))-(mu+tau)*y2(6)];
    [t2,y2]=ode45(g, tspan, xinit,options);
   
    
    
I1=y1(:,3)+y1(:,4)+y1(:,6);
I2=y2(:,3)+y2(:,4)+y2(:,6);
 
HC=ones(size(I1))*4;
  

 plot(y1(:,6),'r')
 legend('S(t)', 'E(t)', 'J(t)', 'I(t)', 'R(t)')
 hold on
 plot(y2(:,6),'b')
 %plot(HC)
Sum1=0;
Sum2=0;
for i=1:100
sn=I1(i);
Sum1=Sum1+sn;
end
for i=1:100
sn=I2(i);
Sum2=Sum2+sn;
end
Sum2-Sum1

