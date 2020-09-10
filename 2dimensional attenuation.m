clear
q=65 %sum constant VERY EXPENSIVE CAREFUL!
nx=20 %x steps
nt=50 %time steps
a=-10 %left endpoint
b=10 %right endpoint
t0=0.01 %time start
tf=100 %time end
Ax=linspace(a,b,nx);
At=linspace(t0,tf,nt);
c=2;
v=0.0009/3;
f=0;
g=@(k) 15*sech(sqrt(20)/2*((k-3.6))).*sech(sqrt(20)/2*((k-1.6))); %initial condition
%g=@(x) e^(-x^2);  %initial condition
syms kk
syms k
syms t
r=exp(-(k-kk)^2/(4*pi*t))*exp(-kk^2);
rr=int(r,kk,-50,50);
rrr=1/2*(exp(-(k-c*t)^2)+exp(-(k+c*t)^2))+(4*v)/(6*sqrt(4*pi*t))*rr;
for i=1:nt
    for j=1:nx
        k=Ax(j);
        t=At(i);
        P(j,i)=(subs(rrr));
    end
end

U=double(P);

%
%for i=1:nt
    %for i=1:nx
        %Ploss(i,j)=4*v/(6*sqrt(4*pi*At(i)).*integral(r,At(i)),-500,500));
    %end
%end

Code for checking solution

clear
syms x
syms t
syms v
syms lambda
syms l
syms c
alpha=lambda^2*16*v^2-4*lambda*c;
syms l;
syms pi;w
lambda=pi^2/l^2;
f=((2*lambda/(l^2*sqrt(-alpha))*sin(sqrt(-alpha)/2*t)+cos(sqrt(-alpha)/2*t))*exp(-lambda*4*v*t))*sin(sqrt(lambda)*x);

r=-diff(f,t,2)+c*diff(f,x,2)+4/3*v*diff(diff(f,x,2),t);

%
x=55
t=1
v=25
lambda=1
l=123
c=1
b=double(subs(r));
