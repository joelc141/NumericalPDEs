%Models pressure wave attenuation due to thermoviscous effects, you give it
%an initial displacement and it will find how it losing its shape in time.
clear
q=65 %sum constant Keep low!!
nx=200 %x steps
nt=500 %time steps
a=0 %left endpoint
b=10 %right endpoint
t0=0 %time start
tf=25 %time end


Ax=linspace(a,b,nx);
At=linspace(t0,tf,nt);
c=0.45;
v=0.0009/3;
f=0;
g=@(k) 15*sech(sqrt(20)/2*((k-3.6))).*sech(sqrt(20)/2*((k-1.6))); %initial condition
%g=@(x) e^(-x^2);  %initial condition

for i=1:q
    fun= @(k) sqrt(c)*sech(sqrt(c)/2*((k-8.6))).*sech(sqrt(c)/2*((k-8.6))).*sin(i*pi*k/(b-a));
    qq = integral(fun,a,b);
    an(i)=2*qq;
    bn(i)=8/((2*i-1)^3*pi^3);
end
    
R=zeros(nx,nt);

for xx=1:nx
    for tt=1:nt
        x=Ax(xx);
        t=At(tt);
        f=linspace(0,0,q);
        for n=1:q
            lambda=(n)^2*pi^2/(b-a)^2;
            alpha=lambda.^2*16*v^2-4*lambda*c^2;
            f(n)=an(n)*((8*n^2*pi^2*v)/((b-a)^2*sqrt(-alpha))*sin(sqrt(-alpha)/2*t)+cos(sqrt(-alpha)/2*t))*sin(sqrt(lambda)*x)*exp(-lambda*4*v*t);
        end
        R(xx,tt)=sum(f);
    end 
end

for xx=1:nx
    for tt=1:nt
        v=0;
        x=Ax(xx);
        t=At(tt);
        f=linspace(0,0,q);
        for n=1:q
            lambda=(n)^2*pi^2/(b-a)^2;
            alpha=lambda.^2*16*v^2-4*lambda*c^2;
            f(n)=an(n)*((8*n^2*pi^2*v)/((b-a)^2*sqrt(-alpha))*sin(sqrt(-alpha)/2*t)+cos(sqrt(-alpha)/2*t))*sin(sqrt(lambda)*x)*exp(-lambda*4*v*t);
        end
        RR(xx,tt)=sum(f);
    end 
end

RRR=abs(R-RR); %energy loss term

for i=1:nt
    AR(i)=sum(RRR(:,i));
end

%for i=2:nx-1
    %for ii=2:nt-1
        %RRR(i,ii)=(RRR(i,ii+1)-RRR(i,ii-1))/(2*(tf-t0));
    %end
%end

v = VideoWriter('oil');
open(v);


for kk = 1:nt
   subplot(2,1,1);
   plot(Ax,transpose(R(:,kk)),Ax,transpose(RR(:,kk))) %,Ax,transpose(RR(:,kk)))
   title 'Oil pressure wave room temp'
   xlabel 'position(meters)'
   ylabel 'pressure(pascals)'
   legend('Viscosity', 'No viscosity')
   xlim([a b])
   ylim([-45 45])
   subplot(2,1,2);
   plot(Ax,transpose(abs(RR(:,kk)-R(:,kk))))
   ylabel 'Pressure loss'
   xlabel 'position(meters)'
   xlim([a b])
   ylim([0 max(max(abs(RR-R)))])
   frame = getframe(gcf);
   writeVideo(v,frame);
end
close(v);
