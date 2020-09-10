clear
dt=1/200
tf=1 %final time
nt=round((tf/dt),0)%or called h
nx=10
dx=1/nx %or k
x=linspace(0,1,nx)
s=1/2*dt/dx^2
t=linspace(0,tf,nt);

for i=1:nx
    if x(i)<=1/2
    U(i,1)=x(i)
    end
    if x(i)>1/2
    U(i,1)=1-x(i)
    end
end

A(1,1)=1-2*s;
A(1,2)=s;

for i=2:nx-1
    A(i,i)=1-2*s;
    A(i,i+1)=s;
    A(i,i-1)=s;
end
A(nx,nx)=1-2*s
for i=2:nt
    U(:,i)=A*U(:,i-1); %matrix multiplication could also use A^i to find the values...
    U(nx,i)=0;  % boundary conditions
    U(1,i)=0;   % boundary conditions
end

%to calculate exact values
syms tt
syms xx
v=0;



for i=1:45
   v=v+4/(i*pi)^2*sin(i*pi/2)*exp(-1/2*(i*pi)^2*tt)*sin(i*pi*xx);
end

g=matlabFunction(v);

for i=1:nt
    for j=1:nx
        UU(j,i)=g(t(i),x(j));
    end
end

for ii=1:nt
    Error(ii)=sqrt(trapz(abs(UU(:,ii)-U(:,ii)).^2));
end
Error=transpose(Error);


O=dx^2+dt
Ov=linspace(O,O,nx);
plot(x,UU(:,4),x,U(:,4))
%[XX YY]=meshgrid(t,x);
%figure(1)
%surf(XX,YY,UU)
%hold on
%plot3(XX,YY,U,'o')
%hold off
%grid on
