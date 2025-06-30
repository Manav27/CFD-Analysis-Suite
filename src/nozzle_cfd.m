clear all
%CONSTRUCTING THE NOZZLE

Thrust=input('Enter the required thrust:')
Pc=input('Enter chamber pressure:')
Tc=input('Enter chamber temperature:')
TR=input('Enter Throat Radius: ')
R=input('Enter universal gas constant: ') 
Alt=input('Enter altitude in meters: ')

if Alt<11000
    Pe=1000*(101.92*((288.14-0.00649*Alt)/(288.08))^5.256);
else if 11000<Alt<25000
        Pe=1000*(22.65*exp(1.73-0.000157*Alt));
else if Alt>25000
            Pe=1000*(2.488*((141.89+0.00299*Alt)/216.6)^(-11.388));
    end
    end
end
Pr=Pe/Pc;
u=1.4;
Tr=(Pr)^((u-1)/u);
TT=(2*u*R*Tc)/(u-1);
Ue=sqrt(TT*(1-Tr));
Te=Tc*Tr;
Ae=sqrt(u*R*Te);
Me=Ue/Ae;
%Prandtl Meyer Function 
A=sqrt((u+1)/(u-1));
B=(u-1)/(u+1);
vpm=@(x) A*atan(sqrt(B*(x^2-1)))-atan(sqrt(x^2-1));

dtor=pi/180;
rtod=180/pi;

Tmax=0.5*vpm(Me)*rtod;
cp=[];%centerline points

dt=90-Tmax -fix(90-Tmax);
t(1)=dt*dtor;
n=2*Tmax;

for i=2:n+1
 t(i)=(dt+(i-1))*dtor;
 func=@(x) t(i)-vpm(x);
 M(i)=fzero(func,1.01*Me);
 p(i)=TR*tan(t(i));
 rr(i)=-TR/p(i);
 lr(i)=tan(t(i)+ asin(1/M(i)));
 sl(i)=-rr(i);
end

lr(1)=[];
rr(1)=[];
sl(1)=[];
p(1)=[]

for j=1:length(p)
    P1=[p(j) 0];
    P2=[0 TR];
    plot(P1,P2,'b')
    xlabel('Centreline');
    ylabel('Throat');
end
hold on;
F=rr(i-1);
for k=1:length(p)-1
    x(k)=(TR+sl(k)*p(k))/(sl(k)-F);
    y(k)=F*x(k)+TR;
    rx=[p(k) x(k)];
    ry=[0 y(k)];
    plot(rx,ry,'b')
end
hold on;
Tm=Tmax*dtor;
s(1)=tan(Tm);
 xw(1)=(TR+sl(1)*p(1))/(sl(1)-s(1));
 yw(1)=s(1)*xw(1) +TR;
 XWS1=[p(1) xw(1)];
 YWS1=[p(2) yw(1)];
 plot(XWS1,YWS1,'g')
 
 dtw=Tm/(length(p)-1);
 b(1)=TR;
 for l =2:length(p)-1
     s(l)=tan(Tm-(l-1)*dtw);
     b(l)=yw(l-1)-s(l)*xw(l-1);
     xw(l)=(b(l)+sl(l)*p(l))/(sl(l)-s(l));
     yw(l)=s(l)*xw(l)+b(l);
     XW=[x(l) xw(l)];
     YW=[y(l) yw(l)];
     plot(XW,YW,'r')
 end
 hold on;
    
xf=(b(length(b))+sl(length(sl))*p(length(p)))/sl(length(sl));
yf=b(length(b));
Xf=[p(length(p)) xf];
Yf=[0 yf];
plot(Xf,Yf,'b')

%NOZZLE CFD

rho1=@(x) 1-0.00203886*x;%initial conditions for density
T1=@(x) 1-0.00217734*x;%initial conditions for temperature
V1=@(x) (0.1+1.09*x)*sqrt(1-0.00217734*x);%initial conditions for velocity

g=1.4;
Lc=0.4*xw(end);%length of converging nozzle
L=xw(end)+Lc;%length of nozzle
Nsteps=100;%number of steps
dx=L/Nsteps;
dt_=[];
rho=[];
T=[];
V=[];
Area=[];
X=[];
ts=3000;%time steps


y1=@(x) yw(end)+15-((yw(end)+15-TR)/Lc)*x;
A2= fit(xw',yw','poly5');%Fitted curve of the nozzle

A1=@(x) (yw(end)+15-((yw(end)+15-TR)/Lc)*x)^2/TR^2;



for i=1:floor(Lc/dx)
    Area(i)=A1((i-1)*dx);
    wall(i)=y1((i-1)*dx);
end
for i=floor(Lc/dx)+1:Nsteps+1
    Area(i)=(A2((i-1-floor(Lc/dx))*dx))^2/TR^2;
    wall(i)=A2((i-1-floor(Lc/dx))*dx);
end
for i=1:Nsteps+1
    rho(i)=rho1((i-1)*dx);
    T(i)=T1((i-1)*dx);
    V(i)=V1((i-1)*dx);
end
for j=1:ts
    Ts(j)=j;
    for i= 1:Nsteps+1
        dt_(i)=0.5*dx/(sqrt(T(i))+V(i));
    end
        Dt=min(dt_);
    %predictor step
    for i=1:Nsteps
        drdt(i) = -rho(i)*((V(i+1)-V(i))/dx)-rho(i)*V(i)*((log(Area(i+1))-log(Area(i)))/dx)-V(i)*((rho(i+1)-rho(i))/dx);
        dvdt(i) = -V(i)*((V(i+1)-V(i))/dx)-(1/g)*(((T(i+1)-T(i))/dx)+(T(i)/rho(i))*((rho(i+1)-rho(i))/dx));
        dTdt(i) = -V(i)*((T(i+1)-T(i))/dx)-(g-1)*T(i)*(((V(i+1)-V(i))/dx)+V(i)*((log(Area(i+1))-log(Area(i)))/dx));

        rho_(i)=rho(i)+drdt(i)*Dt;%predicted density
        V_(i)=V(i)+dvdt(i)*Dt;%predicted velocity
        T_(i)=T(i)+dTdt(i)*Dt;%predicted temperature
    end
        %corrector step
    for i=2:Nsteps
        drdt_(i)=-rho_(i)*((V_(i)-V_((i-1)))/dx)-rho_(i)*V_(i)*((log(Area(i))-log(Area(i-1)))/dx)-V_(i)*((rho_(i)-rho_(i-1))/dx);
        dvdt_(i)=-V_(i)*((V_(i)-V_(i-1))/dx)-(1/g)*(((T_(i)-T_(i-1))/dx)+(T_(i)/rho_(i))*((rho_(i)-rho_(i-1))/dx));
        dTdt_(i)=-V_(i)*((T_(i)-T_(i-1))/dx)-(g-1)*T_(i)*(((V_(i)-V_(i-1))/dx)+V_(i)*((log(Area(i))-log(Area(i-1)))/dx));

        drdtavg(i)=0.5*(drdt(i)+drdt_(i));
        dvdtavg(i)=0.5*(dvdt(i)+dvdt_(i));
        dTdtavg(i)=0.5*(dTdt(i)+dTdt_(i));

        rho(i)=rho(i)+drdtavg(i)*Dt;
        V(i)=V(i)+dvdtavg(i)*Dt;
        T(i)=T(i)+dTdtavg(i)*Dt;
    end    
        %boundary values
    V(1)=2*V(2)-V(3);
    V(Nsteps+1)=2*V(Nsteps)-V(Nsteps-1);
    rho(Nsteps+1)=2*rho(Nsteps)-rho(Nsteps-1);
    T(Nsteps+1)=2*T(Nsteps)-T(Nsteps-1);
   
    Vts(j)=V(Nsteps)+1;
    rhots(j)=rho(Nsteps)+1;
    Tts(j)=T(Nsteps)+1;
    ts=ts+1;
end

X(1)=0;
for i=2:Nsteps+1
    X(i)=X(i-1)+dx;
end
for i=1:Nsteps+1
    Ma(i)=V(i)/(sqrt(T(i)));
    P(i)=rho(i)*T(i);
end

disp('rho=')
disp(rho)
disp('T=')
disp(T)
disp('V=')
disp(V)
disp('P=')
disp(P)
figure();
plot(X,wall)
figure();
plot(X,Ma)
xlabel('distance')
ylabel('Mach number')
figure();
plot(X,rho)
xlabel('distance')
ylabel('density')
figure();
plot(X,T)
xlabel('distance')
ylabel('Temperature')
figure();
plot(X,V)
xlabel('distance')
ylabel('Velocity')
figure();
plot(X,P)
xlabel('distance')
ylabel('Pressure')
figure();
plot(Ts,Vts)
xlabel('time steps')
ylabel('Velocity')