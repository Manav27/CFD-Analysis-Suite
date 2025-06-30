A1=@(x) (-3.642e-8*x^4+1.68e-5*x^3-0.004032*x^2+0.5243*x+35.04)^2/35^2;%area contour of the nozzle

rho1=@(x) 1-0.3146*x;%initial conditions for density
T1=@(x) 1-0.23*x;%initial conditions for temperature
V1=@(x) (0.1+1.09*x)*sqrt(1-0.23*x);%initial conditions for velocity

g=1.4;
R=287;
L=3;%length of nozzle
Nsteps=30;%number of steps
dx=L/Nsteps;
dt=[];
rho=[];
T=[];
V=[];
A=[];
X=[];
ts=1400;%time steps

for i=1:Nsteps+1
    rho(i)=rho1((i-1)*dx);
    T(i)=T1((i-1)*dx);
    V(i)=V1((i-1)*dx);
    A(i)=A1((i-1)*dx);
end
for j=1:ts
    for i= 1:Nsteps+1
        dt(i)=0.5*dx/(sqrt(T(i))+V(i));
    end
        Dt=min(dt);
        %predictor step
    for i=1:Nsteps
        drdt(i) = -rho(i)*((V(i+1)-V(i))/dx)-rho(i)*V(i)*((log(A(i+1))-log(A(i)))/dx)-V(i)*((rho(i+1)-rho(i))/dx);
        dvdt(i) = -V(i)*((V(i+1)-V(i))/dx)-(1/g)*(((T(i+1)-T(i))/dx)+(T(i)/rho(i))*((rho(i+1)-rho(i))/dx));
        dTdt(i) = -V(i)*((T(i+1)-T(i))/dx)-(g-1)*T(i)*(((V(i+1)-V(i))/dx)+V(i)*((log(A(i+1))-log(A(i)))/dx));

        rho_(i)=rho(i)+drdt(i)*Dt;%predicted density
        V_(i)=V(i)+dvdt(i)*Dt;%predicted velocity
        T_(i)=T(i)+dTdt(i)*Dt;%predicted temperature
    end
        %corrector step
    for i=2:Nsteps
        drdt_(i)=-rho_(i)*((V_(i)-V_((i-1)))/dx)-rho_(i)*V_(i)*((log(A(i))-log(A(i-1)))/dx)-V_(i)*((rho_(i)-rho_(i-1))/dx);
        dvdt_(i)=-V_(i)*((V_(i)-V_(i-1))/dx)-(1/g)*(((T_(i)-T_(i-1))/dx)+(T_(i)/rho_(i))*((rho_(i)-rho_(i-1))/dx));
        dTdt_(i)=-V_(i)*((T_(i)-T_(i-1))/dx)-(g-1)*T_(i)*(((V_(i)-V_(i-1))/dx)+V_(i)*((log(A(i))-log(A(i-1)))/dx));

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
   
    ts=ts+1;
end

X(1)=0;
for i=2:Nsteps+1
    X(i)=X(i-1)+dx;
end
for i=1:Nsteps+1
    M(i)=V(i)/(sqrt(T(i)));
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