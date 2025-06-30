Nsteps=20;
dy=1/Nsteps;
%initial conditions
u=[];
for i=1:Nsteps
    u(i)=0;
end
u(21)=1;
E=1;
Red=5000;%Reynold's number
dt=E*Red*(dy)^2;
A=-E/2;
B=1+E;
K=[];
d=[];
a=[];
b=[];
K_=[];
K_(2)=0;
d_=[];
d_(2)=B;
ts=0;%time steps
while ts<240
    %forming equations
    for j=2:Nsteps
       K(j)=(1-E)*u(j)+(E/2)*(u(j+1)+u(j-1));
       d(j)=B;
    end
    for j=2:Nsteps-1
       a(j)=A;
    end
    for j=3:Nsteps
       b(j)=A;
    end
    K(Nsteps)=K(Nsteps)-A;
    
    %solving using Thomas's Algorithm
    %eliminating the lower diagonal elements
    for j=3:Nsteps
       d_(j)=d(j)-(b(j)*a(j-1))/d_(j-1);
       K_(j)=K(j)-(b(j)*K_(j-1))/d_(j-1);
    end
    u(Nsteps)=K_(Nsteps)/d_(Nsteps);
    %finding out u 
    for i=1:Nsteps-2
        u(Nsteps-i)=(K_(Nsteps-i)-a(Nsteps-i)*u(Nsteps-i+1))/d_(Nsteps-i);
    end
    ts=ts+1;
end
disp(u)
y=[];
y(1)=0;
for i=2:Nsteps+1
    y(i)=y(i-1)+dy;
end
plot(u,y);
xlabel('horizaontal velocity')
ylabel('distance ')
