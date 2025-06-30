Imin=1;
Imax=70;
Jmin=1;
Jmax=70;
M_fs=4;%mach number
LHORI=0.00001;%plate length in meters
a_fs=340.28;%free stream speed of sound
T_fs=288.16;%free stream tempertaure
P_fs=101325;%free stream pressure
Tw=T_fs;%wall temperature
g=1.4;%ratio of specific heats
Pr=0.71;%prandtl number
R=287;%universal gas constant
mu0=0.000017894;%refrence dynamic viscosity
ts=0;%time step

Cv=R/(g-1);
Cp=g*Cv;
rho_fs=P_fs/(R*T_fs);%free stream density
u_fs=M_fs*sqrt(g*R*T_fs);%free stream velocity
Rel=rho_fs*u_fs*LHORI/mu0;%Reynolds number for the plate
k_fs=mu0*Cp/Pr;%conductivity constant
e_fs=Cv*T_fs;%internal energy
E_fs=rho_fs*(e_fs+((u_fs)^2)/2);

del=5*LHORI/sqrt(Rel);%boundary layer thickness
LVERT=5*del;

dx=LHORI/(Imax-1);%step size in x direction
dy=LVERT/(Jmax-1);%step size in y direction

%initializing arrays
P=[];
u=[];
v=[];
rho=[];
T=[];
E=[];
e=[];
mu=[];
k=[];
V=[];
dt=[];
a=[];
U1=[];
U2=[];
U3=[];
U4=[];
E1=[];
E2=[];
E3=[];
E4=[];
F1=[];
F2=[];
F3=[];
F4=[];
TAUXX=[];
TAUXY=[];
TAUYY=[];
QX=[];
QY=[];
U1_=[];
U2_=[];
U3_=[];
U4_=[];
E1_=[];
E2_=[];
E3_=[];
E4_=[];
F1_=[];
F2_=[];
F3_=[];
F4_=[];

%initial conditions
%for points except on plate
for i=1:Imax
    for j=2:Jmax
        u(i,j)=u_fs;
        v(i,j)=0;
        T(i,j)=T_fs;
    end
end
%for points on plate
for i=1:Imax
    u(i,1)=0;
    v(i,1)=0;
    T(i,j)=Tw;
end
%for entire grid
for i=1:Imax
    for j=1:Jmax
        rho(i,j)=rho_fs;
        mu(i,j)=mu0;
        k(i,j)=k_fs;
        a(i,j)=sqrt((u(i,j))^2+(v(i,j))^2)*sqrt(g*R*T(i,j));
        e(i,j)=e_fs;
        E(i,j)=rho(i,j)*(e(i,j)+((u(i,j))^2+(v(i,j))^2)/2);
        P(i,j)=P_fs;
    end
end
while ts<4400
    
    %determining time step
    for i=1:Imax
        for j=1:Jmax
            V(i,j)=(4/3)*mu(i,j)*(g*mu(i,j))/(Pr*rho(i,j));
        end
    end
    V_=max(V);
    V__=max(V_);
    for i=1:Imax
        for j=1:Jmax
            dt(i,j)=0.7*(u(i,j)/dx+v(i,j)/dy+a(i,j)*sqrt((1/(dx)^2)+(1/(dy)^2))+2*V__*((1/(dx)^2)+(1/(dy)^2))).^(-1);
        end
    end
    Dt=min(dt);
    Dt_=min(Dt);
    
    %for finding shear stress and thermal conductivity
    %for internal points
    for i=2:Imax-1
        for j=2:Jmax-1
            TAUXY(i,j)=mu(i,j)*((u(i,j+1)-u(i,j-1))/(2*dy)+(v(i+1,j)-v(i-1,j))/(2*dx));
            TAUXX(i,j)=-P(i,j)+2*mu(i,j)*(u(i+1,j)-u(i-1,j))/(2*dx);
            TAUYY(i,j)=-P(i,j)+2*mu(i,j)*(v(i,j+1)-v(i,j-1))/(2*dy);
            QX(i,j)=-k(i,j)*(T(i+1,j)-T(i-1,j))/(2*dx);
            QY(i,j)=-k(i,j)*(T(i,j+1)-T(i,j-1))/(2*dy);
        end
    end
    %for points on plate
    for i=2:Imax
        TAUXY(i,1)=mu(i,1)*((u(i,2)-u(i,1))/dy+(v(i,1)-v(i-1,1))/dx);
        TAUXX(i,1)=-P(i,1)+2*mu(i,1)*(u(i,1)-u(i-1,1))/(dx);
        TAUYY(i,1)=-P(i,1)+2*mu(i,1)*(v(i,2)-v(i,1))/(dy);
        QX(i,1)=-k(i,1)*(T(i,1)-T(i-1,1))/(dx);
        QY(i,1)=-k(i,1)*(T(i,2)-T(i,1))/(dy);
    end
    %for upper boundary
    for i=1:Imax-1
        TAUXY(i,Jmax)=mu(i,Jmax)*((u(i,Jmax)-u(i,Jmax-1))/(dy)+(v(i+1,Jmax)-v(i,Jmax))/(dx));
        TAUXX(i,Jmax)=-P(i,Jmax)+2*mu(i,Jmax)*(u(i+1,Jmax)-u(i,Jmax))/(dx);
        TAUYY(i,Jmax)=-P(i,Jmax)+2*mu(i,Jmax)*(v(i,Jmax)-v(i,Jmax-1))/(dy);
        QX(i,Jmax)=-k(i,Jmax)*(T(i+1,Jmax)-T(i,Jmax))/(dx);
        QY(i,Jmax)=-k(i,Jmax)*(T(i,Jmax)-T(i,Jmax-1))/(dy);
    end
    %for inflow
    for j=1:Jmax-1
        TAUXY(1,j)=mu(1,j)*((u(1,j+1)-u(1,j))/(dy)+(v(2,j)-v(1,j))/(dx));
        TAUXX(1,j)=-P(1,j)+2*mu(1,j)*(u(1,j)-u(2,j))/(dx);
        TAUYY(1,j)=-P(1,j)+2*mu(1,j)*(v(1,j+1)-v(1,j))/(dy);
        QX(1,j)=-k(1,j)*(T(2,j)-T(1,j))/(dx);
        QY(1,j)=-k(1,j)*(T(1,j+1)-T(1,j))/(dy);
    end
    %for outflow
    for j=2:Jmax
        TAUXY(Imax,j)=mu(Imax,j)*((u(Imax,j)-u(Imax,j-1))/(dy)+(v(Imax,j)-v(Imax-1,j))/(dx));
        TAUXX(Imax,j)=-P(Imax,j)+2*mu(Imax,j)*(u(Imax,j)-u(Imax-1,j))/(dx);
        TAUYY(Imax,j)=-P(Imax,j)+2*mu(Imax,j)*(v(Imax,j)-v(Imax,j-1))/(dy);
        QX(Imax,j)=-k(Imax,j)*(T(Imax,j)-T(Imax-1,j))/(dx);
        QY(Imax,j)=-k(Imax,j)*(T(Imax,j)-T(Imax,j-1))/(dy);
    end

    for i=1:Imax
        for j=1:Jmax
            U1(i,j)=rho(i,j);
            U2(i,j)=rho(i,j)*u(i,j);
            U3(i,j)=rho(i,j)*v(i,j);
            U4(i,j)=E(i,j);
            E1(i,j)=rho(i,j)*u(i,j);
            E2(i,j)=rho(i,j)*(u(i,j))^2+P(i,j)-TAUXX(i,j);
            E3(i,j)=rho(i,j)*u(i,j)*v(i,j)-TAUXY(i,j);
            E4(i,j)=(E(i,j)+P(i,j))*u(i,j)-u(i,j)*TAUXX(i,j)-v(i,j)*TAUXY(i,j)+QX(i,j);
            F1(i,j)=rho(i,j)*v(i,j);
            F2(i,j)=rho(i,j)*u(i,j)*v(i,j)-TAUXY(i,j);
            F3(i,j)=rho(i,j)*(v(i,j))^2+P(i,j)-TAUXY(i,j);
            F4(i,j)=(E(i,j)+P(i,j))*v(i,j)-u(i,j)*TAUXY(i,j)-v(i,j)*TAUYY(i,j)+QY(i,j);
        end
    end

    %Mac Cormack's Technique
    %predictor step
    for i=1:Imax-1
        for j=1:Jmax-1
            dU1dt(i,j)=-(E1(i+1,j)-E1(i,j))/dx-(F1(i,j+1)-F1(i,j))/dy;
            dU2dt(i,j)=-(E2(i+1,j)-E2(i,j))/dx-(F2(i,j+1)-F2(i,j))/dy;
            dU3dt(i,j)=-(E3(i+1,j)-E3(i,j))/dx-(F3(i,j+1)-F3(i,j))/dy;
            dU4dt(i,j)=-(E4(i+1,j)-E4(i,j))/dx-(F4(i,j+1)-F4(i,j))/dy;
            %predicted values
            U1_(i,j) = U1(i,j)+ dU1dt(i,j)*Dt_;
            U2_(i,j) = U2(i,j)+ dU2dt(i,j)*Dt_;
            U3_(i,j) = U3(i,j)+ dU3dt(i,j)*Dt_;
            U4_(i,j) = U4(i,j)+ dU4dt(i,j)*Dt_;
        end
    end
    %finding out primitive variables
    for i=1:Imax-1
        for j=1:Jmax-1
            rho_(i,j)=U1_(i,j);
            u_(i,j)=U2_(i,j)/U1_(i,j);
            v_(i,j)=U3_(i,j)/U1_(i,j);
            E_(i,j)=U4_(i,j);
            e_(i,j)=U4_(i,j)/U1_(i,j)-((u_(i,j))^2+(v_(i,j))^2)/2;
            T_(i,j)=e_(i,j)/Cv;
            P_(i,j)=rho_(i,j)*R*T_(i,j);
            mu_(i,j)=mu0*((T_(i,j)/T_fs)^1.5)*(T_fs+110)/(T_(i,j)+110);
            k_(i,j)=mu_(i,j)*Cp*T_(i,j);
        end
    end
    
    %boundary conditions
    %upper boundsary
    for i=1:Imax
        rho_(i,Jmax)=rho_fs;
        u_(i,Jmax)=u_fs;
        v_(i,Jmax)=0;
        P_(i,Jmax)=P_fs;
        T_(i,Jmax)=T_fs;
        E_(i,Jmax)=E_fs;
        e_(i,Jmax)=e_fs;
        mu_(i,Jmax)=mu0*(T_(i,Jmax)/T_fs)^1.5*(T_fs+110)/(T_(i,Jmax)+110);
        k_(i,Jmax)=mu_(i,Jmax)*Cp*T_(i,Jmax);
    end
    %outflow
   for j=1:Jmax
        u_(Imax,j)=2*u_(Imax-1,j)-u_(Imax-2,j);
        v_(Imax,j)=2*v_(Imax-1,j)-v_(Imax-2,j);
        T_(Imax,j)=2*T_(Imax-1,j)-T_(Imax-2,j);
        P_(Imax,j)=2*P_(Imax-1,j)-P_(Imax-2,j);
        rho_(Imax,j)=2*rho_(Imax-1,j)-rho_(Imax-2,j);
        E_(Imax,j)=2*E_(Imax-1,j)-E_(Imax-2,j);
        e_(Imax,j)=2*e_(Imax-1,j)-e_(Imax-2,j);
        mu_(Imax,j)=mu0*(T_(Imax,j)/T_fs)^1.5*(T_fs+110)/(T_(Imax,j)+110);
        k_(Imax,j)=mu_(Imax,j)*Cp*T_(Imax,j);
    end
    
    %at internal points
    for i=2:Imax-1
        for j=2:Jmax-1
            TAUXY_(i,j)=mu_(i,j)*((u_(i,j+1)-u_(i,j-1))/(2*dy)+(v_(i+1,j)-v_(i-1,j))/(2*dx));
            TAUXX_(i,j)=-P_(i,j)+2*mu_(i,j)*(u_(i+1,j)-u_(i-1,j))/(2*dx);
            TAUYY_(i,j)=-P_(i,j)+2*mu_(i,j)*(v_(i,j+1)-v_(i,j-1))/(2*dy);
            QX_(i,j)=-k_(i,j)*(T_(i+1,j)-T_(i-1,j))/(2*dx);
            QY_(i,j)=-k_(i,j)*(T_(i,j+1)-T_(i,j-1))/(2*dy);
        end
    end
    %at inflow
    for j=2:Jmax
         TAUXY_(1,j)=mu_(1,j)*((u_(1,j)-u_(1,j-1))/(dy)+(v_(2,j)-v_(1,j))/(dx));
         TAUXX_(1,j)=-P_(1,j)+2*mu_(1,j)*(u_(2,j)-u_(1,j))/(dx);
         TAUYY_(1,j)=-P_(1,j)+2*mu_(1,j)*(v_(1,j)-v_(1,j-1))/(dy);
         QX_(1,j)=-k_(1,j)*(T_(2,j)-T_(1,j))/(dx);
         QY_(1,j)=-k_(1,j)*(T_(1,j)-T_(1,j-1))/(dy);
    end
    %on the plate
    for i=1:Imax-1
         TAUXY_(i,1)=mu_(i,1)*((u_(i,2)-u_(i,1))/(dy)+(v_(i+1,1)-v_(i,1))/(dx));
         TAUXX_(i,1)=-P_(i,1)+2*mu_(i,1)*(u_(i+1,1)-u_(i,1))/(dx);
         TAUYY_(i,1)=-P_(i,1)+2*mu_(i,1)*(v_(i,2)-v_(i,1))/(dy);
         QX_(i,1)=-k_(1,j)*(T_(i+1,1)-T_(i,1))/(dx);
         QY_(i,1)=-k_(1,j)*(T_(i,2)-T_(i,1))/(dy);
    end
    %on the upper boundary
    for j=1:Jmax-1
            TAUXY_(Imax,j)=mu_(Imax,j)*((u_(Imax,j+1)-u_(Imax,j))/(dy)+(v_(Imax,j)-v_(Imax-1,j))/(dx));
            TAUXX_(Imax,j)=-P_(Imax,j)+2*mu_(Imax,j)*(u_(Imax,j)-u_(Imax-1,j))/(dx);
            TAUYY_(Imax,j)=-P_(Imax,j)+2*mu_(Imax,j)*(v_(Imax,j+1)-v_(Imax,j))/(dy);
            QX_(Imax,j)=-k_(Imax,j)*(T_(Imax,j)-T_(Imax-1,j))/(dx);
            QY_(Imax,j)=-k_(Imax,j)*(T_(Imax,j+1)-T_(Imax,j))/(dy);        
    end
    %on the outflow
    for i=2:Imax-1
            TAUXY_(i,Jmax)=mu_(i,Jmax)*((u_(i,Jmax)-u_(i,Jmax-1))/(dy)+(v_(i,Jmax)-v_(i-1,Jmax))/(dx));
            TAUXX_(i,Jmax)=-P_(i,Jmax)+2*mu_(i,Jmax)*(u_(i,Jmax)-u_(i-1,Jmax))/(dx);
            TAUYY_(i,Jmax)=-P_(i,Jmax)+2*mu_(i,Jmax)*(v_(i,Jmax)-v_(i,Jmax-1))/(dy);
            QX_(i,Jmax)=-k_(i,Jmax)*(T_(i,Jmax)-T_(i-1,Jmax))/(dx);
            QY_(i,Jmax)=-k_(i,Jmax)*(T_(i,Jmax)-T_(i,Jmax-1))/(dy);
    end
    %predicted values
    for i=1:Imax
        for j=1:Jmax
            E1_(i,j)=rho_(i,j)*u_(i,j);
            E2_(i,j)=rho_(i,j)*(u_(i,j))^2+P_(i,j)-TAUXX_(i,j);
            E3_(i,j)=rho_(i,j)*u_(i,j)*v_(i,j)-TAUXY_(i,j);
            E4_(i,j)=(E_(i,j)+P_(i,j))*u_(i,j)-u_(i,j)*TAUXX_(i,j)-v_(i,j)*TAUXY_(i,j)+QX_(i,j);
            F1_(i,j)=rho_(i,j)*v_(i,j);
            F2_(i,j)=rho_(i,j)*u_(i,j)*v_(i,j)-TAUXY_(i,j);
            F3_(i,j)=rho_(i,j)*(v_(i,j))^2+P_(i,j)-TAUXY_(i,j);
            F4_(i,j)=(E_(i,j)+P_(i,j))*v_(i,j)-u_(i,j)*TAUXY_(i,j)-v_(i,j)*TAUYY_(i,j)+QY_(i,j);
        end
    end
    %corrector step
    for i=2:Imax-1
        for j=2:Jmax-1
            dU1dt_(i,j)=-(E1_(i,j)-E1_(i-1,j))/dx-(F1_(i,j)-F1_(i,j-1))/dy;
            dU2dt_(i,j)=-(E2_(i,j)-E2_(i-1,j))/dx-(F2_(i,j)-F2_(i,j-1))/dy;
            dU3dt_(i,j)=-(E3_(i,j)-E3_(i-1,j))/dx-(F3_(i,j)-F3_(i,j-1))/dy;
            dU4dt_(i,j)=-(E4_(i,j)-E4_(i-1,j))/dx-(F4_(i,j)-F4_(i,j-1))/dy;
        end
    end
    
    for i=2:Imax-1
        for j=2:Jmax-1
            dU1dtavg=0.5*(dU1dt(i,j)+dU1dt_(i,j));
            dU2dtavg=0.5*(dU2dt(i,j)+dU2dt_(i,j));
            dU3dtavg=0.5*(dU3dt(i,j)+dU3dt_(i,j));
            dU4dtavg=0.5*(dU4dt(i,j)+dU4dt_(i,j));

            U1(i,j)=U1(i,j)+dU1dtavg*Dt_;
            U2(i,j)=U2(i,j)+dU2dtavg*Dt_;
            U3(i,j)=U3(i,j)+dU3dtavg*Dt_;
            U4(i,j)=U4(i,j)+dU4dtavg*Dt_;
        end
    end
    %finding primitive variables
    for i=2:Imax-1
        for j=2:Jmax-1
            rho(i,j)=U1(i,j);
            u(i,j)=U2(i,j)/U1(i,j);
            v(i,j)=U3(i,j)/U1(i,j);
            E(i,j)=U4(i,j);
            e(i,j)=U4(i,j)/U1(i,j)-((u(i,j))^2+(v(i,j))^2)/2;
            T(i,j)=e(i,j)/Cv;
            P(i,j)=rho(i,j)*R*T(i,j);
            mu(i,j)=mu0*(T(i,j)/T_fs)^1.5*(T_fs+110)/(T(i,j)+110);
            k(i,j)=mu(i,j)*Cp*T_(i,j);
        end
    end
            
    %boundary conditions
    u(1,1)=0;
    v(1,1)=0;
    T(1,1)=T_fs;
    P(1,1)=P_fs;
    %upper boundsary
    for i=1:Imax
        u(i,Jmax)=u_fs;
        v(i,Jmax)=0;
        P(i,Jmax)=P_fs;
        T(i,Jmax)=T_fs;
    end
    %plate
    for i=1:Imax
        u(i,1)=0;
        v(i,1)=0;
        T(i,1)=Tw;
        P(i,1)=2*P(i,2)-P(i,3);
    end
    %outflow
    for j=1:Jmax
        u(Imax,j)=2*u(Imax-1,j)-u(Imax-2,j);
        v(Imax,j)=2*v(Imax-1,j)-v(Imax-2,j);
        T(Imax,j)=2*T(Imax-1,j)-T(Imax-2,j);
        P(Imax,j)=2*P(Imax-1,j)-P(Imax-2,j);
    end
    %inflow
    for j=2:Jmax-1
        u(1,j)=u_fs;
        v(1,j)=0;
        T(1,j)=T_fs;
        P(1,j)=P_fs;
    end
    ts=ts+1;
end

y(1)=0;
    for j=2:Jmax
        y(j)=y(j-1)+dy;
    end
  
plot(P(70,1:70)/P_fs,y/LVERT,'g')
hold on
plot(T(70,1:70)/T_fs,y/LVERT,'r')
hold on
plot(u(70,1:70)/u_fs,y/LVERT,'k')
hold on
xlabel('normalized Pressure,Temperature,Velocity')
ylabel('Normalized y distance')