T=input('Enter the required thrust:')
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
u=1.2;
Tr=(Pr)^((u-1)/u);
TT=(2*u*R*Tc)/(u-1);
Ue=sqrt(TT*(1-Tr))
Te=Tc*Tr;
Ae=sqrt(u*R*Te)
Me=Ue/Ae
%Me = 3.28;
%Prandtl Meyer Function 
A=sqrt((u+1)/(u-1));
B=(u-1)/(u+1);
vpm=@(x) A*atan(sqrt(B*(x^2-1)))-atan(sqrt(x^2-1));

dtor=pi/180;
rtod=180/pi;

Tmax=0.5*vpm(Me)*rtod;
cp=[];%centreline points

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
p(1)=[];

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
    
    
    
    
    
    
    
    
    
    
    
 