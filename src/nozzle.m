Me = input('mach exit number')
g=1.2;
%method of characteristics
TR=57;
rtod=180/pi;
dtor=pi/180;
p=[];%x axis points

%Prandtl meyer function
A=sqrt((g+1)/(g-1));
B=(g-1)/(g+1);
vpm=@(x) A*atan(sqrt(B*(x^2-1)))-atan(sqrt(x^2-1));

%calculate tmax break up into duvision
tmax=0.5*vpm(Me)*rtod;
dt=(90-tmax)-fix(90-tmax);
t(1)=dt*dtor;
n=tmax*2-1;

for m=2:n+1
    t(m)=(dt +(m-1))*dtor;
    %mach number from t(i) using t(i)=vpm(false position)
    func=@(x) t(m)-vpm(x);
    M(m)=fzero(func,Me);
    p(m)=TR*tan(t(m));%x axis points
    %RR characteristics
    rr(m)=-TR/p(m);
    %LR cahracteristics
    lr(m)=tan(t(m)+asin(1/M(m)));
    sl(m)=-rr(m);
end

%plotting
p(1)=[]
l=length(p);
for j=1:l
    p1=[0 TR];
    p2=[p(j) 0];
    plot(p2,p1,'k')
    hold on
    xlabel('centreline')
    ylabel('throat')
end
hold on;
lr(1)=[];
rr(1)=[];
sl(1)=[];
F=rr(m-1);%slope of the last right running wave

for c=1:l-1
    x(c)=(TR+sl(c)*p(c))/(sl(c)-F);%x coordinate of intersection point of 
                                   %left running wave at point c and last
                                   %right running wave
                                   
    %solve:
    %y=F*x+TR and y=sl(c)+y-intercept(y-intercept can be found by
    %substituting x-intercept i.e.p(c) in the equation
    
    y(c)=F*x(c)+TR;%y coordinate of intersection point of left running wave
                   %at point c and last right running wave
    % solved by substituting x(c) in equation y=sl(c)*x+TR 
    
    xp=[p(c) x(c)];
    yp=[0 y(c)];
    plot(xp,yp,'b')
end
hold on   
%first wall section
tm=tmax*dtor;
xw(1)=(TR+sl(1)*p(1))/(sl(1)-tan(tm));%x coordinate of intersection point of
                                      %first left running wave and wall
%solve:
%y=tan(tmax)*x+TR and y=sl(1)*x + y-intercept(y-intercept can be found by
%substituting x-intercept in it i.e.p(1)

yw(1)=tan(tm)*xw(1)+TR;%y-coordinate of intersection point of first left 
                       %running wave and wall
% solve by substituting xw(1) in y=tan(ymax)+TR

xp2=[p(1) xw];
yp2=[p(2) yw];
%plot(xp2,yp2,'g');
%divide delta slopes
dtw=tan(tm)/(length(p)-1);%amount by which slope of wall decreaces along 
                          %the length 
s(1)=tan(tm); %slope of the first wall section    
b(1)=TR;%y-intercept of the first wall section

for k=2:length(p)-1
    s(k)=tan(tm)-(k-1)*dtw;%slope of kth wall section
    b(k)=yw(k-1)-s(k)*xw(k-1);%y int of kth wall section
    
    xw(k)=(b(k)+sl(k)*p(k))/(sl(k)-s(k));%x coordinate of intersection 
                                         %point of kth left running wave
                                         %and kth wall section or end
                                         %point of the kth wall sedction
    %solve:
    %y=s(k)*x+b(k) and y=sl(k)*x +y-intercept(y-intercept can be found by
    %substituting x-intercept in it i.e.p(1))
    
    yw(k)=s(k)*xw(k)+b(k);%y coordinate of intersection point of kth 
                          %left running wave and kth wall section or end
                          %point of the kth wall sedction
    %solve by substituting xw(k) in y=s(k)*x +b(k)
    
    xp3=[x(k) xw(k)];
    yp3=[y(k) yw(k)];
    plot(xp3,yp3,'r');
end
hold on
%last point
xf=(b(length(b))+sl(length(sl))*p(length(p)))/sl(length(sl))-0;
% slope of the last wall section should be zero 

yf=b(length(b));    
    
Xf=[p(length(p)) xf];
Yf=[0 yf];
plot(Xf,Yf,'b')    

xw=[0 xw];
yw=[TR yw];

T=table(transpose(xw),transpose(yw));

writetable(T,'points.txt');


    
    
    
    