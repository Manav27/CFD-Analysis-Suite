theta_c=30;%half cone angle 
M_fs=2.5;%Mach Number
g=1.4;%ratio og specific heats
theta_s=theta_c+12.48;%initial shock angle
rtod=180/pi;%radians to degrees
dtor=pi/180;%degrees to radians
con=true;
sep=100;%seperations between the shock and cone for runge kutta method
if theta_s>theta_c
    while con 
        %oplique shock relations
        Mn_fs=M_fs*sin(dtor*theta_s);
        Mn2=sqrt((1+0.5*(g-1)*Mn_fs^2)/(g*Mn_fs^2-0.5*(g-1)));
        del=atan(2*cot(dtor*theta_s)*((M_fs*(sin(dtor*theta_s))^2-1)/(M_fs^2*(g+(cos(dtor*theta_s))^2)+2)));;
        M2=Mn2/sin(dtor*(theta_s-(del*rtod)));
        alpha=theta_s-(del*rtod);
        %Normalised velocity behind shock
        V_=1/sqrt(2/((g-1)*M2^2)+1);
        Vr_=V_*cos(dtor*alpha);
        Vt_=V_*sin(dtor*alpha);
        %Runge Kutta Method
        func=@(t,y) (y^2*Vr_-0.5*(g-1)*(1-Vr_^2-y^2)*(2*Vr_+y*cot(dtor*t)))/(0.5*(g-1)*(1-Vr_^2-y^2)-y^2);
        theta(1)=theta_s;
        y(1)=Vt_;
        h=-(theta_s-theta_c)/sep;
        for i=1:sep
            K1=func(theta(i),y(i));
            K2=func(theta(i)+h/2,y(i)+K1*h/2);
            K3=func(theta(i)+h/2,y(i)+K2*h/2);
            K4=func(theta(i)+h,y(i)+K3*h);
            y(i+1)=y(i)+(1/6)*(K1+2*K2+2*K3+K4)*h;
            theta(i+1)=theta(1)+i*h;
        end
        if y(sep+1)<=0.00001 && y(sep+1)>=-0.00001
            con=false;
        else
            con=true;
            theta_sn=theta_s+(y(sep+1)*0.0000001/(-y(sep+1)));
            if theta_sn<theta_c
                theta_s=theta_s-0.1;
            else
                theta_s=theta_sn;
            end
        end
    end
end
disp(theta_s*dtor)
