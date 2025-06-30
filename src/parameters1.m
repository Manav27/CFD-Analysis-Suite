%Engine
m_dot=132;
burn_time=162.25;
%vehicle
m_dry=1360.7;
m_payload=5000;
m_nofuel=m_dry+m_payload;
%Initial conditions
v_0=0.01;
g=9.81;
gam_0=pi/2;
R_e=6371e3;
h_0=0;
x_0=0;
%pulse gamma input for gravity turn initiation
t_turn=30;
gam_in=0.1;