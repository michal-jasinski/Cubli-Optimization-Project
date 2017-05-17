clear all;
close all;

phi=1/180*pi;   % initial inclination, rotate around x-axis of the inertial frame
Theta_0_ht=[0.0304   -0.0130   -0.0135
   -0.0130    0.0342   -0.0128
   -0.0135   -0.0128    0.0331];
% Cubli initially at rest

PIK_0=[0.8774
    0.3525
   -0.3018
    0.1213];
PIK_0=PIK_0/norm(PIK_0);
Q_0 = 0;
R_IK_0=eye(3)+2*PIK_0(1)*Skew(PIK_0(2:4))+2*Skew(PIK_0(2:4))^2;
g_0=R_IK_0'*[0;0;-9.81];
% Calculate the corresponding rotation matrix
pww_0=[0;0;0];
pwh_0=[0;0;0];
sample_time = 0.005;
time = 2;
tspan = 0:sample_time:time;

epsilon=0.03;
x0 = [g_0;pwh_0;pww_0;PIK_0;Q_0];
for i=1: length(x0)
    x0_epsilon(i)=x0(i)+epsilon
    [t1,x_epsilon,psi_epsilon] = rk4(@rhs,x0_epsilon,time,sample_time,Theta_0_ht,m);
    [t,x,psi] = rk4(@rhs,x0,time,sample_time,Theta_0_ht,m);
    dQ_dx0(i)=(x_epsilon(end,end)-x(end,end))/epsilon;
    x0_epsilon=x0;
end
x0_epsilon
