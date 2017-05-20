% Implementation of the Cubli-Corner Balancing according to the CDC-Article
% Author: Michael Muehlebach

% IMPORTANT: Note that the controller is implemented as in the CDC-Article. 
% It is therefore not the newest version!

clear all;
close all;
format long e; format compact;

%% Define parameters

% Inertia
Theta_0_ht=[0.0304   -0.0130   -0.0135
   -0.0130    0.0342   -0.0128
   -0.0135   -0.0128    0.0331];

% Inertia of the reaction wheels
Theta_w=5*1e-4*eye(3);

% Quaternion describing the equilibrium orientation
PIK_0=[0.8774
    0.3525
   -0.3018
    0.1213];
PIK_0=PIK_0/norm(PIK_0);

% Calculate the corresponding rotation matrix
R_IK_0=eye(3)+2*PIK_0(1)*Skew(PIK_0(2:4))+2*Skew(PIK_0(2:4))^2;

% Calculate the gravity vector in body coordinate frame (equilibrium
% position)
g_0=R_IK_0'*[0;0;-9.81];

% Calculate the m vector (vector to centre of mass * mass)
m=-g_0/9.81*0.27;

% Parameters of the nonlinear controller
alpha=40;
beta=3;
gamma=10;
delta=.001;


%%

% define initial condtions
phi=1/180*pi;   % initial inclination, rotate around x-axis of the inertial frame

A_IKtmp=expm(Skew([1;1;1])*phi)*R_IK_0;
PIK_0=A_IK2Quat(A_IKtmp);   % form the quaternion from the rotation matrix
g_0=A_IKtmp'*[0;0;-9.81];   % form the gravity vector in the body coordinate frame

% Cubli initially at rest
pww_0=[0;0;0];
pwh_0=[0;0;0];

Q_0 = 0;

x0 = [g_0;pwh_0;pww_0;PIK_0;Q_0];
sample_time = 0.005;
time = 2;
tspan = 0:sample_time:time;

load('d1_f_sprzezona.mat');
global d1_f_sprzezona;
d1_f_sprzezona = rownania_sprzezone;

disp('symulacja');
% [t,x] = ode45(@(t,x)rhs(t,x,Theta_0_ht,m,alpha,beta,gamma,delta),tspan,x0);
[t,x,psi] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m);

% Extract data
PIK=x(:,10:13);
pwh=x(:,4:6);
pww=x(:,7:9);


% express m-vector in inertial coordinate frame
I_m=zeros(length(t),3);
I_pwh=zeros(length(t),3);

for k=1:length(t)
    A_IKtmp=eye(3)+2*PIK(k,1)*Skew(PIK(k,2:4))+2*Skew(PIK(k,2:4))^2;
    I_m(k,:)=(A_IKtmp*m)';
    I_pwh(k,:)=(A_IKtmp*(pwh(k,:)'))';
end

% calculate inclination angle
figure(1)
hold all;
plot(t,acos(I_m(:,3)/norm(m))/pi*180);
xlabel('time [s]')
ylabel('inclination angle [deg]')

% body angulr velocity
% wh=(pwh-pww)*Theta_0_ht^-1;

% wheel velocity
% ww=pww*Theta_w^-1-wh;

% figure
% plot(t,wh)
% ylabel('\omega_h [rad/s]')
% xlabel('time [s]')
% 
% figure
% plot(t,ww/(2*pi)*60)
% ylabel('\omega_w [rpm]')
% xlabel('time [s]')

quality_indicator = x(:,14);
figure(2);
plot(quality_indicator);
