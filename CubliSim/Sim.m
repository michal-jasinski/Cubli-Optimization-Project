% Implementation of the Cubli-Corner Balancing according to the CDC-Article
% Author: Michael Muehlebach

% IMPORTANT: Note that the controller is implemented as in the CDC-Article. 
% It is therefore not the newest version!

clear all;
close all;

format long e; format compact;

%% initialization script
init

u0 = [4,4,4;1,1,1;3,2,3;0 0 0;3,3,3 ]
u0=zeros(5,3);
% for i = 1 : 10
%       u0(i,1) = 1.205787;
%       u0(i,2) = 0.4361158;
%       u0(i,3) = -1.720702;
% end

u = [ 3.759137302933496e-02     7.957694477885314e-03    -7.051632208034343e-02
     1.191776768252162e-02    -1.057997480541581e-02    -1.707723219940835e-02
     3.237050511267915e-03    -1.366680914984241e-02    -7.231339670213214e-04
    -2.189740416313436e-03    -1.000065197291657e-02     7.579525565884447e-03
     1.175117177712085e-03    -2.718723914399226e-03     1.257035932692729e-03];
global Q_hist;
Q_hist = zeros(10000,1);
global iter;
iter = 1;
options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region');
retval = fminunc(@bfgs,u0,options)


% simulation
[t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,retval);

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
plot(t,acos(I_m(:,1)/norm(m))/pi*180);
xlabel('time [s]')
ylabel('inclination angle [deg]')
figure(2)
hold all;
plot(t,acos(I_m(:,2)/norm(m))/pi*180);
xlabel('time [s]')
ylabel('inclination angle [deg]')
figure(3)
hold all;
plot(t,acos(I_m(:,3)/norm(m))/pi*180);
xlabel('time [s]')
ylabel('inclination angle [deg]')
% body angulr velocity
wh=(pwh-pww)*Theta_0_ht^-1;

% wheel velocity
ww=pww*Theta_w^-1-wh;

% figure
% plot(t,wh)
% ylabel('\omega_h [rad/s]')
% xlabel('time [s]')

% figure
% plot(t,ww/(2*pi)*60)
% ylabel('\omega_w [rpm]')
% xlabel('time [s]')

quality_indicator = x(:,14);
figure(4);
plot(quality_indicator);
