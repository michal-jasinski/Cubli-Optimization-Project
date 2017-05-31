% Implementation of the Cubli-Corner Balancing according to the CDC-Article
% Author: Michael Muehlebach

% IMPORTANT: Note that the controller is implemented as in the CDC-Article. 
% It is therefore not the newest version!

clear all;
close all;

format long e; format compact;

%% initialization script
init;

u0=zeros(3*nodes,1);
LB = [-0.5; -0.5; -0.5];
UB = [0.5; 0.5; 0.5];

options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region',...
                        'Display','iter-detailed');
retval = fminunc(@get_quality_indicator,u0,options)

%  retval = [
%      6.445541347944299e-02     2.850171357943448e-02    -9.653822335520916e-02
%      1.625001297773137e-02    -3.354873469253270e-03    -1.402923768392790e-02
%     -1.255918671396476e-02    -7.460134345625878e-03     1.965638068591700e-02
%     -7.342160234901270e-03     2.091839125267733e-03     7.010486629167909e-03
%      3.132421821444559e-03     8.531401979116125e-04    -5.491867853553403e-03
%     -5.540889298720768e-04     1.212774450472074e-03    -3.940894016203839e-05
%      6.664233550356398e-04     3.090097643139585e-05    -7.345738477607324e-04
%      8.925533366707238e-04    -7.136477506431300e-03     5.717558671291215e-03
%      6.572060195940508e-05     4.061753763550837e-03    -3.748590736695321e-03
%     -7.973789563574593e-04     1.114378733291850e-03    -2.824600455530665e-04];
 
% simulation
[t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,retval,nodes);

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

figure(4)
plot(t,wh)
ylabel('\omega_h [rad/s]')
xlabel('time [s]')

% figure
% plot(t,ww/(2*pi)*60)
% ylabel('\omega_w [rpm]')
% xlabel('time [s]')

quality_indicator = x(:,14);
figure(5);
plot(t,quality_indicator);
