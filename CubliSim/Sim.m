% Implementation of the Cubli-Corner Balancing according to the CDC-Article
% Author: Michael Muehlebach

% IMPORTANT: Note that the controller is implemented as in the CDC-Article. 
% It is therefore not the newest version!

clear all;
close all;

format long e; format compact;

%% initialization script
init;

u0=zeros(nodes,3);

options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region',...
                        'Display','iter-detailed');
retval = fminunc(@bfgs,u0,options)

% retval = [
%      7.735571067387728e-02     2.967611502000844e-02    -1.110622332274811e-01
%     -1.498368789357186e-02    -7.405833694822156e-03     2.106634914182178e-02
%     -9.345755394289819e-03    -3.225632194952382e-03     1.522620423314977e-02
%      2.649579858168301e-03     7.044889982018520e-04    -4.452646728321130e-03
%      3.376866417163700e-04     5.028336823312962e-05    -3.518435708294290e-04];
 
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
