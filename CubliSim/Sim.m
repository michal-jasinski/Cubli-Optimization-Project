% Implementation of the Cubli-Corner Balancing according to the CDC-Article
% Main Author: Michael Muehlebach
% Script has been modified by Jakub Wasik i Michal Jasinski
% in order to find optimal control by solving optimization problem

clear all;
close all;

format long e; format compact;

%% initialization script
init;

u0 = zeros(3*nodes,1);
lb = -2;
ub = 2;
epsilon = 2e-4;
maxIter = 100;

% options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region',...
%                         'Display','iter-detailed');
% optimal_control = fminunc(@get_quality_indicator,u0,options)

% optimal_control = bfgs(u0,epsilon,lb,ub,maxIter);
 
 optimal_control = [
     8.038229602639262e-01
     2.081528775122536e-01
    -1.066730030711171e+00
    -3.324912476177477e-01
    -9.212437975265829e-02
     4.455989354232792e-01
     3.608807086097594e-02
     6.736143862505707e-04
    -3.905948199267195e-02
    -9.858984023725736e-03
     3.337695128563527e-03
     8.386190754135091e-03
    -3.127088776896491e-02
    -4.475588568594685e-03
     3.703750421403575e-02];

% simulation
[t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,...
                         Theta_0_ht,m,optimal_control,nodes);

% Extract data
PIK = x(:,10:13);
pwh = x(:,4:6);
pww = x(:,7:9);
quality_indicator = x(:,14);

% express m-vector in inertial coordinate frame
I_m = zeros(length(t),3);
I_pwh = zeros(length(t),3);

for k = 1 : length(t)
    A_IKtmp    = eye(3)+2*PIK(k,1)*Skew(PIK(k,2:4))+2*Skew(PIK(k,2:4))^2;
    I_m(k,:)   = (A_IKtmp*m)';
    I_pwh(k,:) = (A_IKtmp*(pwh(k,:)'))';
end

% figure(1); hold on; grid on;
% plot(t,acos(I_m(:,1)/norm(m))/pi*180);
% plot(t,acos(I_m(:,2)/norm(m))/pi*180);
% plot(t,acos(I_m(:,3)/norm(m))/pi*180);
% xlabel('Czas [s]')
% ylabel('K�t nachylenia kostki [deg]');
% legend('X', 'Y', 'Z');


% body angular velocity
wh = (pwh-pww)*Theta_0_ht^-1;

% wheel velocity
ww = pww*Theta_w^-1-wh;

% figure(2); hold on; grid on;
% plot(t,wh);
% ylabel('Pr�dko�� k�towa kostki [rad/s]');
% xlabel('Czas [s]');
% legend('X','Y','Z');
% 
% figure(3); hold on; grid on;
% plot(t,ww/(2*pi)*60);
% ylabel('Pr�dko�� obrotowa k� reakcyjnych [rpm]');
% xlabel('Czas [s]');
% legend('X','Y','Z');
% 
% figure(4);
% plot(t,quality_indicator);
% xlabel('Czas [s]');
% ylabel('Warto�� wska�nika jako�ci');
% grid on;

figure(5); hold all; 

subplot(2,2,1); hold on; grid on;
plot(t,quality_indicator);
xlabel('Czas [s]');
ylabel('Warto�� wska�nika jako�ci');

subplot(2,2,2); hold on; grid on;
plot(t,acos(I_m(:,1)/norm(m))/pi*180);
plot(t,acos(I_m(:,2)/norm(m))/pi*180);
plot(t,acos(I_m(:,3)/norm(m))/pi*180);
xlabel('Czas [s]')
ylabel('K�t nachylenia kostki [deg]');
legend('X', 'Y', 'Z');

subplot(2,2,3); hold on; grid on;
plot(t,wh);
ylabel('Pr�dko�� k�towa kostki [rad/s]');
xlabel('Czas [s]');
legend('X','Y','Z');

subplot(2,2,4); hold on; grid on;
plot(t,ww/(2*pi)*60);
ylabel('Pr�dko�� obrotowa k� reakcyjnych [rpm]');
xlabel('Czas [s]');
legend('X','Y','Z');

hold off;
