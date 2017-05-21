% Implementation of the Cubli-Corner Balancing according to the CDC-Article
% Author: Michael Muehlebach

% IMPORTANT: Note that the controller is implemented as in the CDC-Article. 
% It is therefore not the newest version!

clear all;
close all;

format long e; format compact;

%% initialization script
init

u0 = zeros(10,3);
for i = 1 : 10
      u0(i,1) = 1.205787;
      u0(i,2) = 0.4361158;
      u0(i,3) = -1.720702;
end

u = [4.373708296524337e+00     5.922144628611976e+00     4.067068252559787e+00;
     3.500507727691882e-01    -6.602550087301692e-01    -1.105923798422186e+00;
    -1.082807051052882e-01    -1.129203562429491e+00    -8.340948931193852e-01;
    -4.238753713801007e-02    -9.672587460129165e-01    -3.021832763793145e-01;
    -5.561735751136453e-01    -8.546500134290863e-02    -1.419593419378375e-01;
    -6.652087597678327e-01    -5.818740342414767e-01     6.263123647229792e-01;
    -7.401912392012053e-01    -3.227300336945300e-01     1.758238393939733e-01;
    -1.136063254582653e-01     2.509284589290787e-01     8.319370514404725e-03;
     5.915408210206344e-01     9.411657490632475e-02     1.379199886681730e-01;
     2.576740803369788e+00     2.704385809540259e+00     2.517690986867580e+00];
% global Q_hist;
% Q_hist = zeros(10000,1);
% global iter;
% iter = 1;
% options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region');
% retval = fminunc(@bfgs,u0,options)


% simulation
[t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u);

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
wh=(pwh-pww)*Theta_0_ht^-1;

% wheel velocity
ww=pww*Theta_w^-1-wh;

figure
plot(t,wh)
ylabel('\omega_h [rad/s]')
xlabel('time [s]')

figure
plot(t,ww/(2*pi)*60)
ylabel('\omega_w [rpm]')
xlabel('time [s]')

quality_indicator = x(:,14);
figure(2);
plot(quality_indicator);
