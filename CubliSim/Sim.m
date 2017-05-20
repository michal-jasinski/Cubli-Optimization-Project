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

u = [6.076814227052029e+00     5.615751902410058e+00     3.857374217866485e+00;
     1.530025558502412e+00     4.019225264973242e-01    -1.442684045935703e+00;
     7.469138642942980e-01     3.669637219037544e-02    -2.127308544612713e+00;
     9.111888303822799e-01     3.415254653678650e-02    -1.910714790871892e+00;
     9.713902278496580e-01     3.431464402297310e-01    -1.952542608602176e+00;
     1.271400402448310e+00     3.931807805072060e-01    -1.542432411286200e+00;
     1.407036749905861e+00     7.972608715652225e-01    -1.355409838356558e+00;
     1.645781707124391e+00     7.545150979420012e-01    -1.368753829776202e+00;
     1.841261850690880e+00     1.113846290879535e+00    -9.971467945291672e-01;
     2.043183623749663e+00     1.329691617576480e+00    -8.698385650049383e-01];

% options = optimoptions('fminunc','GradObj','on');
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
