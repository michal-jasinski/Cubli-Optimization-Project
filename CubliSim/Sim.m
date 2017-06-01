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
lb = -2;
ub = 2;
epsilon = 2e-4;
maxIter = 500;

% options = optimoptions('fminunc','GradObj','on','Algorithm','trust-region',...
%                         'Display','iter-detailed');
% retval = fminunc(@get_quality_indicator,u0,options)

retval = bfgs(u0,epsilon,lb,ub,maxIter)

 
% retval = [
%      6.934161042971873e-02
%      2.563401714063781e-02
%     -9.941398766136766e-02
%      2.234805925923642e-03
%      1.162939769305259e-03
%     -3.544536592805359e-03
%     -6.408200313467240e-03
%     -2.581676789982574e-03
%      9.350003109377073e-03
%     -3.970712885263781e-03
%     -1.793001106040176e-03
%      5.987394154127943e-03
%     -1.196767260802018e-03
%     -5.829287977509501e-04
%      1.840861075819404e-03];
%  
% % simulation
% [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,retval,nodes);
% 
% % Extract data
% PIK=x(:,10:13);
% pwh=x(:,4:6);
% pww=x(:,7:9);
% 
% % express m-vector in inertial coordinate frame
% I_m=zeros(length(t),3);
% I_pwh=zeros(length(t),3);
% 
% for k=1:length(t)
%     A_IKtmp=eye(3)+2*PIK(k,1)*Skew(PIK(k,2:4))+2*Skew(PIK(k,2:4))^2;
%     I_m(k,:)=(A_IKtmp*m)';
%     I_pwh(k,:)=(A_IKtmp*(pwh(k,:)'))';
% end
% 
% % calculate inclination angle
% figure(1)
% hold all;
% plot(t,acos(I_m(:,1)/norm(m))/pi*180);
% xlabel('time [s]')
% ylabel('inclination angle [deg]')
% figure(2)
% hold all;
% plot(t,acos(I_m(:,2)/norm(m))/pi*180);
% xlabel('time [s]')
% ylabel('inclination angle [deg]')
% figure(3)
% hold all;
% plot(t,acos(I_m(:,3)/norm(m))/pi*180);
% xlabel('time [s]')
% ylabel('inclination angle [deg]')
% % body angulr velocity
% wh=(pwh-pww)*Theta_0_ht^-1;
% 
% % wheel velocity
% ww=pww*Theta_w^-1-wh;
% 
% figure(4)
% plot(t,wh)
% ylabel('\omega_h [rad/s]')
% xlabel('time [s]')
% 
% % figure
% % plot(t,ww/(2*pi)*60)
% % ylabel('\omega_w [rpm]')
% % xlabel('time [s]')
% 
% quality_indicator = x(:,14);
% figure(5);
% plot(t,quality_indicator);
