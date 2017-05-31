function [ Q, gradient] = bfgs(u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global Q_hist;
    global iter;
    init;
    disp(u)
    [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u);
    
%     g=x(1:3);
% p_wh=x(4:6);
% p_ww=x(7:9);
% PIK=x(10:13);
% Theta_0_ht=[0.0304   -0.0130   -0.0135
%    -0.0130    0.0342   -0.0128
%    -0.0135   -0.0128    0.0331];

% angular velocity
% wh=Theta_0_ht^-1*(p_wh(:,end)-p_ww(:,end));
    Q = x(end,end)%+wh(1)^2+wh(2)^2+wh(3)^2;
    Q_hist(iter) = Q;
    iter = iter + 1
    Q
    
end

