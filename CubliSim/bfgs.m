function [ Q, gradient] = bfgs(u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    global Q_hist;
    global iter;
    init;
    disp(u)
    [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u);
    
    Q = x(end,end)
    Q_hist(iter) = Q;
    iter = iter + 1
    
end

