function [ Q, gradient] = bfgs(u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    init;
    disp(u)
    [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u);

    Q = x(end,end)

end

