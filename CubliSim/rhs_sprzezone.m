function [ d_psi ] = rhs_sprzezone(t,psi,x,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global d1_f_sprzezona;
    matrix = subs(d1_f_sprzezona,[wh,pwh,pww,g],x);
    d_psi=-matrix*psi;
end