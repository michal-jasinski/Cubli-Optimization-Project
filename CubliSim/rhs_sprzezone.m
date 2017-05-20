function [ d_psi ] = rhs_sprzezone(t,psi,x,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    global d1_f_sprzezona;
    syms Q
    syms g pwh pw m1 m2 m3
    syms g1 g2 g3
    syms wh1 wh2 wh3
    syms pwh1 pwh2 pwh3
    syms pww1 pww2 pww3
    syms PIK1 PIK2 PIK3 PIK4
    
    matrix = double(subs(d1_f_sprzezona,[g1,g2,g3,pwh1,pwh2,pwh3,pww1,pww2,pww3,PIK1,PIK2,PIK3,PIK4],x(1:13)'));
    d_psi=-matrix*psi;
    
end