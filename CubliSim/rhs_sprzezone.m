function [ d_psi,gradient ] = rhs_sprzezone(t,psi,x,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    
    g1 = x(1); g2 = x(2); g3 = x(3);
    pwh1 = x(4); pwh2 = x(5); pwh3 = x(6);
    pww1 = x(7); pww2 = x(8); pww3 = x(9); 
    PIK1 = x(10); PIK2 = x(11); PIK3 = x(12); PIK4 = x(13);

    T1 = u(1); T2 = u(2); T3 = u(3);
    
    get_f_and_f_sprzezona;
    
    H = transpose(psi)*f - x(14);
    d1_f_sprzezona = double(d1_f_sprzezona);
    d_psi=-d1_f_sprzezona*psi;
    gradient=[ psi(7); psi(8); psi(9)];

end