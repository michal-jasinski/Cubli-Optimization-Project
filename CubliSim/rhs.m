function dx = rhs(t,x,Theta_0_ht,m,u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
g=x(1:3);
p_wh=x(4:6);
p_ww=x(7:9);
PIK=x(10:13);
Q=x(14);

% angular velocity
wh=Theta_0_ht\(p_wh-p_ww);

% control

% Cubli dynamis:
wh_skew = [0, -wh(3), wh(2); wh(3), 0, -wh(1); -wh(2), wh(1), 0];

m_skew  = [0, -m(3), m(2); m(3), 0, -m(1); -m(2), m(1), 0];

PIK_tmp = PIK(2:4);
PIK_skew = [0, -PIK_tmp(3), PIK_tmp(2); PIK_tmp(3), 0, -PIK_tmp(1); -PIK_tmp(2), PIK_tmp(1), 0];

dot_g=-wh_skew*g;
dot_pwh=-wh_skew*p_wh + m_skew*g;
dot_pww=u;

phi_ht=1/2*[-PIK(2:4), PIK(1)*eye(3)-PIK_skew];

% kinematics
PIK_dot=phi_ht'*(Theta_0_ht\(p_wh-p_ww)); % --> slower;
% PIK_dot = phi_ht'*wh;

% quality indicator
A_IKtmp=eye(3)+2*PIK(1)*PIK_skew+2*PIK_skew^2;
I_m=(A_IKtmp*m)';

norm_m = norm(m);
d_Q = (I_m(1)/norm_m)^2 + (I_m(2)/norm_m)^2 + (I_m(3)/norm_m-1)^2 + 0.01*(wh(1)^2+wh(2)^2+wh(3)^2);
dx=[dot_g;dot_pwh;dot_pww;PIK_dot;d_Q];
end

