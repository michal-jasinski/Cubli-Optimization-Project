function dx = rhs(t,x,Theta_0_ht,m,u)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
g=x(1:3);
p_wh=x(4:6);
p_ww=x(7:9);
PIK=x(10:13);
Q=x(14);

% angular velocity
wh=Theta_0_ht^-1*(p_wh-p_ww);

% control


% Cubli dynamis:
dot_g=-cross(wh,g);
dot_pwh=-cross(wh,p_wh)+cross(m,g);
dot_pww=u;

phi_ht=1/2*[-PIK(2:4), PIK(1)*eye(3)-Skew(PIK(2:4))];

% kinematics
PIK_dot=phi_ht'*Theta_0_ht^(-1)*(p_wh-p_ww);

% quality indicator
A_IKtmp=eye(3)+2*PIK(1)*Skew(PIK(2:4))+2*Skew(PIK(2:4))^2;
I_m=(A_IKtmp*m)';

d_Q = (I_m(1)/norm(m))^2 + (I_m(2)/norm(m))^2 + (I_m(3)/norm(m)-1)^2;%+(wh(1)^2+wh(2)^2+wh(3)^2)/1000;
dx=[dot_g;dot_pwh;dot_pww;PIK_dot;d_Q];
end

