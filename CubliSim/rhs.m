function dx = rhs(t,x,Theta_0_ht,m,alpha,beta,gamma,delta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
g=x(1:3);
p_wh=x(4:6);
p_ww=x(7:9);
PIK=x(10:13);
Q=x(14);
%[PIK_dot, d_x15] = kinematics(p_wh,p_ww,PIK,Theta_0_ht,m)
%[dot_g, dot_pwh, dot_pww]= fcn(T,g,p_wh,p_ww,Theta_0_ht,m)

% angular velocity
wh=Theta_0_ht^-1*(p_wh-p_ww);

%control
p_wh_perp=p_wh-g*g'*p_wh/(9.81^2);

% constants for nonlinear controller
K1=(1+beta*gamma+delta)*eye(3)+alpha*Theta_0_ht;
K2=alpha*Theta_0_ht*Skew(p_wh_perp)+beta*Skew(m)*Skew(g)+Skew(p_wh);
K3=gamma*(eye(3)+alpha*Theta_0_ht*(eye(3)-g*g'/(9.81^2)));
K4=gamma*eye(3);
global control;
% control law
u=K1*Skew(m)*g+K2*wh+K3*p_wh-K4*p_ww;
% u=[1;1;1];
control=[control; u'];
% Cubli dynamics:
dot_g=-cross(wh,g);
dot_pwh=-cross(wh,p_wh)+cross(m,g);
dot_pww=u;

phi_ht=1/2*[-PIK(2:4), PIK(1)*eye(3)-Skew(PIK(2:4))];

% kinematics
PIK_dot=phi_ht'*Theta_0_ht^(-1)*(p_wh-p_ww);

% quality indicator
A_IKtmp=eye(3)+2*PIK(1)*Skew(PIK(2:4))+2*Skew(PIK(2:4))^2;
I_m=(A_IKtmp*m)';

d_Q = (I_m(1)/norm(m))^2 + (I_m(2)/norm(m))^2 + (I_m(3)/norm(m)-1)^2;
dx=[dot_g;dot_pwh;dot_pww;PIK_dot;d_Q];
end

