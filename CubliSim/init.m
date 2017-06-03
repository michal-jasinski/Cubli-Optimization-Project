%% Define parameters

% Inertia
Theta_0_ht=[0.0304   -0.0130   -0.0135
   -0.0130    0.0342   -0.0128
   -0.0135   -0.0128    0.0331];

% Inertia of the reaction wheels
Theta_w=5*1e-4*eye(3);

% Quaternion describing the equilibrium orientation
PIK_0=[0.8774
    0.3525
   -0.3018
    0.1213];
PIK_0=PIK_0/norm(PIK_0);

% Calculate the corresponding rotation matrix
R_IK_0=eye(3)+2*PIK_0(1)*Skew(PIK_0(2:4))+2*Skew(PIK_0(2:4))^2;

% Calculate the gravity vector in body coordinate frame (equilibrium
% position)
g_0=R_IK_0'*[0;0;-9.81];

% Calculate the m vector (vector to centre of mass * mass)
m=-g_0/9.81*0.27;

%% Initial conditions

% define initial condtions
phi=1/20*pi;   % initial inclination, rotate around x-axis of the inertial frame

A_IKtmp=expm(Skew([1;1;1])*phi)*R_IK_0;
PIK_0=A_IK2Quat(A_IKtmp);   % form the quaternion from the rotation matrix
g_0=A_IKtmp'*[0;0;-9.81];   % form the gravity vector in the body coordinate frame

% Cubli initially at rest
pww_0=[0;0;0];
pwh_0=[0;0;0];
Q_0 = 0;

x0 = [g_0;pwh_0;pww_0;PIK_0;Q_0];
sample_time = 0.005;
time = 0.5;
nodes = 5;