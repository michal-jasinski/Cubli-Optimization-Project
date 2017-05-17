% Convert rotation matrix to quaternion
function [P_IK]=A_IK2Quat(A_IK)
Log=logm(A_IK);
v=[-Log(2,3);Log(1,3);-Log(1,2)];
phi=norm(v);
n=v/phi;
P_IK=[cos(phi/2);n*sin(phi/2)];
end