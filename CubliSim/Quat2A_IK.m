function A_IK=Quat2A_IK(P_IK)
% Get rotation matrix A_IK from Quaternion P_IK
A_IK=eye(3)+2*P_IK(1)*Skew(P_IK(2:4))+2*Skew(P_IK(2:4))^2;
end
