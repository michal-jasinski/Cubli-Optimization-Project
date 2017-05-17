clear all;
close all;
syms Q
syms g wh pwh pw m1 m2 m3
syms g1 g2 g3
syms wh1 wh2 wh3
syms pwh1 pwh2 pwh3
syms pww1 pww2 pww3
syms T1 T2 T3
syms PIK1 PIK2 PIK3 PIK4
syms t0h11 t0h12 t0h13 t0h21 t0h22 t0h23 t0h31 t0h32 t0h33
PIK=[PIK1; PIK2; PIK3; PIK4];
Theta_0_ht=[t0h11 t0h12 t0h13;t0h21 t0h22 t0h23; t0h31 t0h32 t0h33];
phi_ht=1/2*[-PIK(2:4), PIK(1)*eye(3)-Skew(PIK(2:4))];
T=[T1; T2; T3];
m=[m1;m2;m3];
g=[g1; g2; g3];
pwh=[pwh1; pwh2; pwh3];
wh=[wh1; wh2; wh3];
pww=[pww1; pww2; pww3];

dg=cross(-wh,g);
dpwh=cross(-wh,pwh)+cross(m,g);
dpww=T;
dPIK=transpose(phi_ht)*Theta_0_ht^-1*(pwh-pww);
A_IKtmp=eye(3)+2*PIK(1)*Skew(PIK(2:4))+2*Skew(PIK(2:4))^2;
I_m=transpose(A_IKtmp*m);
dQ=(I_m(1)/norm(m))^2+(I_m(2)/norm(m))^2+(I_m(3)/norm(m)-1)^2;

X=[g; pwh;pww; PIK; Q];
dX=[dg;dpwh;dpww;dPIK;dQ];
for i=1:14
    f(i)=dX(i);
end

for i =  1 : 14
    for j = 1 : 14
        dF(i,j) = diff(f(j),X(i));
    end
end


d1_f_sprzezona = subs(dF,[m,Theta_0_ht],[[0.1661; 0.1473; 0.1537],...
         [0.0304 -0.0130 -0.0135; -0.013 0.0342 -0.0128; -0.0135 -0.0128 0.0331]]);
     
matrix2 = subs(d1_f_sprzezona,[wh,pwh,pww,g],[[1;1;1],[1;1;1],[1;1;1],[1;1;1]]);
matrix2 = double(subs(matrix2,[PIK],[1;1;1;1]))

