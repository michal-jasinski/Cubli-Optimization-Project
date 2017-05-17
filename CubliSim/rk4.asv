function [t,x] = rk4(rhs,x0,tf,sample_time,Theta_0_ht,m)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   @params:
%   rhs : right side of the differential equation
%   x0  : initial condition vector
%   tf  : duration of the simulation
%   I, m, g, T : params used in rhs

load('optimal_control');

% czas
T = tf;

% okres próbkowania
h0 = sample_time;

% wêz³y strukturalne
tau = linspace(0,T,11);

% odleg³oœæ czasowa pomiêdzy tymi wêz³ami
dtau = diff(tau);
s = size(dtau);

% sterowanie
u = zeros(s(2),3);
for i = 1 : s(2)
%     u(i,:) = control(2,:);
      u(i,:) = 1;
end
% iloœæ próbek na ka¿dy wêze³ strukturalny
n = ceil(dtau/h0);

% wêz³y dyskretyzacji cn(i) numer wêz³a dyskretyzacji
% bêd¹cego i-tym wêz³em strukturalnym
cn = cumsum([1 n]);

x = zeros(cn(end),length(x0));   % matrix for the data (nt x n)
t = zeros(cn(end),1);

tmp  = zeros(length(x0),1);  % (3 x 1) vector of zeros for calculations
x(1,:) = x0';       % first element in the data matrix -> x0


% Solution of equation of state
for j = 1 : length(dtau)
    h = dtau(j)/n(j);
    h_2=h/2; h_6=h/6; h_26=2*h_6;
    % variables used in calculations
    dx1=zeros(length(x0),1);dx2=dx1;dx3=dx1;dx4=dx1;
    for i = cn(j) : cn(j+1) - 1
        xtmp = x(i,:)';
        dx1  = rhs(t,xtmp,Theta_0_ht,m,u(j,:)');
        tmp  = xtmp+h_2*dx1;
        dx2  = rhs(t,tmp,Theta_0_ht,m,u(j,:)');
        tmp  = xtmp+h_2*dx2;
        dx3  = rhs(t,tmp,Theta_0_ht,m,u(j,:)');
        tmp  = xtmp+h*dx3;
        dx4  = rhs(t,tmp,Theta_0_ht,m,u(j,:)');
        x(i+1,:) = (xtmp+h_6*(dx1+dx4)+h_26*(dx2+dx3))';
        t(i+1) = t(i) + h;
    end
end

Q = x(end,end);
psi = zeros(size(x));
psi(end,:)=[zeros(13,1);-1];

% Solution of equation of state
for j = length(dtau): -1 : 1
    h = dtau(j)/n(j);
    h_2=h/2; h_6=h/6; h_26=2*h_6;
    % variables used in calculations
    dpsi1=zeros(length(x0),1);dpsi2=dpsi1;dpsi3=dpsi1;dpsi4=dpsi1;
    for i = cn(j+1) : cn(j) + 1
        psi_tmp = psi(i,:)';
        dpsi1  = rhs(t,psi_tmp,Theta_0_ht,m,u(j,:)');
        tmp  = psi_tmp+h_2*dpsi1;
        dpsi2  = rhs(t,tmp,Theta_0_ht,m,u(j,:)');
        tmp  = psi_tmp+h_2*dpsi2;
        dpsi3  = rhs(t,tmp,Theta_0_ht,m,u(j,:)');
        tmp  = psitmp+h*dpsi3;
        dpsi4  = rhs(t,tmp,Theta_0_ht,m,u(j,:)');
        x(i+1,:) = (psitmp+h_6*(dpsi1+dpsi4)+h_26*(dpsi2+dpsi3))';
        t(i+1) = t(i) + h;
    end
end

end

