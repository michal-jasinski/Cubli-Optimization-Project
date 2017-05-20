function [t,x,psi] = rk4(rhs,rhs_sprzezone,x0,tf,sample_time,Theta_0_ht,m)
%UNTITLED Summary of this function goes here

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

% Solution of state equations
for j = 1 : length(dtau)
    h = dtau(j)/n(j);
    h_2=h/2; h_6=h/6; h_3=h/3;
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
        x(i+1,:) = (xtmp+h_6*(dx1+dx4)+h_3*(dx2+dx3))';
        t(i+1) = t(i) + h;
    end
end

psi = zeros(size(x));
psi(end,:)=[zeros(13,1);-1]; % do przemyœlenia

% Solution of equilibrium equations
for j = length(dtau): -1 : 1
    h = dtau(j)/n(j);
    h_2=h/2; h_6=h/6; h_3=h/3; h_8=h/8;
    % variables used in calculations
    dpsi1=zeros(length(x0),1);dpsi2=dpsi1;dpsi3=dpsi1;dpsi4=dpsi1;
    dxp = rhs(t,x(cn(j+1),:)',Theta_0_ht,m,u(j,:)');
    for i = cn(j+1) : -1 : cn(j) + 1
        psi_tmp = psi(i,:)';
        dpsi1   = rhs_sprzezone(t,psi_tmp,x(i,:)',u(j,:)');
        tmp     = psi_tmp-h_2*dpsi1;
        dxl     = rhs(t,x(i-1,:)',Theta_0_ht,m,u(j,:)');
        xp      = (x(i-1,:)'+x(i,:)')/2 + (dxl-dxp)*h_8;
        dxp     = dxl;
        dpsi2   = rhs_sprzezone(t,tmp,xp,u(j,:)');
        tmp     = psi_tmp-h_2*dpsi2;
        dpsi3   = rhs_sprzezone(t,tmp,xp,u(j,:)');
        tmp     = psi_tmp-h*dpsi3;
        dpsi4   = rhs_sprzezone(t,tmp,x(i-1,:)',u(j,:)'); 
        psi(i-1,:) = psi_tmp - h_3*(dpsi2+dpsi3) - h_6*(dpsi1+dpsi4);
    end
end

end

