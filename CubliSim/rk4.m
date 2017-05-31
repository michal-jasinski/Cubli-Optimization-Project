function [t,x,psi,gradient] = rk4(rhs,rhs_sprzezone,x0,tf,sample_time,Theta_0_ht,m,u)
%UNTITLED Summary of this function goes here

% czas
T = tf;

% okres próbkowania
h0 = sample_time;

% wêz³y strukturalne
tau = linspace(0,T,6);

% odleg³oœæ czasowa pomiêdzy tymi wêz³ami
dtau = diff(tau);
% s = size(dtau);

% sterowanie
% u = zeros(s(2),3);
% for i = 1 : s(2)
%       u(i,:) = 1;
% end

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
g1 = x(end,1); g2 = x(end,2); g3 = x(end,3);
    pwh1 = x(end,4); pwh2 = x(end,5); pwh3 = x(end,6);
    pww1 = x(end,7); pww2 = x(end,8); pww3 = x(end,9); 
    PIK1 = x(end,10); PIK2 = x(end,11); PIK3 = x(end,12); PIK4 = x(end,13);
psi = zeros(size(x));
%psi(end,:)=[ 0, 0, 0, 0, 0, 0, 0, 0, 0, (200000000*((1537*PIK2)/5000 - (1661*PIK4)/5000)*((1661*PIK1*PIK4)/5000 - (1537*PIK1*PIK2)/5000 + (1661*PIK2*PIK3)/5000 + (1537*PIK3*PIK4)/5000 - (1473*PIK2^2)/5000 - (1473*PIK4^2)/5000 + 1473/10000))/7291019 - (200000000*((1537*PIK3)/5000 - (1473*PIK4)/5000)*((1537*PIK1*PIK3)/5000 - (1473*PIK1*PIK4)/5000 + (1473*PIK2*PIK3)/5000 + (1537*PIK2*PIK4)/5000 - (1661*PIK3^2)/5000 - (1661*PIK4^2)/5000 + 1661/10000))/7291019 - (20000*7291019^(1/2)*((10000*7291019^(1/2)*((1473*PIK1*PIK2)/5000 - (1661*PIK1*PIK3)/5000 + (1661*PIK2*PIK4)/5000 + (1473*PIK3*PIK4)/5000 - (1537*PIK2^2)/5000 - (1537*PIK3^2)/5000 + 1537/10000))/7291019 - 1)*((1473*PIK2)/5000 - (1661*PIK3)/5000))/7291019, (200000000*((1537*PIK1)/5000 + (1473*PIK2)/2500 - (1661*PIK3)/5000)*((1661*PIK1*PIK4)/5000 - (1537*PIK1*PIK2)/5000 + (1661*PIK2*PIK3)/5000 + (1537*PIK3*PIK4)/5000 - (1473*PIK2^2)/5000 - (1473*PIK4^2)/5000 + 1473/10000))/7291019 - (200000000*((1473*PIK3)/5000 + (1537*PIK4)/5000)*((1537*PIK1*PIK3)/5000 - (1473*PIK1*PIK4)/5000 + (1473*PIK2*PIK3)/5000 + (1537*PIK2*PIK4)/5000 - (1661*PIK3^2)/5000 - (1661*PIK4^2)/5000 + 1661/10000))/7291019 - (20000*7291019^(1/2)*((10000*7291019^(1/2)*((1473*PIK1*PIK2)/5000 - (1661*PIK1*PIK3)/5000 + (1661*PIK2*PIK4)/5000 + (1473*PIK3*PIK4)/5000 - (1537*PIK2^2)/5000 - (1537*PIK3^2)/5000 + 1537/10000))/7291019 - 1)*((1473*PIK1)/5000 - (1537*PIK2)/2500 + (1661*PIK4)/5000))/7291019, (20000*7291019^(1/2)*((10000*7291019^(1/2)*((1473*PIK1*PIK2)/5000 - (1661*PIK1*PIK3)/5000 + (1661*PIK2*PIK4)/5000 + (1473*PIK3*PIK4)/5000 - (1537*PIK2^2)/5000 - (1537*PIK3^2)/5000 + 1537/10000))/7291019 - 1)*((1661*PIK1)/5000 + (1537*PIK3)/2500 - (1473*PIK4)/5000))/7291019 - (200000000*((1661*PIK2)/5000 + (1537*PIK4)/5000)*((1661*PIK1*PIK4)/5000 - (1537*PIK1*PIK2)/5000 + (1661*PIK2*PIK3)/5000 + (1537*PIK3*PIK4)/5000 - (1473*PIK2^2)/5000 - (1473*PIK4^2)/5000 + 1473/10000))/7291019 - (200000000*((1537*PIK1)/5000 + (1473*PIK2)/5000 - (1661*PIK3)/2500)*((1537*PIK1*PIK3)/5000 - (1473*PIK1*PIK4)/5000 + (1473*PIK2*PIK3)/5000 + (1537*PIK2*PIK4)/5000 - (1661*PIK3^2)/5000 - (1661*PIK4^2)/5000 + 1661/10000))/7291019, (200000000*((1473*PIK1)/5000 - (1537*PIK2)/5000 + (1661*PIK4)/2500)*((1537*PIK1*PIK3)/5000 - (1473*PIK1*PIK4)/5000 + (1473*PIK2*PIK3)/5000 + (1537*PIK2*PIK4)/5000 - (1661*PIK3^2)/5000 - (1661*PIK4^2)/5000 + 1661/10000))/7291019 - (200000000*((1661*PIK1)/5000 + (1537*PIK3)/5000 - (1473*PIK4)/2500)*((1661*PIK1*PIK4)/5000 - (1537*PIK1*PIK2)/5000 + (1661*PIK2*PIK3)/5000 + (1537*PIK3*PIK4)/5000 - (1473*PIK2^2)/5000 - (1473*PIK4^2)/5000 + 1473/10000))/7291019 - (20000*7291019^(1/2)*((10000*7291019^(1/2)*((1473*PIK1*PIK2)/5000 - (1661*PIK1*PIK3)/5000 + (1661*PIK2*PIK4)/5000 + (1473*PIK3*PIK4)/5000 - (1537*PIK2^2)/5000 - (1537*PIK3^2)/5000 + 1537/10000))/7291019 - 1)*((1661*PIK2)/5000 + (1473*PIK3)/5000))/7291019, -1];


psi(end,:)=[zeros(13,1);-1]; % do przemyœlenia
gradient=zeros(length(dtau),3);
% Solution of equilibrium equations

for j = length(dtau): -1 : 1
    h = dtau(j)/n(j);
    h_2=h/2; h_6=h/6; h_3=h/3; h_8=h/8;
    % variables used in calculations
    dpsi1=zeros(length(x0),1);dpsi2=dpsi1;dpsi3=dpsi1;dpsi4=dpsi1;
    d_gradient1=zeros(3,1);d_gradient2=d_gradient1;d_gradient3=d_gradient1;d_gradient4=d_gradient1;
    dxp = rhs(t,x(cn(j+1),:)',Theta_0_ht,m,u(j,:)');
    gradient_tmp = [0;0;0];
    for i = cn(j+1) : -1 : cn(j) + 1
        psi_tmp = psi(i,:)';
        [dpsi1, d_gradient1] = rhs_sprzezone(t,psi_tmp,x(i,:)',u(j,:)');
        tmp     = psi_tmp-h_2*dpsi1;
        dxl     = rhs(t,x(i-1,:)',Theta_0_ht,m,u(j,:)');
        xp      = (x(i-1,:)'+x(i,:)')/2 + (dxl-dxp)*h_8;
        dxp     = dxl;
        [dpsi2, d_gradient2] = rhs_sprzezone(t,tmp,xp,u(j,:)');
        tmp     = psi_tmp-h_2*dpsi2;
        [dpsi3, d_gradient3] = rhs_sprzezone(t,tmp,xp,u(j,:)');
        tmp     = psi_tmp-h*dpsi3;
        [dpsi4, d_gradient4] = rhs_sprzezone(t,tmp,x(i-1,:)',u(j,:)'); 
        psi(i-1,:) = psi_tmp - h_3*(dpsi2+dpsi3) - h_6*(dpsi1+dpsi4);
        gradient_tmp = gradient_tmp - h_3*(d_gradient2+d_gradient3) - h_6*(d_gradient1+d_gradient4);
        
    end
    gradient(j,:)=gradient_tmp;
end

end

