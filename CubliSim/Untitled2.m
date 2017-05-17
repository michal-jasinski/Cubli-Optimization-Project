clear all;

load('optimal_control');

% czas
T = 10;

% okres próbkowania
h0 = 0.005;

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

x = zeros(cn(end),3);   % matrix for the data (nt x n)
t = zeros(cn(end),1);



