clear all;

load('optimal_control');

% czas
T = 10;

% okres pr�bkowania
h0 = 0.005;

% w�z�y strukturalne
tau = linspace(0,T,11);

% odleg�o�� czasowa pomi�dzy tymi w�z�ami
dtau = diff(tau);
s = size(dtau);

% sterowanie
u = zeros(s(2),3);
for i = 1 : s(2)
%     u(i,:) = control(2,:);
      u(i,:) = 1;
end

% ilo�� pr�bek na ka�dy w�ze� strukturalny
n = ceil(dtau/h0);

% w�z�y dyskretyzacji cn(i) numer w�z�a dyskretyzacji
% b�d�cego i-tym w�z�em strukturalnym
cn = cumsum([1 n]);

x = zeros(cn(end),3);   % matrix for the data (nt x n)
t = zeros(cn(end),1);



