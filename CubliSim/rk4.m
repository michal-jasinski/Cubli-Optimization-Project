function [t,x] = rk4(rhs,x0,tf,T,Theta_0_ht,m,alpha,beta,gamma,delta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   @params:
%   rhs : right side of the differential equation
%   x0  : initial condition vector
%   tf  : duration of the simulation
%   I, m, g, T : params used in rhs 


nt = length(T);     % length of the control vector, number of nodes 
n  = length(x0);    % size of the initial condition vector
h  = tf/nt;         % step width (sample time)       
x  = zeros(nt,n);   % matrix for the data (nt x n)

tmp  = zeros(n,1);  % (3 x 1) vector of zeros for calculations
xtmp = x0;
x(1,:) = x0';       % first element in the data matrix -> x0
t = 0;

% variables used in calculations
dx1=zeros(n,1);dx2=dx1;dx3=dx1;dx4=dx1;
h_2=h/2; h_6=h/6; h_26=2*h_6;

for i = 1 : nt
    dx1  = rhs(t,xtmp,Theta_0_ht,m,alpha,beta,gamma,delta);
    tmp  = xtmp+h_2*dx1; t = t+h_2;
    dx2  = rhs(t,tmp,Theta_0_ht,m,alpha,beta,gamma,delta);
	tmp  = xtmp+h_2*dx2;
	dx3  = rhs(t,tmp,Theta_0_ht,m,alpha,beta,gamma,delta);
	tmp  = xtmp+h*dx3; t = t+h_2;
	dx4  = rhs(t,tmp,Theta_0_ht,m,alpha,beta,gamma,delta);
 	xtmp = xtmp+h_6*(dx1+dx4)+h_26*(dx2+dx3);
   	x(i,:)=xtmp';
end

t = linspace(0,tf,nt)';

end

