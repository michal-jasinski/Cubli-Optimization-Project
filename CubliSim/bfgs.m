function [ Q, gradient] = bfgs(u)

    init;
    [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u,nodes);
    Q = x(end,end);
  
end

