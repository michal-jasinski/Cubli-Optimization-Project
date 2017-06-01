function [ optimal_control ] = bfgs(u,epsilon,lb,ub,maxIter)
    
    % brief - zak³adamy, ¿e u jest wektorem o jednej kolumnie
    % initialize parameters
    init;
    
    % reneval parameter
    Reneval = 1;
    contraction = 0.1;
    error = 0;
    optimal_control = zeros(size(u));
    
    % first gradient calculation
    [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u,nodes);
    Q = x(end,end);
    
    % algorithm iterations
    for i = 1 : maxIter
        
        fprintf('Iteration - %d  Q = %.10f  norm = %.10f\n',i,Q,norm(gradient));
        % check stop condition
        
        if(norm(gradient) < epsilon || error == 1) 
            break;
        end
        
        if Reneval == 1
            Weight = eye(length(gradient));
        else
            s_t = transpose(s);
            r_t = transpose(r);
            Weight = Weight + (r*r_t)/(s_t*r) - ...
                    (Weight*s*s_t*Weight)/(s_t*Weight*s); 
        end

        direction = -Weight\gradient;
        if transpose(direction)*gradient >= 0
            if Reneval == 1
                % there is no improvement in the direction of the gradient
                break;
            else
                Reneval = 1;
                continue;
            end
        end
        
        % line search
        step = 2;
        for j = 1 : 100
            u_temp = u;
            % get new control with saturation
            for k = 1 : length(u_temp)
                value = u_temp(k) + step*direction(k);
                if value > ub
                    u_temp(k) = ub;
                elseif value < lb
                     u_temp(k) = lb;
                else
                    u_temp(k) = value;
                end
            end
            % calculate quality indicator in this point
            [t,x,psi,gradient_temp] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u_temp,nodes);
            Q_temp = x(end,end);
            
            % check stop conditions
            if Q_temp < Q
                Q = Q_temp;
                Reneval = 0;
                s = u_temp - u;
                r = gradient_temp - gradient;
                u = u_temp;
                gradient = gradient_temp;
                break;
            else
                if step > 1e-15
                    step = step * contraction;
                else
                    Reneval = 1;
                    if i == 1
                        error = 1;
                    end
                    break;
                end
            end   
        end
        Reneval = 1;
    end
    optimal_control = u;
end

