function [ optimal_control ] = bfgs(u,epsilon,lb,ub,maxIter)
    
    % brief - zak³adamy, ¿e u jest wektorem o jednej kolumnie
    % initialize parameters
    init;
    
    % reneval parameter
    Reneval = 1;
    contraction = 0.5;
    
    % first gradient calculation
    [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u,nodes);
    Q = x(end,end);
    
    % algorithm iterations
    for i = 1 : maxIter
        
        disp('Iteration - ');disp(i);
        disp(Q);
        % check stop condition
        if(norm(gradient) < epsilon) 
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
        step = 1;
        for j = 1 : 50
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
            gradient_prev = gradient;
            [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u_temp,nodes);
            Q_temp = x(end,end)
            
            % check stop conditions
            if Q_temp < Q
                Q = Q_temp;
                Reneval = 0;
                s = u_temp - u;
                r = gradient - gradient_prev;
                u = u_temp;
                break;
            else
                if step > 1e-15
                    step = step * contraction;
                else
                    Reneval = 1;  
                    break;
                end
            end   
        end
    end

end

