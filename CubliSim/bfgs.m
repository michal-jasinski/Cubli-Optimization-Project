function [ optimal_control ] = bfgs(u,epsilon,lb,ub,maxIter)
    
    % brief - zak³adamy, ¿e u jest wektorem o jednej kolumnie
    % initialize parameters
    init;
    
    u_temp = 0;
    Q_prev = 0;
    gradient_temp = 0;
    
    % reneval parameter
    Reneval = 1;
    contraction = 0.5;
    
    % algorithm iterations
    for i = 1 : maxIter
        
        % simulation
        [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u,nodes);
        Q = x(end,end);

        if(norm(gradient) < epsilon) 
            break;
        end

        if Reneval == 1
            Weight = eye(length(gradient));
        else
            Weight = Weight; % dopytac, jaki wymiar wektora W sterowanie jest (10 na 3)
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
            [t,x,psi,gradient_temp] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u_temp,nodes);
            Q_new = x(end,end);
            
            % check stop conditions
            if Q_new < Q
                Reneval = 0;
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

