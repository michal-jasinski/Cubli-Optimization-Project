function [ optimal_control ] = bfgs(u,epsilon,lb,ub,maxIter)
    
    % brief - zak�adamy, �e u jest wektorem o jednej kolumnie
    % initialize parameters
    init;
    
    u_prev = 0;
    Q_prev = 0;
    gradient_prev = 0;
    
    % reneval parameter
    Renewal = 1;
    step = 1;
    contraction = 0.5;
    
    % algorithm iterations
    for i = 1 : maxIter
        
        % simulation
        [t,x,psi,gradient] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u,nodes);
        Q = x(end,end);

        if(norm(gradient) < epsilon) 
            break;
        end

        if(Renewal == 1)
            Weight = eye(length(gradient));
        else
            Weight = Weight; % dopytac, jaki wymiar wektora W sterowanie jest (10 na 3)
        end

        direction = -Weight\gradient;

        if transpose(d)*gradient >= 0
            if i == 1
                % na kierunku gradientu nie ma poprawy - KONIEC
                break;
            else
                Reneval = 1;
                continue;
            end
        end

        % line search
        for j = 1 : 50
            u_prev = u;
            for k = 1 : length(u)
                value = u_prev(k) + s*d(k);
                if value > ub
                    u(k) = ub;
                elseif value < lb
                     u(k) = lb;
                else
                    u(k) = value;
                end
            end
            
            [t,x,psi,gradient_prev] = rk4(@rhs,@rhs_sprzezone,x0,time,sample_time,Theta_0_ht,m,u,nodes);
            Q_new = x(end,end);
            if(Q_new < Q)
                Reneval = 0;
                u_prev = u;
                continue;
            else
                step = step * contraction;
            end   
        end
    end


end

