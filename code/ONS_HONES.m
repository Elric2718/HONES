function output = ONS_HONES(x1, data, f, df, T, delta, beta_ons, eps, init_eps)
% This function returns the result at each step from
    % HONES algorithm
    %  Args: 
    %   x1: initial point
    %   data: the data for ONS
    %   f: the loss/payoff function
    %   df: the gradient of f
    %   T: steps to run
    %   delta, beta_ons: parameters in ONS. See (Hazan 2007).
    %   eps: an tolerance to avoid problems induced by numerical accuracy
    %   init_eps: the scale of the initial matrix A0.   
    %  Returns:
    %   x: a vector of size m, the projected restuls at every step t
    %   turn_count: a vetor of size m, the number of turning points at every step t
    %   time: a vector of size m, the time taken at every step t
    %   S: the final support
    %   S_set: a struct of size m, the numbers of non-zero entries in the support at every step t
    %   S_star: the union of the supports of all the steps
    %   nu: the vector containing both nonzero x (S) and nonzero mu (Sc)
    %   mu0: a value in KKT
    %   sparsity: a vector, the size of support at each step t
    %   func_val: a vector of size m, the objection function value at every step t
    %%%%%%%% Initialization %%%%%%%%
    d = size(data, 1);
    T = min(T, size(data, 2));
    time = [];
    turn_count_A = [];
    turn_count_r = [];
    S_set = [];
    sparsity = [];

    % A
    A0 = diag(ones(1,d));
    A = init_eps * A0;
    
    % basic parameters
    func_val = f(x1, data, 1);
    g = df(x1, data, 1);
    r_equiv = zeros(d, 1);
    y_equiv = zeros(d,1);
           
    % support
    S = find(x1);
    Sc = setdiff(1:d, S);
    len_S = length(S);
    len_Sc = d - len_S;
    
    % parameters
    M = zeros(d, d);
    M(S, S) = inv(A(S, S));
    M(Sc, S) = -A(Sc, S)/A(S,S);
    
    eta = zeros(d, 1);
    eta_tilde = zeros(d, 1);
    eta(S) = A(S, S) \ g(S);
    eta(Sc) = g(Sc) + M(Sc, S) * g(S);
    eta_tilde(S) = sum(M(S, S), 2);
    eta_tilde(Sc) = ones(len_Sc, 1) + sum(M(Sc, S), 2);
    
    D = sum(eta_tilde(S));
    Dg = sum(eta(S));
    Dgg = g(S)' * eta(S);
    Dgy = -eta(S)' * r_equiv(S);
    
    % x, mu, mu0, nu
    x = x1;
    mu0 = (1 - sum(M(S, S) * r_equiv(S)))/D;
    
    mu = zeros(d, 1);
    mu(Sc) = -mu0 * (ones(len_Sc, 1) + sum(M(Sc, S),2)) - (r_equiv(Sc) + M(Sc, S) * r_equiv(S));
    
    nu = zeros(d, 1);
    nu(S) = x1(S);
    nu(Sc) = -mu(Sc);
    %%%%%%%% Initialization %%%%%%%%
    

    t = 1;
    while true
        %%%%%%%%%%%%%%%%%%%%%%%%%%% HONES_UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% DirectUpdate %%%%%%%%%
        if t > 1
            Sc = setdiff(1:d, S);
            eta(S) = M(S, S) * g(S);
            eta(Sc) = g(Sc) + M(Sc, S) * g(S);
     
            Dg = sum(eta(S));
            Dgg = eta(S)' * g(S);
            Dgy = eta(Sc)' * y_equiv(Sc) - eta(S)' * r_equiv(S);
        end
        %%%%%%%%% DirectUpdate %%%%%%%%%
        
        %%%%%%%%% Update A part %%%%%%%%%
        lambda = 0;
        count1 = 0;        
        tic;          
        while lambda < 1 - eps
            count1 = count1 + 1;
            upper_increment = 1 - lambda;
        
            %%% FIND_LAMBDA %%%
            rm_entry = [];
            add_entry = [];
            %S_new = S;
    
            u = (Dg * mu0 - Dgy) * (Dg * eta_tilde - D * eta);
            alpha_vec = (D * nu)./(Dg^2 * nu - u);
            index1 = find(alpha_vec >= eps);
            [~, index2] = sort(alpha_vec(index1));
            index2 = index1(index2);
    
    
            if isempty(index2)
                lambda_increment = upper_increment;
            else 
                alpha = alpha_vec(index2(1));
                lambda_increment = alpha/(1 - alpha * Dgg);
        
                if lambda_increment < eps || lambda_increment > upper_increment + eps
                    lambda_increment = upper_increment;
                else
                    %error('debug');
                    index3 = index2((alpha_vec(index2) <= alpha + eps));
                    rm_entry = index3(ismember(index3, S));
                    add_entry = setdiff(index3, rm_entry);

                    %S_new = sort(union(setdiff(S, rm_entry), add_entry));
                end
            end
                        
            lambda = lambda + lambda_increment;
            %disp(lambda)
            %%% FIND_LAMBDA %%%
            
            %%% UPDATE_BY_LAMBDA %%%
            u = (Dg * mu0 - Dgy) * (Dg * eta_tilde - D * eta); % redundant calculation in FIND_LAMBDA, 0.006s
            alpha0 = 1/(1 + lambda_increment * Dgg);
            alpha = lambda_increment * alpha0;
            alpha_tilde = alpha/(D - alpha * Dg^2);%0.017s
            nu = nu + alpha_tilde * u;%0.0064s
            if ~isempty([rm_entry; add_entry])
                nu([rm_entry; add_entry]) = 0;
            end
            mu0 = mu0 + alpha_tilde * (Dg * mu0 - Dgy) * Dg;%0.005s
            M(:, S) = M(:, S) - (alpha * eta) * eta(S)'; %0.05s
            eta_tilde = eta_tilde - alpha * Dg * eta;
            eta = alpha0 * eta;
            D = D - alpha * Dg^2;
            Dg = alpha0 * Dg;
            Dgg = alpha0 * Dgg;
            Dgy = alpha0 * Dgy;%0.026s
            %%% UPDATE_BY_LAMBDA %%%                        
        
            %%% UPDATE_SHRINK_SUPPORT %%%
            while ~isempty(rm_entry)
                new_entry = rm_entry(1);
                Sc = setdiff(1:d, S);
                S_new = setdiff(S, new_entry);
                %Sc_new = setdiff(1:d, S_new);
                Sc_new = [new_entry Sc];
            
                beta = zeros(d, 1);
                beta_tilde = zeros(d, 1);
                beta(S_new) = M(new_entry, S_new)';
                beta_tilde(S_new) = beta(S_new);
                beta(Sc_new) = [-1; M(Sc, new_entry)];
            
                c = y_equiv(Sc_new)' * beta(Sc_new) - r_equiv(S_new)' * beta(S_new) -r_equiv(new_entry) * M(new_entry, new_entry);%-y(new_entry) + M(Sc, new_entry)' * y(Sc);
                Mjj = M(new_entry, new_entry);

            
                eta_j = eta(new_entry);
                eta(new_entry) = 0;
                eta = eta - eta_j/Mjj * beta;
            
                eta_tilde_j = eta_tilde(new_entry);
                eta_tilde(new_entry) = 0;
                eta_tilde = eta_tilde - eta_tilde_j/Mjj * beta;
            
                D = D - eta_tilde_j^2/Mjj;
                Dg = Dg - eta_j * eta_tilde_j/Mjj;
                Dgg = Dgg - eta_j^2/Mjj;
                Dgy = Dgy - eta_j * c/Mjj;

                M(new_entry, :) = 0;
                M(:, S_new) = M(:, S_new) - 1/Mjj * beta * beta_tilde(S_new)';
                M(:, new_entry) = 0;
            
                S = S_new;
                rm_entry = setdiff(rm_entry, new_entry);        
            end            
            %%% UPDATE_SHRINK_SUPPORT %%%    
            
        
            %%% UPDATE_EXPAND_SUPPORT %%%
            while ~isempty(add_entry)
                new_entry = add_entry(1);
                S_new = [S; new_entry];
                Sc_new = setdiff(1:d, S_new);
                [~, index_order] = sort(S_new);% need to think more abou the order
            
                eta_j = eta(new_entry);
                eta_tilde_j = eta_tilde(new_entry);
                Ajj_tilde = A(new_entry, new_entry) + M(new_entry, S) * A(S, new_entry) + lambda * g(new_entry) * eta_j;
            
                gamma = zeros(d, 1);
                gamma_tilde = zeros(d, 1);
                gamma(S_new) = [M(new_entry, S)'; 1];
                gamma_tilde(S_new) = gamma(S_new);
                gamma(Sc_new) = -A(Sc_new, new_entry) - A(Sc_new, S) * M(new_entry, S)' - lambda * eta_j * g(Sc_new);
            
                b = y_equiv(Sc_new)' * gamma(Sc_new) - r_equiv(S_new)' * gamma(S_new);
           
            
                D = D + eta_tilde_j^2/Ajj_tilde;
                Dg = Dg + eta_j * eta_tilde_j/Ajj_tilde;
                Dgg = Dgg + eta_j^2/Ajj_tilde;
                Dgy = Dgy - y_equiv(new_entry) * eta_j + eta_j * b/Ajj_tilde;
       
                M(new_entry, :) = 0;
                M(:, S_new) = M(:, S_new) + gamma * gamma_tilde(S_new)'/Ajj_tilde; 
           
                eta(new_entry) = 0;
                eta = eta + eta_j/Ajj_tilde * gamma;
           
                eta_tilde(new_entry) = 0;
                eta_tilde = eta_tilde + eta_tilde_j/Ajj_tilde * gamma;
            
                S = S_new(index_order);
                add_entry = setdiff(add_entry, new_entry);
            end 
            %%% UPDATE_EXPAND_SUPPORT %%%        
        end  
        time_recorder1 = toc;
        %%%%%%%%% Update A part %%%%%%%%%
        
        
        %%%%%%%%% update A, l, r_euiqv %%%%%%%%%
        A = A + g * g';
        
        %% ONS
        %l = A * x + delta * (1 + 1/beta_ons) * g - r_equiv;
        %% ONS variant        
        l = delta * (1 + 1/beta_ons) * g;
        
        r_equiv = l + r_equiv;
        %%%%%%%%% update A, l, r_euiqv %%%%%%%%%
   
        %%%%%%%%% DirectUtildeUpdate %%%%%%%%%
        Sc = setdiff(1:d, S);
     
        xi = zeros(d, 1);
        xi(S) = -M(S, S) * l(S);
        xi(Sc) = -l(Sc) - M(Sc, S) * l(S);
     
        Dl = sum(xi(S));
        %%%%%%%%% DirectUtildeUpdate %%%%%%%%%
   
        %%%%%%%%% Update r part %%%%%%%%%
        lambda_utilde = 0;
        count2 = 0;
        
        tic;  
        while lambda_utilde < 1 - eps
            count2 = count2 + 1;
            upper_increment = 1 - lambda_utilde;
       
            %%% FIND_UTILDE_LAMBDA %%%
            denominator = xi - Dl / D * eta_tilde;
            nu_xi = nu ./ denominator;
            index1 = find(nu_xi >= eps);
            [~, index2] = sort(nu_xi(index1));
            index2 = index1(index2);
    
            rm_entry = [];
            add_entry = [];
            %S_new = S;  
    
            if isempty(index2)
                lambda_utilde_increment = upper_increment;
            else 
                lambda_utilde_increment = nu_xi(index2(1));
                if lambda_utilde_increment > upper_increment + eps
                    lambda_utilde_increment = upper_increment;
                else
                %error('debug');
                index3 = index2((nu_xi(index2) <= lambda_utilde_increment + eps));
                rm_entry = index3(ismember(index3, S));
                add_entry = setdiff(index3, rm_entry);

                %S_new = union(setdiff(S, rm_entry), add_entry);
                end
            end    
            lambda_utilde = lambda_utilde + lambda_utilde_increment;
            %disp(lambda_utilde)
            %%% FIND_UTILDE_LAMBDA %%%
   
            %%% UPDATE_BY_UTILDE_LAMBDA %%%
            nu = nu - (xi - Dl / D * eta_tilde) * lambda_utilde_increment;
            if ~isempty([rm_entry; add_entry])
                nu([rm_entry; add_entry]) = 0;
            end
            mu0 = mu0 + Dl/D * lambda_utilde_increment;
            %%% UPDATE_BY_UTILDE_LAMBDA %%%
       
            %%% UPDATE_UTILDE_SHRINK_SUPPORT %%%
            while ~isempty(rm_entry)
                %error('debug');
                new_entry = rm_entry(1);
                Sc = setdiff(1:d, S);
                S_new = setdiff(S, new_entry);
                %Sc_new = setdiff(1:d, S_new);
                Sc_new = [new_entry Sc];
            
                %error('debug');
                beta = zeros(d, 1);
                beta_tilde = zeros(d, 1);
                beta(S_new) = M(new_entry, S_new)';
                beta_tilde(S_new) = beta(S_new);
                beta(Sc_new) = [-1; M(Sc, new_entry)];
            
                Mjj = M(new_entry, new_entry);
                eta_tilde_j = eta_tilde(new_entry);
                eta_tilde(new_entry) = 0;
                eta_tilde = eta_tilde - eta_tilde_j/Mjj * beta;                        

                K_tilde = l(S_new)' * beta(S_new) + l(new_entry) * Mjj;
                N_tilde = sum(beta(S_new));
            
                D = D - eta_tilde_j^2/Mjj;
                Dl = Dl + K_tilde * N_tilde/Mjj - xi(new_entry);
            
                xi(new_entry) = 0;
                xi = xi + K_tilde/Mjj * beta;           
            
                M(new_entry, :) = 0;
                M(:, S_new) = M(:, S_new) - 1/Mjj * beta * beta_tilde(S_new)';
                M(:, new_entry) = 0;           
            
                S = S_new;
                rm_entry = setdiff(rm_entry, new_entry);
            end
            %%% UPDATE_UTILDE_SHRINK_SUPPORT %%%
       
            %%% UPDATE_UTILDE_EXPAND_SUPPORT %%%
            while ~isempty(add_entry)
                new_entry = add_entry(1);
                S_new = [S; new_entry];
                Sc_new = setdiff(1:d, S_new);
                [~, index_order] = sort(S_new);% need to think more abou the order
                  
                Ajj_tilde = A(new_entry, new_entry) + M(new_entry, S) * A(S, new_entry);
                gamma = zeros(d, 1);
                gamma_tilde = zeros(d, 1);
                gamma(S_new) = [M(new_entry, S)'; 1];
                gamma(Sc_new) = -A(Sc_new, new_entry) - A(Sc_new, S) * M(new_entry, S)';
                gamma_tilde(S_new) = gamma(S_new);
                      
                eta_tilde_j = eta_tilde(new_entry);
                eta_tilde(new_entry) = 0;
                eta_tilde = eta_tilde + eta_tilde_j/Ajj_tilde * gamma;
            
                K = l(S_new)' * gamma(S_new);
                N = sum(gamma(S_new));
       
                D = D + N^2/Ajj_tilde;
                Dl = Dl - K * N/Ajj_tilde;
            
                xi(new_entry) = 0;
                xi = xi - K/Ajj_tilde * gamma;
            
                M(new_entry, :) = 0;
                M(:, S_new) = M(:, S_new) + gamma * gamma_tilde(S_new)'/Ajj_tilde; 
            
                S = S_new(index_order);
                add_entry = setdiff(add_entry, new_entry);
            end 
            %%% UPDATE_UTILDE_EXPAND_SUPPORT %%%
               
        end
        time_recorder2 = toc;
        %%%%%%%%% Update r part %%%%%%%%%   
        t = t + 1;   
        
        S_set(t - 1).S = S;
        sparsity = [sparsity, length(S)];
        %time = [time; time_recorder];
        time = [time; time_recorder1 + time_recorder2];
        turn_count_A = [turn_count_A count1];
        turn_count_r = [turn_count_r count2];
        disp([t count1 count2]); 
        if t > T
            break;
        end
        
        x = zeros(d, 1);
        x(S) = ProjToSimp(nu(S));
        g = df(x, data, t);                      
        func_val = [func_val, f(x, data, t)]; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%% HONES_UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    nu(S) = ProjToSimp(nu(S));
    x = zeros(d, 1);
    x(S) = nu(S);
    output = struct('S', S, 'nu', nu, 'mu0', mu0, 'x', x, 'turn_count_A', turn_count_A, 'turn_count_r', turn_count_r, 'time', time, 'S_set', S_set, 'sparsity', sparsity, 'func_val', func_val);
end