function [cost,y3,y2] = model_hyst_HP(xx,print_index)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3 = 0;                     % no constraints exist
    y2 = 0;                     % no constraints exist
    a1_1 = xx(1);               % 
    a1_2 = xx(2);               % 
    alpha1 = xx(3);             %
    alpha2 = xx(4);             %
    mu_0 = 4*pi*10^-7 ; 
    % range of frequency
    global f_min_a1; 
    global f_max_a1; 
    global J_min_a1; 
    global J_max_a1; 
    
    % basis depth
    global depth_basis_local_loss_parameter
    % declare the global variables
    global cut_edge_model 
    global width_names
    global width_vector
    global frequency_names
    global frequency_vector
    global A_parameters
    global struct_name
    eval(['global ' struct_name]) % measured data will be changed to be globals


    % the smallest frequency will be used to calculate the local hysteresis
    % parameter
    number_of_width = length(width_vector);
    
    % only part of frequency measurement are used to fit data
    index_frequency = (frequency_vector >= f_min_a1 & frequency_vector <= f_max_a1);  % filter by frequency
    frequency_for_a1 = frequency_vector(index_frequency); % only these frequencies are used to fit a_5 parameters
    frequency_name_filtered = frequency_names(index_frequency);      

    % initialisation of cost function
    cost = 0; % fitness function
    
    % loops for frequency
    for k = 1: length(frequency_name_filtered) 
        frequency_name = frequency_name_filtered{k};
        frequency_tmp = frequency_for_a1(k);
        
        % mu vs H of undamaged material
        H_no_degration_interp_vector = cut_edge_model.(frequency_name).H_no_degreation_pchip; 
        mu_no_degradation_interp_vector = cut_edge_model.(frequency_name).mu_no_degreation_pchip; 
        
        % delta_mu_c(H)   
        delta_mu_interp_vector = cut_edge_model.(frequency_name).delta_mu_interp_pchip_final  ;
        delta_H_interp_vector = cut_edge_model.(frequency_name).H_interp_pchip_final  ; 
        
        % eta(x) import, there is only more than frequency used for excess loss, and
        % therefore more than eta(x) and more than one Einflusstiefe
        eta_x = cut_edge_model.(frequency_name).eta.fitting_function;
        parameters = coeffvalues(eta_x);
        depth = 1/ abs(parameters(2));   
        eta_x = @(x) exp(-x/depth);

        % definition of local loss parameter functions  
        % a1(x) to do if the same depth is used in a1_x and alpha_x
        
        a1_u = A_parameters.(width_names{end}).a1;
        a1_x = @ (x) a1_u.*(1+a1_1.*exp(-a1_2.*x/depth_basis_local_loss_parameter));
        % alpha(x)
        alpha_u =  A_parameters.(width_names{end}).alpha; 
        alpha_x = @ (x) alpha_u .* (1-alpha1.* exp(-alpha2.*x/depth_basis_local_loss_parameter)) ;      
        for i = 1: number_of_width-1 % 120mm is excluded

            width_tmp = width_vector(i); 
            Jmax_T_cur = eval(struct_name);
            Jmax_T_cur= Jmax_T_cur.(width_names{i}).(frequency_name).Jmax_T;
            number_of_measurement = length(Jmax_T_cur); 

            a1_messen = A_parameters.(width_names{i}).a1 ; 
            alpha_messen = A_parameters.(width_names{i}).alpha ; 
            HystFunktion = @(J) a1_messen .* J .^(alpha_messen) .* frequency_tmp;
            for j = 1 : number_of_measurement
                J_tmp = eval(struct_name);
                J_tmp = J_tmp.(width_names{i}).(frequency_name).Jmax_T(j);
                if (J_tmp < J_min_a1) || (J_tmp > J_max_a1)
                    continue;
                end     
                
                P_hyst_measurement = HystFunktion(J_tmp); 
                H_tmp = eval(struct_name);
                H_tmp = H_tmp.(width_names{i}).(frequency_name).Hmax_A_m(j);

                if H_tmp < min(H_no_degration_interp_vector) || H_tmp > max(H_no_degration_interp_vector)
                    continue;
                end

                mu_relative_undamaged =  interp1(H_no_degration_interp_vector,mu_no_degradation_interp_vector,H_tmp, 'pchip','extrap'); % mu_r_j at H_j

                if isnan(mu_relative_undamaged)
                    disp('nan')
                end

                if mu_relative_undamaged < 0 
                     disp('negative')
                end
                delta_mu_c = interp1(delta_H_interp_vector,delta_mu_interp_vector,H_tmp, 'pchip','extrap'); 
                if isnan(mu_relative_undamaged)
                    disp('nan')
                end
                if delta_mu_c < 0 
                     disp('negative')
                end
                P_hyst_func = @(x) a1_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^alpha_x(x).*frequency_tmp;

                P_hyst_mean = 2*integral(P_hyst_func, 0 , width_tmp/2)/width_tmp; 

                cost_tmp = abs(P_hyst_measurement - P_hyst_mean).^2;
                cost = cost + cost_tmp; 

            end

        end
    
    end  

    % constraint definition
%     if max(max(f_H_b_matrix)) > 1
%         y3 = 1; 
%         y2 =  max(max(f_H_b_matrix)) -1; 
%     end
    
    if print_index == 0
        return;
    end
    str = sprintf('RESULTS:\n');
    disp(str);
    str = sprintf('cost function:\n');
    disp(str);
    cost 
    str = sprintf('a1,1:\n');
    disp(str);
    a1_1 
    str = sprintf('a1,2:\n');
    disp(str);
    a1_2     
    str = sprintf('alpha1:\n');
    disp(str);
    alpha1 
    str = sprintf('alpha2:\n');
    disp(str);
    alpha2       
end