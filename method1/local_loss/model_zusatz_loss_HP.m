function [cost,y3,y2] = model_zusatz_loss_HP(xx,print_index)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3 = 0;                     % no constraints exist
    y2 = 0;                     % no constraints exist
    a5_1 = xx(1);               % 
    a5_2 = xx(2);               % 
    mu_0 = 4*pi*10^-7 ; 

    % declare the global variables
    
    global f_min_a5; 
    global f_max_a5; 
    global J_min_a5; 
    global J_max_a5;     
    % basis depth
    global depth_basis_local_loss_parameter
    %
    global cut_edge_model    
    global width_names
    global width_vector
    global frequency_names
    global frequency_vector
    global A_parameters
    global struct_name
    eval(['global ' struct_name]) % measured data will be changed to be globals


    % preparation of data
    number_of_width = length(width_vector);
    % only part of frequency measurement are used to fit data
    index_frequency = (frequency_vector >= f_min_a5 & frequency_vector <= f_max_a5);  % filter by frequency
    frequency_for_a5 = frequency_vector(index_frequency); % only these frequencies are used to fit a_5 parameters
    frequency_name_filtered = frequency_names(index_frequency);
    
    % initialisation of cost function
    cost = 0; % fitness function

    % loops for frequency
    for k = 1: length(frequency_name_filtered) 
        frequency_name = frequency_name_filtered{k};
        frequency_tmp = frequency_for_a5(k);
        
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
        % a5(x) to do if the same depth is used in a5_x
        a5_u = A_parameters.(width_names{end}).a5   ;
        a5_x = @ (x) a5_u.*(1+a5_1.*exp(-a5_2.*x/depth_basis_local_loss_parameter));   
        
        % loops for width
        for i = 1: number_of_width-1 % 120mm is excluded

            width_tmp = width_vector(i); 
            number_of_measurement = eval(struct_name);
            number_of_measurement = length(number_of_measurement.(width_names{i}).(frequency_name).Jmax_T ); 
            a5_messen = A_parameters.(width_names{i}).a5; 
            Excess_Funktion = @(J) a5_messen .* J .^(1.5) .* frequency_tmp.^(1.5) ;
            for j = 1 : number_of_measurement
                J_tmp = eval(struct_name);
                J_tmp =J_tmp.(width_names{i}).(frequency_name).Jmax_T(j);                
                if (J_tmp < J_min_a5) || (J_tmp > J_max_a5)
                    continue;
                end              
                P_excess_measurement = Excess_Funktion(J_tmp); 
                H_tmp = eval(struct_name);
                H_tmp = H_tmp.(width_names{i}).(frequency_name).Hmax_A_m(j);
                mu_relative_undamaged =  interp1(H_no_degration_interp_vector,mu_no_degradation_interp_vector,H_tmp, 'pchip','extrap'); % mu_r_j at H_j
                
                if H_tmp < min(H_no_degration_interp_vector) || H_tmp > max(H_no_degration_interp_vector)
                    continue;
                end
                
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
                % P_excess_func = @(x) a5_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x)')).^(1.5).* frequency_tmp.^(1.5) ;
                P_excess_func = @(x) a5_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^(1.5).* frequency_tmp.^(1.5) ;

                P_excess_mean = 2*integral(P_excess_func, 0 , width_tmp/2)/width_tmp ;

                cost_tmp = abs(P_excess_measurement - P_excess_mean).^2;
                % update of cost
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
    str = sprintf('a5,1:\n');
    disp(str);
    a5_1 
    str = sprintf('a5,2:\n');
    disp(str);
    a5_2          
end