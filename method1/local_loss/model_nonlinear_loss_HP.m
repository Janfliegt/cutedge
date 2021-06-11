function [cost,y3,y2] = model_nonlinear_loss_HP(xx,print_index)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3 = 0;                     % no constraints exist
    y2 = 0;                     % no constraints exist
    a3_1 = xx(1);               % 
    a3_2 = xx(2);               %    
    a4_1 = xx(3);               % 
    a4_2 = xx(4);               %     
    mu_0 = 4*pi*10^-7 ; 

    % declare the global variables
    
    global f_min_a3_a4; 
    global f_max_a3_a4; 
    global J_min_a3_a4; 
    global J_max_a3_a4;     
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
    index_frequency = (frequency_vector >= f_min_a3_a4 & frequency_vector <= f_max_a3_a4);  % filter by frequency
    frequency_for_a3_a4 = frequency_vector(index_frequency); % only these frequencies are used to fit a_5 parameters
    frequency_name_filtered = frequency_names(index_frequency);
    
    
    % initialisation of cost function
    cost = 0; % fitness function

    % loops for frequency
    for k = 1: length(frequency_name_filtered) 
        frequency_name = frequency_name_filtered{k};
        frequency_tmp = frequency_for_a3_a4(k);
        
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
        % a3(x) to do if the same depth is used in a3_x and a4_x
        a3_u = A_parameters.(width_names{end}).a3   ;
        a3_x = @ (x) a3_u.*(1+a3_1.*exp(-a3_2.*x/depth_basis_local_loss_parameter));

        % a4(x)
        a4_u = A_parameters.(width_names{end}).a4   ;
        a4_x = @ (x) a4_u.*(1-a4_1.*exp(-a4_2.*x/depth_basis_local_loss_parameter));        
        
        % loops for width
        for i = 1: number_of_width-1 % 120mm is excluded
            width_tmp = width_vector(i); 
            number_of_measurement=eval(struct_name);
            number_of_measurement = length(number_of_measurement.(width_names{i}).(frequency_name).Jmax_T ); 
            a3_messen = A_parameters.(width_names{i}).a3 ; 
            a4_messen = A_parameters.(width_names{i}).a4 ;      
            a2_messen =  A_parameters.(width_names{i}).a2 ;
            nonlinear_loss_Funktion = @(J) a2_messen.*a3_messen.* (J.^(2+a4_messen)) .* frequency_tmp.^2 ; 
            for j = 1 : number_of_measurement
                J_tmp = eval(struct_name);
               J_tmp=J_tmp.(width_names{i}).(frequency_name).Jmax_T(j);                
                if (J_tmp < J_min_a3_a4) || (J_tmp > J_max_a3_a4)
                    continue;
                end              
                P_nonlinear_measurement = nonlinear_loss_Funktion(J_tmp); 
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
               % P_nonlinear_func = @(x) a2_messen.* a3_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x)')).^(2+a4_x(x)).* frequency_tmp.^(2) ;
                P_nonlinear_func = @(x) a2_messen.* a3_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^(2+a4_x(x)).* frequency_tmp.^(2) ;
                P_nonlinear_mean = 2*integral(P_nonlinear_func, 0 , width_tmp/2)/width_tmp ;

                cost_tmp = abs(P_nonlinear_measurement - P_nonlinear_mean).^2;
                % update of cost function
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
    str = sprintf('a3,1:\n');
    disp(str);
    a3_1 
    str = sprintf('a3,2:\n');
    disp(str);
    a3_2 
    str = sprintf('a4,1:\n');
    disp(str);
    a4_1 
    str = sprintf('a4,2:\n');
    disp(str);
    a4_2       
end