function [cost,y3,y2] = modeling_corrected_data(xx,print_DE)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3 = 0;     % no constraints exist
    y2 = 0;     % no constraints exist
    F_b = xx;               % from 4,5, ..., 60mm
    % declare the global variables
    global f_in_Hz         
    global width_names
    global width_vector
    global H_interpolation
    global mu_no_degradation
    
    % define mu without degradation, corresponding to 120mm
    % measurement points in H are used
    mu_no_degradation = eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'mu_relative_interp_2'));
            
    % definition of interpolation points in H
    H_interpolation = eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'H_interp_2'));  % here can be changed, if more interpolated H are desired to approximate the F(b)

    % length 
    number_width = length(width_names);
    number_interpolated_H = length(H_interpolation);
    % define matrix for f(H,bi)
    global f_H_b_matrix
    f_H_b_matrix = zeros(number_interpolated_H, number_width-1);  % each row for a H, each column corresponds to a width
    % filling f_H_b_matrix
    for k = 1 : number_interpolated_H   % iteration along H
        for g = 1: number_width-1       % iteration along width
            F_b_tmp = F_b(g);
            mu_without_degradation_tmp = mu_no_degradation(k);
            % mu under other width will be interplolated
            mu_tmp = eval(strcat('f_in_Hz', '.', width_names{g}, '.', 'mu_relative_interp_2'));
            H_tmp = eval(strcat('f_in_Hz', '.', width_names{g}, '.', 'H_interp_2')); 
            mu_with_degradation_vector = interp1(H_tmp, mu_tmp, H_interpolation, 'pchip');
            % mu under the same H is saved in mu_measure_with_degradation_tmp
            mu_measure_with_degradation_tmp = mu_with_degradation_vector(k); 
            % 
            f_H_b_tmp = (mu_without_degradation_tmp - mu_measure_with_degradation_tmp)/ F_b_tmp / mu_without_degradation_tmp; % normalized, different from silas formula
            %  f_H_b_tmp = (mu_without_degradation_tmp - mu_measure_with_degradation_tmp)/ F_b_tmp ; % non-normalized, the silas formula
            %  f_H_b_tmp = min(f_H_b_tmp,1);                                  % f_H_b_tmp can not over 1
            %  matrix filling
            f_H_b_matrix(k,g) = f_H_b_tmp;
        end
    end
    cost = mean(var(f_H_b_matrix,0,2));
    if max(max(f_H_b_matrix)) > 1  % f_H_b_tmp should not be greater than 1, which is a constraint to be considered
        y3 = 1; 
        y2 =  max(max(f_H_b_matrix)) -1; 
    end
    
    if print_DE == 0
        return;
    end
    str = sprintf('RESULTS:\n');
    disp(str);
    str = sprintf('F_b:\n');
    disp(str);
    F_b;
    str = sprintf('cost function:\n');
    disp(str);
    cost;
    f_H_b_matrix;
    global variance_for_different_H
    variance_for_different_H = var(f_H_b_matrix,0,2);
    % delta_mu_c_H calculate 
    global delta_mu_c_H
    delta_mu_c_H =  zeros(length(H_interpolation),1); 
    % preprocess f_H_b_matrix so that no element bigger than 1 due to numeric error
    f_H_b_matrix(f_H_b_matrix > 1) = 1; 
    %
    for k = 1:length(H_interpolation)
        delta_mu_c_H_tmp = mean(f_H_b_matrix(k,:)) .* mu_no_degradation(k); % average f_H_b used to calculate delta_mu_c_H 
        delta_mu_c_H(k) = delta_mu_c_H_tmp; 
    end
end