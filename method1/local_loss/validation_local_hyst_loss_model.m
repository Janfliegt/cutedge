% valiateion of the local iron loss parameters determined by the different evoltuon
function [] = validation_local_hyst_loss_model(processed_data_dir, validate_width, frequency_max, Polarisation_J_max, main_dir)

% Polarisation_J_max: validation range 
%% import of local loss model developed without einzelne width
local_loss_paras_dir = strcat('local_loss_parameters','_without_', num2str(validate_width),'mm');
path_of_local_para_base = fullfile(processed_data_dir, local_loss_paras_dir, 'hyst'); 

path_of_cut_edge_model =fullfile(processed_data_dir, 'plot_cut_edge_model_cor_data', 'cut_edge_model.mat');
path_of_A_parameters = fullfile(processed_data_dir, 'A_parameters','A_parameters.mat');
load(fullfile(processed_data_dir, 'simplified_data.mat')); 
load(path_of_cut_edge_model);
load(path_of_A_parameters)
%%
mu_0 = 4*pi*10^-7 ;
%%
width_name = strcat('B', num2str(validate_width), 'mm');
save_dir = fullfile(strcat('validate_result','_', width_name), 'hyst');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

struct_name = 'L_messung_simplified';
eval(['global ' struct_name]) % measured data will be changed to be globals
% width_names = fieldnames(eval(struct_name));    % width names 

%%
width_struct = strcat(struct_name, '.', width_name); 
global frequency_names
frequency_names = fieldnames(eval(width_struct));       % measured frequencies under each width
number_frequency = length(frequency_names); 
frequency_vector = zeros(number_frequency,1);
num_frequency_used = 0;
for j = 1 : number_frequency
    frequency_name_txt =  frequency_names{j};
    frequency_tmp = regexp(frequency_name_txt, '\d*\d*', 'match');
    if length(frequency_tmp) == 1
        frequency_tmp = frequency_tmp{:};
        frequency_tmp = str2double(frequency_tmp);
        if frequency_tmp <= frequency_max
            frequency_vector(j) = frequency_tmp;
            num_frequency_used = num_frequency_used + 1;
        end
    else       
        frequency_tmp= strcat(frequency_tmp{1}, '.', frequency_tmp{end});
        frequency_tmp = str2double(frequency_tmp);
%         frequency_vector(j) = frequency_tmp ; %  width vector are saved in width_vector
        if frequency_tmp <= frequency_max
            frequency_vector(j) = frequency_tmp;
            num_frequency_used = num_frequency_used+ 1;
        end
    end      

end
number_frequency = num_frequency_used;
frequency_vector = frequency_vector(1:num_frequency_used); % the frequencies based on which the loss model will be validated

%% comparison between measurement and calculated loss by using the cut edge model
P_hyst_mean = [];
P_hyst_measurement = [];

for index_frequency = 1 : number_frequency  
        P_hyst_mean_cur_fz = [];
        P_hyst_measurement_fz = [];
        Jmax_T_cur_plot = []; 
        P_difference = [];
        frequency_name = frequency_names{index_frequency};
        path_of_alpha = fullfile(path_of_local_para_base, frequency_name, 'alpha_x', 'alpha_distance.mat');
        path_of_a1 = fullfile(path_of_local_para_base, frequency_name, 'a1_x','a1_distance.mat');
        load(path_of_alpha)
        load(path_of_a1)
        a1_x = a1_distance.a1_x;
        alpha_x = alpha_distance.alpha_x;
        
        H_no_degration_interp_vector = cut_edge_model.(frequency_name).H_no_degreation_pchip; 
        mu_no_degradation_interp_vector = cut_edge_model.(frequency_name).mu_no_degreation_pchip; 
        
        % delta_mu_c(H)   
        delta_mu_interp_vector = cut_edge_model.(frequency_name).delta_mu_interp_pchip_final  ;
        delta_H_interp_vector = cut_edge_model.(frequency_name).H_interp_pchip_final  ; 
        
        % eta(x) import, there is only more than frequency used for zusatz loss, and
        % therefore more than eta(x) and more than one Einflusstiefe
        eta_x = cut_edge_model.(frequency_name).eta.fitting_function;
        parameters = coeffvalues(eta_x);
        depth = 1/ abs(parameters(2));   
        eta_x = @(x) exp(-x/depth);  
        L_messung = eval(struct_name);
        Jmax_T_cur= L_messung.(width_name).(frequency_name).Jmax_T;
        Ps_W_kg_cur = L_messung.(width_name).(frequency_name).Ps_W_kg;
        number_of_measurement = length(Jmax_T_cur); 

        a1_messen = A_parameters.(width_name).a1 ; 
        alpha_messen = A_parameters.(width_name).alpha ; 
        frequency_tmp = frequency_vector(index_frequency); 
        HystFunktion = @(J) a1_messen .* J .^(alpha_messen) .* frequency_tmp;
        for j = 1 : number_of_measurement
            J_tmp = Jmax_T_cur(j);  
            if J_tmp >= Polarisation_J_max % validation in low J area. to do
                continue; 
            end
            P_hyst_measurement_cur = HystFunktion(J_tmp); 
            H_tmp = L_messung.(width_name).(frequency_name).Hmax_A_m(j);
            
            if (H_tmp > max(H_no_degration_interp_vector) || H_tmp < min(H_no_degration_interp_vector))
                continue; 
            end
 
            if (H_tmp > max(delta_H_interp_vector) || H_tmp < min(delta_H_interp_vector))
                continue; 
            end
            
            mu_relative_undamaged =  interp1(H_no_degration_interp_vector,mu_no_degradation_interp_vector,H_tmp, 'pchip','extrap'); % mu_r_j at H_j
            delta_mu_c = interp1(delta_H_interp_vector,delta_mu_interp_vector,H_tmp, 'pchip','extrap'); 
            % loss local calculatetion
            P_hyst_func = @(x) a1_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^alpha_x(x).*frequency_tmp;
            % integration to calculate average
            P_hyst_mean_cur = 2*integral(P_hyst_func, 0 , validate_width/2)/validate_width; 
            % update of the loss for calculated and measrued 
            P_hyst_mean_cur_fz = [P_hyst_mean_cur_fz P_hyst_mean_cur];
            P_hyst_measurement_fz = [P_hyst_measurement_fz, P_hyst_measurement_cur];
            % difference between measured values and calculated values in
            % square of the value difference
            P_difference_tmp = (P_hyst_measurement_cur - P_hyst_mean_cur);
            % update of the variablel P_difference
            P_difference = [P_difference P_difference_tmp]; 
            
            % update of the Jmax_T_cur_plot for the plots later
            Jmax_T_cur_plot = [Jmax_T_cur_plot J_tmp]; 
        end  
        figure
        % plot(Jmax_T_cur, P_hyst_measurement_fz, Jmax_T_cur, P_hyst_mean_cur_fz); % by Li
        % plot(Jmax_T_cur_plot, P_hyst_measurement_fz, Jmax_T_cur_plot, P_hyst_mean_cur_fz); % by Li
        bar(Jmax_T_cur_plot , [P_hyst_measurement_fz ; P_hyst_mean_cur_fz  ; P_difference])
        xlabel('$J_{\rm max}$ in T',  'Interpreter', 'latex', 'Fontsize', 16)
        ylabel('$P_{\rm hyst-mean}$ in W/kg',  'Interpreter', 'latex', 'Fontsize', 16)
        title(strcat('validation of P_{hys} ', ' of', '--', width_name,  '--under--',  frequency_name), 'Fontsize', 16 )
        legend('measured', 'calculated', 'difference', 'Location', 'Northwest', 'Fontsize', 16)
        
        img = gcf;
        set(gcf,'position',[0.2,0.2,720,405]); 

        save_dir_fz = [save_dir,'\' frequency_name,];
        if ~exist(save_dir_fz,'dir')
            mkdir(save_dir_fz);
        end
        filename = strcat('P_hyst_compare');
        print(img, '-dpng', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.jpg']);
        print(img, '-depsc', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.eps']);
        print(img, '-dmeta', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.emf']); 
        saveas(img,[save_dir_fz,'\',convertStringsToChars(filename),'.fig']);
        pause(2)
        close 
        save_matfile_frequency = strcat(save_dir_fz,'\','P_hyst_compare.mat');
        Phys.measured = P_hyst_measurement_fz; 
        Phys.calculated = P_hyst_mean_cur_fz; 
        Phys.J_measured = Jmax_T_cur_plot; 
        Phys.P_difference = P_difference; 
        save(save_matfile_frequency,'Phys');
        % save(save_matfile_frequency,'P_hyst_mean_cur_fz');  
        % save(save_matfile_frequency,'P_hyst_measurement_fz', '-append');
        
        % why are the following variables needed? saving the data for all
        % frequency?
        % P_hyst_mean = [P_hyst_mean; P_hyst_mean_cur_fz];
        % P_hyst_measurement = [P_hyst_measurement;P_hyst_measurement_fz];
        
        %% update in all frequency
        eval(strcat('allfrequency', '.', frequency_name, '.', 'Phys', '=', 'Phys'))

end
save_matfile_allfrequency =  fullfile(save_dir, 'allfrequency.mat');    
save(save_matfile_allfrequency, 'allfrequency');
close all
%% return to the main dierectory
eval(['cd ' main_dir])
