% valiateion of the local iron loss parameters determined by the different evoltuon

function [] = validatation_of_local_loss_model_all_losses(processed_data_dir, validate_width, frequency_max, Polarisation_J_max, main_dir)

%% import of local loss model developed without einzelne width
local_loss_paras_dir = strcat('local_loss_parameters','_without_', num2str(validate_width),'mm');
path_of_local_para_base_hyst = fullfile(processed_data_dir, local_loss_paras_dir, 'hyst');
path_of_local_para_base_excess = fullfile(processed_data_dir, local_loss_paras_dir, 'zusatz');
path_of_local_para_base_nonlinear = fullfile(processed_data_dir, local_loss_paras_dir, 'nonlinear');

load(fullfile(processed_data_dir, 'simplified_data.mat')); 
path_of_cut_edge_model =fullfile(processed_data_dir, 'plot_cut_edge_model_cor_data', 'cut_edge_model.mat');
path_of_A_parameters = fullfile(processed_data_dir, 'A_parameters', 'A_parameters.mat');
load(path_of_cut_edge_model);
load(path_of_A_parameters)
%% 
mu_0 = 4*pi*10^-7 ;
width_name = strcat('B', num2str(validate_width), 'mm');
save_dir = fullfile(processed_data_dir, strcat('validate_local_loss','_', width_name), 'all_losses');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
% save_dir_figures = [save_dir,'\' ,'figures'];
% if ~exist(save_dir_figures,'dir')
%     mkdir(save_dir_figures);
% end

% global struct_name
struct_name= 'L_messung_simplified';
eval(['global ' struct_name]) % measured data will be changed to be globals
% width_names = fieldnames(eval(struct_name));    % width names 
width_struct = strcat(struct_name, '.', width_name); 
global frequency_names
frequency_names = fieldnames(eval(width_struct));       % measured frequencies under each width
number_frequency = length(frequency_names); 
frequency_vector = zeros(number_frequency,1);
num_frequency_used = 0;
%% 
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
for index_frequency = 1 : number_frequency  
        %% hysteresis
        P_hyst_mean_cur_fz = [];
        P_hyst_measurement_fz = [];
        % hysteresis
        P_eddy_mean_cur_fz = [];
        P_eddy_measurement_fz = [];        
        % excess loss
        P_zusatz_mean_cur_fz = [];
        P_zusatz_measurement_fz = [];
        % nonlinear loss
        P_nonlinear_mean_cur_fz = [];
        P_nonlinear_measurement_fz = [];
        % all loss
        P_all_loss_mean_cur_fz = [];
        P_all_loss_meas_measurement_cur_fz = [];
        P_all_loss_meas_original_fz = [];
        % saving J
        Jmax_T_cur_plot = [];
        % power difference between 
        P_difference_local_global_fz = [];
        P_difference_local_original_fz = [];
        % frequency determination
        frequency_name = frequency_names{index_frequency};
        
        %% import of local hysteris loss model
        path_of_alpha = fullfile(path_of_local_para_base_hyst, frequency_name, 'alpha_x', 'alpha_distance.mat');
        path_of_a1 = fullfile(path_of_local_para_base_hyst, frequency_name, 'a1_x', 'a1_distance.mat');
        load(path_of_alpha)
        load(path_of_a1)
        a1_x = a1_distance.a1_x;
        alpha_x = alpha_distance.alpha_x;
        
        %% import of local excess loss model
        path_of_a5 = fullfile(path_of_local_para_base_excess, frequency_name, 'a5_x','a5_distance.mat');
        load(path_of_a5)
        a5_x = a5_distance.a5_x;
        
        %% import of local nonlinear loss model
        path_of_a3 = fullfile(path_of_local_para_base_nonlinear, frequency_name, 'a3_x','a3_distance.mat');
        path_of_a4 = fullfile(path_of_local_para_base_nonlinear, frequency_name, 'a4_x', 'a4_distance.mat');
        load(path_of_a3)
        load(path_of_a4)
        a3_x = a3_distance.a3_x;
        a4_x = a4_distance.a4_x;
        
        %% import of permeability model of undamaged material
        H_no_degration_interp_vector = cut_edge_model.(frequency_name).H_no_degreation_pchip; 
        mu_no_degradation_interp_vector = cut_edge_model.(frequency_name).mu_no_degreation_pchip; 
        
        %% import of delta_mu_c(H)   
        delta_mu_interp_vector = cut_edge_model.(frequency_name).delta_mu_interp_pchip_final  ;
        delta_H_interp_vector = cut_edge_model.(frequency_name).H_interp_pchip_final  ; 
        
        %% eta(x) import, which is frequency dependent
        eta_x = cut_edge_model.(frequency_name).eta.fitting_function;
        parameters = coeffvalues(eta_x);
        depth = 1/ abs(parameters(2));   
        eta_x = @(x) exp(-x/depth);  
        % measured data import "simplified_data.mat", and saved in L_messung
        L_messung = eval(struct_name);
        %% import of the original measured data
        Jmax_T_cur= L_messung.(width_name).(frequency_name).Jmax_T;
        Ps_W_kg = L_messung.(width_name).(frequency_name).Ps_W_kg;
        Ps_W_kg = Ps_W_kg';
        number_of_measurement = length(Jmax_T_cur); 
        %% import of average gloabl a1 and alpha
        a1_messen = A_parameters.(width_name).a1 ; 
        alpha_messen = A_parameters.(width_name).alpha ; 
        %% import of average gloabl a5
        a5_messen = A_parameters.(width_name).a5 ; 
        %% import of average gloabl a5
        a3_messen = A_parameters.(width_name).a3 ; 
        a4_messen = A_parameters.(width_name).a4 ; 
        a2_messen = A_parameters.(width_name).a2 ;         
        %% frequency determination
        frequency_tmp = frequency_vector(index_frequency); 
        %% various of loss functions are defined by using the global average
        % parameters
        HystFunktion = @(J) a1_messen .* J .^(alpha_messen) .* frequency_tmp;
        Excess_Funktion = @(J) a5_messen .* J .^(1.5) .* frequency_tmp.^(1.5) ;
        nonlinear_loss_Funktion = @(J) a2_messen.*a3_messen.* (J.^(2+a4_messen)) .* frequency_tmp.^2 ; 
        eddyFunktion = @(J) a2_messen.*J.^2.*frequency_tmp.^2; 
        %% loop begins
        for j = 1 : number_of_measurement
            J_tmp = Jmax_T_cur(j);
            Ps_W_kg_cur = Ps_W_kg(j);
            if J_tmp >= Polarisation_J_max % validation in low J area. to do
                continue; 
            end
            P_hyst_measurement_cur = HystFunktion(J_tmp); 
            P_zusatz_measurement_cur = Excess_Funktion(J_tmp); 
            P_nonlinear_measurement_cur = nonlinear_loss_Funktion(J_tmp); 
            P_eddy_measurement_cur = eddyFunktion(J_tmp); 
            
            H_tmp = L_messung.(width_name).(frequency_name).Hmax_A_m(j);          
            if (H_tmp > max(H_no_degration_interp_vector) || H_tmp < min(H_no_degration_interp_vector))
                continue; 
            end
            if (H_tmp > max(delta_H_interp_vector) || H_tmp < min(delta_H_interp_vector))
                continue; 
            end
            
            mu_relative_undamaged =  interp1(H_no_degration_interp_vector,mu_no_degradation_interp_vector,H_tmp, 'pchip','extrap'); % mu_r_j at H_j
            delta_mu_c = interp1(delta_H_interp_vector,delta_mu_interp_vector,H_tmp, 'pchip','extrap'); 
            % hysteresis loss  calculatetion by using local model
            P_hyst_func = @(x) a1_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^alpha_x(x).*frequency_tmp;
            % integration to calculate average
            P_hyst_mean_cur = 2*integral(P_hyst_func, 0 , validate_width/2)/validate_width; 
            % update of the loss for calculated and measrued 
            P_hyst_mean_cur_fz = [P_hyst_mean_cur_fz P_hyst_mean_cur];
            P_hyst_measurement_fz = [P_hyst_measurement_fz, P_hyst_measurement_cur];
            
            % the eddy current loss is calculated by the same way as the the gloable model 
            P_eddy_mean_cur = P_eddy_measurement_cur; 
            % update of the loss for calculated and measrued 
            P_eddy_mean_cur_fz = [P_eddy_mean_cur_fz P_eddy_mean_cur]; 
            P_eddy_measurement_fz = [P_eddy_measurement_fz, P_eddy_measurement_cur ];
            
            % excess loss calculatetion by using the local model
            P_zusatz_func = @(x) a5_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^(1.5).* frequency_tmp.^(1.5) ;
            % integration to calculate average
            P_zusatz_mean_cur = 2*integral(P_zusatz_func, 0 , validate_width/2)/validate_width;            
            % update of the loss for calculated and measrued 
            P_zusatz_mean_cur_fz = [P_zusatz_mean_cur_fz P_zusatz_mean_cur];
            P_zusatz_measurement_fz = [P_zusatz_measurement_fz, P_zusatz_measurement_cur];
            
            % nonlinear loss calculatetion by using local model
            P_nonlinear_func = @(x) a2_messen.* a3_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^(2+a4_x(x)).* frequency_tmp.^(2) ;
            % integration to calculate average
            P_nonlinear_mean_cur = 2*integral(P_nonlinear_func, 0 , validate_width/2)/validate_width; 
            % update of the loss for calculated and measrued 
            P_nonlinear_mean_cur_fz = [P_nonlinear_mean_cur_fz P_nonlinear_mean_cur];
            P_nonlinear_measurement_fz = [P_nonlinear_measurement_fz, P_nonlinear_measurement_cur];           
            
            % sum of various losses
            P_all_loss_measurement_cur = P_hyst_measurement_cur + P_zusatz_measurement_cur + P_nonlinear_measurement_cur + P_eddy_measurement_cur;
            P_all_loss_mean_cur = P_hyst_mean_cur + P_zusatz_mean_cur + P_nonlinear_mean_cur + P_eddy_mean_cur;
            % update of the loss for calculated and measrued
            P_all_loss_mean_cur_fz = [P_all_loss_mean_cur_fz P_all_loss_mean_cur];
            P_all_loss_meas_measurement_cur_fz = [P_all_loss_meas_measurement_cur_fz P_all_loss_measurement_cur];
            P_all_loss_meas_original_fz = [P_all_loss_meas_original_fz Ps_W_kg_cur];
            
            % difference between measured values and calculated values in
            % square of the value difference
            P_difference_tmp1 = (P_all_loss_measurement_cur - P_all_loss_mean_cur); % abs or no abs
            P_difference_tmp2 = (Ps_W_kg_cur - P_all_loss_mean_cur);                % abs or no abs
            % update of the variablel P_difference
            P_difference_local_global_fz = [P_difference_local_global_fz P_difference_tmp1]; 
            P_difference_local_original_fz = [P_difference_local_original_fz P_difference_tmp2];
            % update of the Jmax_T_cur_plot for the plots later
            Jmax_T_cur_plot = [Jmax_T_cur_plot J_tmp]; 
        end  
%%
        figure
%         plot(Jmax_T_cur, P_hyst_measurement_fz, Jmax_T_cur, P_hyst_mean_cur_fz); % by Li
%         plot(Jmax_T_cur_plot, P_hyst_measurement_fz, Jmax_T_cur_plot, P_hyst_mean_cur_fz); % by Li
%         bar(Jmax_T_cur_plot , [P_all_loss_meas_measurement_cur_fz ; P_all_loss_mean_cur_fz  ; P_all_loss_meas_original_fz; P_difference_local_original_fz])
        
bar(Jmax_T_cur_plot, P_all_loss_meas_measurement_cur_fz);
hold
bar(Jmax_T_cur_plot, P_all_loss_mean_cur_fz)
bar(Jmax_T_cur_plot, P_all_loss_meas_original_fz)
bar(Jmax_T_cur_plot, P_difference_local_original_fz)
        xlabel('$J_{\rm max}$ in T',  'Interpreter', 'latex', 'Fontsize', 16)
        ylabel('$P_{\rm iron loss}$ in W/kg',  'Interpreter', 'latex', 'Fontsize', 16)
        title(strcat('validation of P_{Fe} ', ' of', '--', width_name,  '--under--',  frequency_name), 'Fontsize', 16 )
%         legend('measured fitted', 'calculated', 'original', 'diff-origi-cal', 'Location', 'Northwest', 'Fontsize', 16)
        
        img = gcf;
        set(gcf,'position',[0.2,0.2,720,405]); 

        save_dir_fz = [save_dir,'\' frequency_name,];
        if ~exist(save_dir_fz,'dir')
            mkdir(save_dir_fz);
        end
        filename = strcat('P_ironloss_compare');
        print(img, '-dpng', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.jpg']);
        print(img, '-depsc', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.eps']);
        print(img, '-dmeta', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.emf']); 
        saveas(img,[save_dir_fz,'\',convertStringsToChars(filename),'.fig']);
        pause(2)
        close 
        %%
        save_matfile_frequency = strcat(save_dir_fz,'\','P_ironloss_compare.mat');
        Phys.measured = P_hyst_measurement_fz; 
        Phys.calculated = P_hyst_mean_cur_fz; 
        Phys.J_measured = Jmax_T_cur_plot; 
        Phys.P_difference = abs(Phys.measured - Phys.calculated);
        
        Pzusatz.measured = P_zusatz_measurement_fz; 
        Pzusatz.calculated = P_zusatz_mean_cur_fz; 
        Pzusatz.J_measured = Jmax_T_cur_plot; 
        Pzusatz.P_difference = abs(Pzusatz.measured -Pzusatz.calculated); 
 
        Pnl.measured = P_nonlinear_measurement_fz; 
        Pnl.calculated = P_nonlinear_mean_cur_fz; 
        Pnl.J_measured = Jmax_T_cur_plot; 
        Pnl.P_difference = abs(Pnl.measured - Pnl.calculated); 
        
        Pfe.global_cal = P_all_loss_meas_measurement_cur_fz; 
        Pfe.local_cal = P_all_loss_mean_cur_fz; 
        Pfe.original = P_all_loss_meas_original_fz; 
        Pfe.diff_local_origi= P_difference_local_original_fz; 
        Pfe.diff_local_global = P_difference_local_global_fz; 
        Pfe.J_measured = Jmax_T_cur_plot;
        
        Peddy.Peddy = P_eddy_mean_cur_fz;
        Peddy.P_difference = zeros(1,length(Jmax_T_cur_plot));
        Peddy.J_measured =  Jmax_T_cur_plot;
        
        P_loss_all.Pfe = Pfe; 
        P_loss_all.Phys = Phys; 
        P_loss_all.Pzusatz = Pzusatz;
        P_loss_all.Pnl = Pnl; 
        P_loss_all.Peddy = Peddy; 
        
        save(save_matfile_frequency,'P_loss_all');

        
        %% update in all frequency
        eval(strcat('allfrequency', '.', frequency_name, '.', 'P_loss_all', '=', 'P_loss_all'))

end
save_matfile_allfrequency = fullfile(save_dir, 'allfrequency.mat');    
save(save_matfile_allfrequency, 'allfrequency');
close all
%% return to the main dierectory
eval(['cd ' main_dir])