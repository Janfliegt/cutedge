% valiateion of the local iron loss parameters determined by the different evoltuon

clear all 
clc;
close all 
%% config parameters
validate_width =15;  % validation by using the measurement data of validate_width
frequency_max = 200;    % the frequencies smaller than frequency_max will be used for validation

%% polarization_J_max
Polarisation_J_max = 1.5; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% import of local loss model developed without einzelne width
% path_of_local_para_base = 'D:\Materials\thesis\schnittkanteeffekt\AW__local_parameter_a1_1_a1_2_alpha1_alpha2\local_loss_parameters_10\zusatz\';
% path_of_local_para_base = 'P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\local_loss_parameters_without_10mm\zusatz\';
% path_of_local_para_base = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\local_loss_parameters_without_',num2str(validate_width), 'mm\zusatz\');
% path_of_local_para_base = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L_test\local_loss_parameters_without_',num2str(validate_width), 'mm\zusatz\');
path_of_local_para_base = 'local_loss_parameters\zusatz\';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% import of measurement data
% load('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\simplified_data.mat');
load('simplified_data.mat'); 

% path_of_cut_edge_model = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\plot_cut_edge_model_cor_data\','cut_edge_model.mat');
% load('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L_test\simplified_data.mat');

% import of the cut edge model
% path_of_cut_edge_model = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\plot_cut_edge_model_cor_data\','cut_edge_model.mat');
% path_of_cut_edge_model = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L_test\plot_cut_edge_model_cor_data\','cut_edge_model.mat');
path_of_cut_edge_model = 'plot_cut_edge_model_cor_data\cut_edge_model.mat';

% import of averaged A_paramters corresponding to the width
% path_of_A_parameters = 'P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\A_parameters\A_parameters.mat';
% path_of_A_parameters = 'P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\A_parameters_adjusted_nl\A_parameters.mat';
% path_of_A_parameters = 'P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L_test\A_parameters\A_parameters.mat';
 path_of_A_parameters = 'A_parameters\A_parameters.mat';

% mu_0
mu_0 = 4*pi*10^-7 ;
%%
width_name = strcat('B', num2str(validate_width), 'mm');
save_dir = strcat('validate_result','_', width_name, '\zusatz\');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
% save_dir_figures = [save_dir,'\' ,'figures'];
% if ~exist(save_dir_figures,'dir')
%     mkdir(save_dir_figures);
% end
C = who; 
global struct_name
struct_name= C{:};            % todo
eval(['global ' struct_name]) % measured data will be changed to be globals
width_names = fieldnames(eval(struct_name));    % width names 

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
P_zusatz_mean = [];
P_zusatz_measurement = [];

load(path_of_cut_edge_model);
load(path_of_A_parameters)
for index_frequency = 1 : number_frequency  
        P_zusatz_mean_cur_fz = [];
        P_zusatz_measurement_fz = [];
        Jmax_T_cur_plot = []; 
        P_difference = [];
        frequency_name = frequency_names{index_frequency};
                
        path_of_a5 = strcat(path_of_local_para_base, frequency_name, '\a5_x\a5_distance.mat');
        load(path_of_a5)
        a5_x = a5_distance.a5_x;
        
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

        a5_messen = A_parameters.(width_name).a5 ; 
        frequency_tmp = frequency_vector(index_frequency); 
        Excess_Funktion = @(J) a5_messen .* J .^(1.5) .* frequency_tmp.^(1.5) ;

        for j = 1 : number_of_measurement
            J_tmp = Jmax_T_cur(j);          
            if J_tmp >= Polarisation_J_max % validation in low J area. to do
                continue; 
            end
            P_zusatz_measurement_cur = Excess_Funktion(J_tmp); 
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
            P_zusatz_func = @(x) a5_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^(1.5).* frequency_tmp.^(1.5) ;

            % integration to calculate average
            P_zusatz_mean_cur = 2*integral(P_zusatz_func, 0 , validate_width/2)/validate_width; 
            % update of the loss for calculated and measrued 
            P_zusatz_mean_cur_fz = [P_zusatz_mean_cur_fz P_zusatz_mean_cur];
            P_zusatz_measurement_fz = [P_zusatz_measurement_fz, P_zusatz_measurement_cur];
            % difference between measured values and calculated values in
            % square of the value difference
            P_difference_tmp = (P_zusatz_measurement_cur - P_zusatz_mean_cur);
            % update of the variablel P_difference
            P_difference = [P_difference P_difference_tmp]; 
            
            % update of the Jmax_T_cur_plot for the plots later
            Jmax_T_cur_plot = [Jmax_T_cur_plot J_tmp]; 
        end  
        figure
        % plot(Jmax_T_cur, P_zusatz_measurement_fz, Jmax_T_cur, P_zusatz_mean_cur_fz); % by Li
        % plot(Jmax_T_cur_plot, P_zusatz_measurement_fz, Jmax_T_cur_plot, P_zusatz_mean_cur_fz); % by Li
        bar(Jmax_T_cur_plot , [P_zusatz_measurement_fz ; P_zusatz_mean_cur_fz  ; P_difference])
        xlabel('$J_{\rm max}$ in T',  'Interpreter', 'latex', 'Fontsize', 16)
        ylabel('$P_{\rm zusatz-mean}$ in W/kg',  'Interpreter', 'latex', 'Fontsize', 16)
        title(strcat('validation of P_{excess} ', ' of', '--', width_name,  '--under--',  frequency_name), 'Fontsize', 16 )
        legend('measured', 'calculated', 'diff:meas-cal', 'Location', 'Northwest', 'Fontsize', 16)
        
        img = gcf;
        set(gcf,'position',[0.2,0.2,720,405]); 

        save_dir_fz = [save_dir,'\' frequency_name,];
        if ~exist(save_dir_fz,'dir')
            mkdir(save_dir_fz);
        end
        filename = strcat('P_zusatz_compare');
        print(img, '-dpng', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.jpg']);
        print(img, '-depsc', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.eps']);
        print(img, '-dmeta', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.emf']); 
        saveas(img,[save_dir_fz,'\',convertStringsToChars(filename),'.fig']);
        pause(2)
        close 
        save_matfile_frequency = strcat(save_dir_fz,'\','P_zusatz_compare.mat');
        Pzusatz.measured = P_zusatz_measurement_fz; 
        Pzusatz.calculated = P_zusatz_mean_cur_fz; 
        Pzusatz.J_measured = Jmax_T_cur_plot; 
        Pzusatz.P_difference = P_difference; 
        save(save_matfile_frequency,'Pzusatz');
        % save(save_matfile_frequency,'P_zusatz_mean_cur_fz');  
        % save(save_matfile_frequency,'P_zusatz_measurement_fz', '-append');
        
        % why are the following variables needed? saving the data for all
        % frequency?
        % P_zusatz_mean = [P_zusatz_mean; P_zusatz_mean_cur_fz];
        % P_zusatz_measurement = [P_zusatz_measurement;P_zusatz_measurement_fz];
        
        %% update in all frequency
        eval(strcat('allfrequency', '.', frequency_name, '.', 'Pzusatz', '=', 'Pzusatz'))

end
save_matfile_allfrequency =  strcat(save_dir,'\', 'allfrequency','.mat');    
save(save_matfile_allfrequency, 'allfrequency');
close all
