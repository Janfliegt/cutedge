% valiateion of the local iron loss parameters determined by the different evoltuon

clear all 
clc;
close all 
%% config parameters
validate_width = 30;
frequency_max = 20;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 这里要改一下
path_of_local_para_base = 'D:\Materials\thesis\schnittkanteeffekt\AW__local_parameter_a1_1_a1_2_alpha1_alpha2\local_loss_parameters_10\hyst\';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\simplified_data.mat');
% path_of_cut_edge_model = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\plot_cut_edge_model_cor_data\','cut_edge_model.mat');
path_of_cut_edge_model = strcat('P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\plot_cut_edge_model_cor_data\','cut_edge_model.mat');
path_of_A_parameters = 'P:\intern_2020_schnittkanten\04-arbeitspakete\00_Messdaten\Testdaten\Material1_Handschlagschere\L\A_parameters\A_parameters.mat';
mu_0 = 4*pi*10^-7 ;
%%
width_name = strcat('B', num2str(validate_width), 'mm');
save_dir = 'validate_result';
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
% save_dir_figures = [save_dir,'\' ,'figures'];
% if ~exist(save_dir_figures,'dir')
%     mkdir(save_dir_figures);
% end
C = who; 
global struct_name
struct_name= C{:};
eval(['global ' struct_name]) % measured data will be changed to be globals
width_names = fieldnames(eval(struct_name));    % width names 

width_struct = strcat(struct_name, '.', width_name); 
global frequency_names
frequency_names = fieldnames(eval(width_struct));       % measured frequencies under each width
number_frequency = length(frequency_names); 
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
frequency_vector = frequency_vector(1:num_frequency_used);
cost = []; 
P_hyst_mean = [];
P_hyst_measurement = [];
load(path_of_cut_edge_model);
load(path_of_A_parameters)
for index_frequency = 1 : number_frequency  
        P_hyst_mean_cur_fz = [];
        P_hyst_measurement_fz = [];
        frequency_name = frequency_names{index_frequency};
        path_of_alpha = strcat(path_of_local_para_base, frequency_name, '\alpha_x\alpha_distance.mat');
        path_of_a1 = strcat(path_of_local_para_base, frequency_name, '\a1_x\a1_distance.mat');
        load(path_of_alpha)
        load(path_of_a1)
        a1_x = a1_distance.a1_x;
        alpha_x = alpha_distance.alpha_x;
        
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
        L_messung = eval(struct_name);
        Jmax_T_cur= L_messung.(width_name).(frequency_name).Jmax_T;
        Ps_W_kg_cur = L_messung.(width_name).(frequency_name).Ps_W_kg;
        number_of_measurement = length(Jmax_T_cur); 

        a1_messen = A_parameters.(width_name).a1 ; 
        alpha_messen = A_parameters.(width_name).alpha ; 
        HystFunktion = @(J) a1_messen .* J .^(alpha_messen) .* frequency_tmp;
        for j = 1 : number_of_measurement
            J_tmp = Jmax_T_cur(j);          
            P_hyst_measurement_cur = HystFunktion(J_tmp); 
            H_tmp = L_messung.(width_name).(frequency_name).Hmax_A_m(j);
            mu_relative_undamaged =  interp1(H_no_degration_interp_vector,mu_no_degradation_interp_vector,H_tmp, 'pchip','extrap'); % mu_r_j at H_j
            delta_mu_c = interp1(delta_H_interp_vector,delta_mu_interp_vector,H_tmp, 'pchip','extrap'); 
            
            P_hyst_func = @(x) a1_x(x).* (mu_0.* H_tmp .* (mu_relative_undamaged - delta_mu_c .* eta_x(x))).^alpha_x(x).*frequency_tmp;
            P_hyst_mean_cur = 2*integral(P_hyst_func, 0 , validate_width/2)/validate_width; 
            P_hyst_mean_cur_fz = [P_hyst_mean_cur_fz P_hyst_mean_cur];
            P_hyst_measurement_fz = [P_hyst_measurement_fz, P_hyst_measurement_cur];
            cost_tmp = abs(P_hyst_measurement_cur - P_hyst_mean_cur).^2;
            cost = [cost cost_tmp]; 
        end  
        figure
        plot(Jmax_T_cur, P_hyst_measurement_fz, Jmax_T_cur, P_hyst_mean_cur_fz);
        xlabel('$J_{max}$',  'Interpreter', 'latex')
        ylabel('$P_{hyst\_mean}$',  'Interpreter', 'latex')
        legend('measured', 'calculated')
        
        img = gcf;
        save_dir_fz = [save_dir,'\' frequency_name,];
        if ~exist(save_dir_fz,'dir')
            mkdir(save_dir_fz);
        end
        filename = strcat('P_hyst_mean');
        print(img, '-dpng', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.jpg']);
        print(img, '-depsc', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.eps']);
        print(img, '-dmeta', '-r600', [save_dir_fz,'\',convertStringsToChars(filename),'.emf']); 
        saveas(img,[save_dir_fz,'\',convertStringsToChars(filename),'.fig']);
        pause(2)
        close 
        save_matfile_frequency = strcat(save_dir_fz,'\','P_hyst_mean.mat');
        save(save_matfile_frequency,'P_hyst_mean_cur_fz');  
        save(save_matfile_frequency,'P_hyst_measurement_fz', '-append');
        P_hyst_mean = [P_hyst_mean; P_hyst_mean_cur_fz];
        P_hyst_measurement = [P_hyst_measurement;P_hyst_measurement_fz];
end

