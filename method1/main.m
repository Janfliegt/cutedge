% The whole process to calculate the cut edge model
clear all;
clc 
close all;
%% setup parameter 
stage = 8;   % from which stage to run the programm, default starting from reading data
original_data_dir = 'D:\Materials\thesis\cut_edge_model\data\L'; % data dir cntains the excel tables of the data
processed_data_dir = 'D:\Materials\thesis\cut_edge_model\method1\data'; % absolute path dir
number_bad_lines = 0;  % some lines at the beginning are bad, they should be deleted
frequency_max = [100]; % 
validate_width = 15;   % don't use data of this width, it is used to validate the result
threshold_F_0 = 0.2;   %in stage 3, abs(F_0-1)<=threshold_F_0, then stops the calculation
hyst_art = 1;          % stage 5, hyst_art = 1 means the formula of Silas, a1.*J.^alpha.*f % hyst_art = 2 means the formula of Silas, a1.*J.^(alpha+beta.*J).*f
Polarisation_J_max = 1.5;
%%
cur_dir = pwd;
addpath(genpath(pwd));
if ~exist(processed_data_dir, 'dir')
    mkdir(processed_data_dir)
end
%% read data 
% output -> raw_data, simplified_data, plot_original_data
if stage == 1 
    data_import_simplification_plot_check(original_data_dir, number_bad_lines, cur_dir,processed_data_dir)
    stage = stage + 1;
end
%% preprocess data
% output -> plot_corrected_data
if stage == 2 
   preprocess_before_schnittkanten_spring_ball_for_all_frequency(cur_dir, processed_data_dir)
   stage = stage + 1;
end
%% optimization to find Fb and some other parameters
% output -> plot_cut_edge_model_cor_data
if stage == 3
    optimierung_main_find_Fb_corrected_data_for_all_frequency(processed_data_dir, frequency_max, threshold_F_0, validate_width, cur_dir)
    stage = stage + 1;
end
%% validate the fitting of cdge edge model
% output -> validate_cut_edge_xxmm
if stage == 4
    validation_cut_edge_model(processed_data_dir, validate_width, cur_dir)
    stage = stage + 1;
end
%% determine A parameters
% output -> A_parameters
if stage == 5
    A_parameter_determination_without_interface(processed_data_dir, hyst_art, cur_dir)
    stage = stage + 1;
end
%% permeability validation
% output -> local_loss_parameters_without_xxmm
if stage == 6
    local_hyst_loss_parameter_modified_without_einzelne_width(processed_data_dir, validate_width, frequency_max, cur_dir);
    local_nonlinear_loss_parameter_modified_without_einzelne_width(processed_data_dir, validate_width, frequency_max, cur_dir)
    local_zusatz_loss_parameter_modified_without_einzelne_width(processed_data_dir, validate_width, frequency_max, cur_dir)
    
    stage = stage + 1;  
end
%% validate local loss hyst
% output -> validate_local_loss_xxmm
if stage == 7
    validatation_of_local_loss_model_all_losses(processed_data_dir, validate_width, frequency_max, Polarisation_J_max, cur_dir);
    
    stage = stage + 1;
end
%% use all widths to calculate the local loss parameters
if stage == 8
    optimierung_main_find_Fb_corrected_data_for_all_frequency(frequency_max, cur_dir, processed_data_dir)
    A_parameter_determination_without_interface(processed_data_dir, hyst_art, cur_dir)
    local_hyst_loss_parameter(processed_data_dir, frequency_max, cur_dir);
    local_nonlinear_loss_parameter(processed_data_dir, frequency_max, cur_dir)
    local_zusatz_loss_parameter(processed_data_dir, frequency_max, cur_dir)  
    stage = stage + 1;
end
%% output model for use in FEM simulation
% a matrix(2-dim) of mu(x, H)
% alpha(x), a1(x), a2(x), a3(x), a4(x), a5(x)



