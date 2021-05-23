% The whole process to calculate the cut edge model
clear all;
clc 
close all;
%% setup parameter
original_data_dir = 'D:\Materials\thesis\cut_edge_model\data\L'; % data dir cntains the excel tables of the data
number_bad_lines = 0; % some lines at the beginning are bad, they should be deleted
used_frequency = [100];
%% Part one: find Fb
cur_dir = pwd;
data_import_simplification_plot_check(original_data_dir, number_bad_lines, cur_dir)
preprocess_before_schnittkanten_spring_ball_for_all_frequency(original_data_dir, cur_dir)
find_Fb(original_data_dir, used_frequency, cur_dir)
