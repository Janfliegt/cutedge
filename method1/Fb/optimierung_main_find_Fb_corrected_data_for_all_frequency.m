function [] = optimierung_main_find_Fb_corrected_data_for_all_frequency(processed_data_dir, frequency_max, threshold_F_0, validate_width, main_dir)
%threshold_F_0: a threshold for F_0, which should equal 1, if abs(F_0-1)<=threshold_F_0, then stops the calculation

%% setups 
search_min = 0.6;
search_max = 1;
%%
load(fullfile(processed_data_dir, 'plot_corrected_data', 'allfrequency.mat'));
%% remove the unused frequency data and unused width 
frequency_fields = fieldnames(allfrequency);
number_of_frequency = length(frequency_fields);
validate_width_str = strcat('B',num2str(validate_width), 'mm');
remove_validate_width = true;  % whether to remove the data of the validate width 
if ~ismember(validate_width_str, width_fields)
    remove_validate_width = false;
end

for i = 1 : number_of_frequency
    frequency_name = frequency_fields{i};
    cur_frequency = str2num(frequency_name(2:end-2));
    if cur_frequency > frequency_max
        allfrequency = rmfield(allfrequency, frequency_name);
        continue
    end
    if remove_validate_width
        tmp = rmfield(allfrequency.(frequency_name), validate_width_str);
        allfrequency.(frequency_name) = tmp;
    end
end
%%
struct_name= 'allfrequency'; 
struct_cut_edge = 'cut_edge_model'; 
save_dir = fullfile(processed_data_dir, 'plot_cut_edge_model_cor_data');
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% frequency name vector extraction
frequency_names = fieldnames(eval(struct_name));
number_of_frequency = length(frequency_names);
%%
% figure number 
figure_num = 1;
global width_names
global width_vector
global f_in_Hz
global f_H_b_matrix
global variance_for_different_H
global delta_mu_c_H
global H_interpolation
global mu_no_degradation
RES = 4;        % restarten mal              % 10-mal Optimierungsversuche
ITE = 12;       % Anzahl der Iterationen
for i_frequency = 1 : number_of_frequency
    frequency_name = frequency_names{i_frequency}; 
    var_name = strcat(struct_name, '.', frequency_name);                       % frequency struct, darunter verschiedene breite
    % width vector for each frequency determine
    %global width_names
    width_names = fieldnames(eval(var_name));                                  % width names under a frequency
    number_width = length(width_names);                                        % number of width, for each width there is a corresponding F(bi)
    %global width_vector
    width_vector = zeros(number_width,1);
    %global f_in_Hz
    f_in_Hz = eval(var_name);
    for j = 1 : number_width
        width_name_txt =  width_names{j} ;
        width_tmp = regexp(width_name_txt, '\d*\d*', 'match');
        if length(width_tmp) == 1
            width_tmp = width_tmp{:};
            width_tmp = str2double(width_tmp);
            width_vector(j) = width_tmp;
        else       
            width_tmp= strcat(width_tmp{1}, '.', width_tmp{end});
            width_tmp = str2double(width_tmp);
            width_vector(j) = width_tmp ;           %  width vector are saved in width_vector
        end      

    end
    width_vector_fit = width_vector(1:end-1);
    % fitting method for F(b) after discrete F(b) are determined
    Fb_fitting_method = 2; % 1 for smoothing spline; 2 for a*exp(x*b); 3 for a*x^b+c
    while 1
      upper = 0.5 * (search_min + search_max);
      [Best_Result, Best_Solution, Fitness_xbest_evolution] = evolutionary(upper, number_width, ITE, RES);
      
    print_DE = 1;
    [z,y] = min(Best_Result); % Auswahl der besten Lösung innherhalb 10-mal Optimierungsversuche
    xx = Best_Solution(y,:);
    %% fit the curve to get Fb(0), which should equal one
    fun = @(params)sseval(params, width_vector_fit', xx);
    params0 = rand(3, 1); % 
    best_params = fminsearch(fun, params0);
    
    yfit = best_params(1) * exp(-best_params(2)* width_vector_fit+ best_params(3));
    % plot the fitted figure and the original data together to check the
    % fitting quality
    figure
    plot(width_vector_fit', xx, '.',width_vector_fit, yfit); 
    xlabel('b in mm', 'Interpreter', 'none', 'Fontsize', 20 );
    ylabel( 'F(b)', 'Interpreter', 'none' );
    
    F_0 = best_params(1) * exp(-best_params(2) * 0 + best_params(3));
    if abs(F_0 - 1) > threshold_F_0
        if F_0 > 1
            search_max = upper;         
        else
            search_min = upper;          
        end
        close all
    else
        break
    end   
      
    end

   %% best results for F(b1)....F(bn-1) bestimmen, und die entsprechende delta_mu_c_H
    
    
    modeling_corrected_data(xx, print_DE);
    figure;
    for n = 1:RES
        x = 1:1:ITE;
        % x = x*P_Size;
        % Label axes
        xlabel( 'number of iterations', 'Interpreter', 'none' );
        ylabel( 'min. mean variance of f(H,b) under different H after each iteration', 'Interpreter', 'none' ); 
        hold on;
        plot(x,Fitness_xbest_evolution(n,:),'blue');   % Verläufe der Zielgröße bei jeder Optimeirungsversuche
    end
    hold off; 
    % saving
    frequency_filename = frequency_name; 
    save_frequency = [save_dir,'\',convertStringsToChars(frequency_filename)];
    if ~exist(save_frequency,'dir')
        mkdir(save_frequency);
    end     
     save_dir_DE_results = [save_frequency,'\' ,'DE_results'];
    if ~exist(save_dir_DE_results,'dir')
        mkdir(save_dir_DE_results);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('convergence_of_DE');
    img =gcf;  
    print(img, '-dpng', '-r600', [save_dir_DE_results,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_dir_DE_results,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_dir_DE_results,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_dir_DE_results,'\',convertStringsToChars(filename),'.fig'])    ;
    % save data
    save_matfile = strcat(save_dir_DE_results,'\','F_b','.mat');
    F_b_values = xx; 
    save(save_matfile,'F_b_values');
    

    %% fitting discrete F(b) to curves 
    Fb_value = [xx, 0]; % 0 for the maximal width case: e.g. 120mm in our case
    [xData, yData] = prepareCurveData( width_vector, Fb_value );

    % Set up fittype and options.
    if Fb_fitting_method == 1
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( 'Method', 'SmoothingSpline' );
        opts.SmoothingParam = 0.121355375301033;
    elseif Fb_fitting_method == 2
        ft = fittype( 'exp1' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [2.844830437318 -0.130834509580792];    
    elseif  Fb_fitting_method == 3
        ft = fittype( 'power2' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [1.6938495127251 -0.4620752608313 0.0424692213487243];
    else
        ft = fittype( '2*a/x*(1-1*exp(-x/(2*b)))', 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = 0.836104367833986;      
    end 
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    figure;
    plot( fitresult); hold on; 
    plot( xData, yData, 'ko','MarkerFaceColor', 'black' ,'MarkerSize', 5)
    hold off;
    xlim([min(width_vector) max(width_vector)]);
    ylim([0 1]);
    txt_legend = [string(strcat('Fb vs width for', '_', frequency_name)), "discrete F(b)" ];
    legend( txt_legend, 'Location', 'NorthEast', 'Interpreter', 'none' );
    xlabel( 'width in mm', 'Interpreter', 'none' );
    ylabel( 'F(b)', 'Interpreter', 'none' );

 % saving plots
    frequency_filename = frequency_name; 
    save_frequency = [save_dir,'\',convertStringsToChars(frequency_filename)];
    if ~exist(save_frequency,'dir')
        mkdir(save_frequency);
    end     
     save_dir_Fb_results = [save_frequency,'\' ,'Fb_results'];
    %save_dir = 'plot_normalized_mu_vs_width';
    if ~exist(save_dir_Fb_results,'dir')
        mkdir(save_dir_Fb_results);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('Fb_curve');
    img =gcf;  
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_dir_Fb_results,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_dir_Fb_results,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_dir_Fb_results,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_dir_Fb_results,'\',convertStringsToChars(filename),'.fig'])    ;
    
    %% save data struct
    F_b.width = width_vector;
    F_b.Fb_value = Fb_value'; 
    F_b.Fb_fitting_function = fitresult;
    F_b.Fb_fitting_goodness = gof;
    % F_b.H = eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'H_final_vector'));
    F_b.H = H_interpolation; 
    % f_H_b_matrix saving for different H und B, page 73, Silas Diss
    %global f_H_b_matrix
    F_b.f_H_b_matrix = f_H_b_matrix ;  % normalized f(H,bi) by the relative mu of 120mm material, only for b ~= 120mm, because that value at b == 120mm is equal to 0
    % variance of F(H,bi) for different b saving, under different H
    %global variance_for_different_H
    F_b.variance_for_different_H = variance_for_different_H;
    % delta_mu_c_H saving
    % global delta_mu_c_H
    F_b.delta_mu_c_H = delta_mu_c_H; 
    % F_b.mu_no_degradation = eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'mu_final_vector')); 
    F_b.mu_no_degradation = mu_no_degradation;  % here mu_no_degradation uses fitted data, not the original 120mm data
    %% plot delta mu_c(H) vs the original mu_undamaged(H)
    % smooth the curves by pchip interpolation
    % 4000 points are to be interpolated based on the qualitativ corrected measurement data
    H_interp_pchip_final = linspace(min(F_b.H), max(F_b.H), 4000)';
    % mu interpolation
    delta_mu_interp_pchip_final = interp1(F_b.H, F_b.delta_mu_c_H, H_interp_pchip_final, 'pchip');     % piesewise interpolation in the 4000 Points 
    % saving data in F_b
    eval(strcat('F_b', '.', 'H_interp_pchip_final', '=', 'H_interp_pchip_final'));   % from now on, the interpolation will use H_interp_pchip 
    eval(strcat('F_b', '.', 'delta_mu_interp_pchip_final', '=', 'delta_mu_interp_pchip_final'));   
    % plot delta mu_c_H and undamaged mu_H in function of H
    figure;
    semilogx(F_b.H_interp_pchip_final, F_b.delta_mu_interp_pchip_final); 
    hold on;
    % semilogx(eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'H_interp_pchip_final')) , eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'mu_interp_pchip_final')) )
    H_no_degreation_pchip = linspace( eval(strcat('min(','f_in_Hz', '.', width_names{end}, '.', 'H_interp_2', ')')) , eval(strcat('max(','f_in_Hz', '.', width_names{end}, '.', 'H_interp_2', ')')), 4000);
    mu_no_degreation_pchip = interp1(eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'H_interp_2')), eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'mu_relative_interp_2')), H_no_degreation_pchip, 'pchip');     % piesewise interpolation in the 4000 Points 
    H_no_degreation_pchip = H_no_degreation_pchip'; 
    mu_no_degreation_pchip = mu_no_degreation_pchip';
    % plot
    semilogx( H_no_degreation_pchip ,  mu_no_degreation_pchip) ; % corrected 120mm data
    xlabel( 'H in A/m', 'Interpreter', 'none', 'Fontsize', 20 );
    ylabel( '$\Delta \mu_c (H)$ und $\mu_{r,u} (H)$',  'Interpreter', 'latex', 'Fontsize', 20 );
    hold off; 
    gcf
    txt_legend = ["mu-decrease", "mu-undamaged"];
    legend(txt_legend, 'Location', 'Northeast', 'Fontsize', 20)  
    
    % saving the relative mu_c together in F_b    
    eval(strcat('F_b', '.', 'H_no_degreation_pchip', '=', 'H_no_degreation_pchip'));
    eval(strcat('F_b', '.', 'mu_no_degreation_pchip', '=', 'mu_no_degreation_pchip' ));   
  % saving plots
    frequency_filename = frequency_name; 
    save_frequency = [save_dir,'\',convertStringsToChars(frequency_filename)];
    if ~exist(save_frequency,'dir')
        mkdir(save_frequency);
    end     
     save_dir_delta_mu_c_results = [save_frequency,'\' ,'delta_mu_c'];
    %save_dir = 'plot_normalized_mu_vs_width';
    if ~exist(save_dir_delta_mu_c_results,'dir')
        mkdir(save_dir_delta_mu_c_results);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('delta_mu_c');
    img =gcf;  
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_dir_delta_mu_c_results,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_dir_delta_mu_c_results,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_dir_delta_mu_c_results,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_dir_delta_mu_c_results,'\',convertStringsToChars(filename),'.fig'])    ;  
    %% determine eta(x) and saving plots
    fb_derivative = differentiate(F_b.Fb_fitting_function, F_b.width ); % F'(b)
  % eta.values = [1 ; fb_derivative.*  F_b.width + F_b.Fb_value; 0]; % eta(b/2) = F'(b)*b+F(b)  und   1, 0 are added to ensure boundary condition of eta(x)
    eta.values = [1 ; fb_derivative.*  F_b.width + F_b.Fb_fitting_function(F_b.width) ; 0]; % eta(b/2) = F'(b)*b+F(b), interpolated value instead of optimized values; for eta < 1
  % eta.values = [ fb_derivative.*  F_b.width + F_b.Fb_fitting_function(F_b.width) ; 0];      % for eta > 1
    eta.distance =[0; F_b.width/2; max(F_b.width)];                  % 0 for distance = 0; max(F_b.width) for endless distance;  for eta < 1
  % eta.distance =[ F_b.width/2; max(F_b.width)];  % for eta > 1
  % data preparing
    [xData, yData] = prepareCurveData( eta.distance, eta.values );

    % Set up fittype and options. exponential function with one parameter:
    % influence depth of the cut edge effect
    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [1 -Inf];
    opts.StartPoint = [1.16993070349615 -0.338661326885645];
    opts.Upper = [1 Inf];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    figure;
    plot( fitresult, xData, yData );
    legend('$\eta (x)$ vs. distance', 'discrete values $\eta(x)$ ', 'Location', 'NorthEast', 'Interpreter', 'latex' );
    % Label axes
    xlabel( 'distance in mm', 'Interpreter', 'none' );
    ylabel( '$\eta (x)$', 'Interpreter', 'latex' );

  % saving plots
    frequency_filename = frequency_name; 
    save_frequency = [save_dir,'\',convertStringsToChars(frequency_filename)];
    if ~exist(save_frequency,'dir')
        mkdir(save_frequency);
    end     
     save_dir_eta_x_results = [save_frequency,'\' ,'eta_x'];
    %save_dir = 'plot_normalized_mu_vs_width';
    if ~exist(save_dir_eta_x_results,'dir')
        mkdir(save_dir_eta_x_results);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('eta_x');
    img =gcf;  
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_dir_eta_x_results,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_dir_eta_x_results,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_dir_eta_x_results,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_dir_eta_x_results,'\',convertStringsToChars(filename),'.fig'])    ;     
    % saving data
    eta.fitting_function = fitresult; 
    eta.goodness = gof; 
    F_b.eta = eta;  
    
    %% relative mu as 2D-lookuptable in function of distance and H
    distance = 0:1:20; % 0 mm to 20 mm with step of 1 mm; can be changed; to do
    H_tmp = F_b.H_no_degreation_pchip ; % H  vector
    number_distance = length(distance); % number of Distance
    number_H = length(H_tmp);           % number of H
    relative_mu_lookup = zeros(number_H, number_distance); % 2d-lookup table for 
    eta_function = F_b.eta.fitting_function; 
    delta_mu_c = F_b.delta_mu_interp_pchip_final;
    mu_undamaged = F_b.mu_no_degreation_pchip;
    for i_distance = 1 : number_distance
        for j_H = 1 : number_H
            mu_tmp = -eta_function(distance(i_distance)).* delta_mu_c(j_H) + mu_undamaged(j_H); 
            relative_mu_lookup(j_H, i_distance) = mu_tmp;
        end
    end
    % plot mu in function of H under difference distance
    figure; 
    hold on;
    for i_distance = 1 : number_distance
        plot(H_tmp, relative_mu_lookup(:,i_distance))
    end    
    hold off;
    xlabel( 'H in A/m', 'Interpreter', 'none', 'Fontsize', 20 );
    ylabel( '$\mu_{r,u} (H)$',  'Interpreter', 'latex', 'Fontsize', 20 );
    txt_legend = string(strcat('D',  string(distance), 'mm')) ;
    legend( txt_legend, 'Location', 'NorthEast', 'Interpreter', 'none' );   
    % save plots
    frequency_filename = frequency_name; 
    save_frequency = [save_dir,'\',convertStringsToChars(frequency_filename)];
    if ~exist(save_frequency,'dir')
        mkdir(save_frequency);
    end     
     save_dir_mu_von_distance_results = [save_frequency,'\' ,'mu_distance_H'];
    %
    if ~exist(save_dir_mu_von_distance_results,'dir')
        mkdir(save_dir_mu_von_distance_results);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('mu_distance_H');
    img =gcf;  
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_dir_eta_x_results,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_dir_eta_x_results,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_dir_eta_x_results,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_dir_mu_von_distance_results,'\',convertStringsToChars(filename),'.fig'])    ;      
    % save data
    mu_distance_H.distance = distance'; 
    mu_distance_H.H = H_tmp; 
    mu_distance_H.mu = relative_mu_lookup; 
    F_b.mu_distance_H = mu_distance_H;     
    %% save daten under each frequency
    save_matfile_frequency = strcat(save_frequency,'\',frequency_name,'.mat');
    eval(strcat(frequency_name,'_cut_edge',  '=',  'F_b'))   
    save(save_matfile_frequency,strcat(frequency_name, '_cut_edge'));
    
    % save daten under the struct with all frequency
    eval(strcat(struct_cut_edge, '.', frequency_name, '=', strcat(frequency_name,'_cut_edge')))
    close all   
end
save_matfile_allfrequency =  fullfile(save_dir, strcat(struct_cut_edge,'.mat'));    
save(save_matfile_allfrequency, struct_cut_edge);
close all
%% different delta mu_c in one figure
figure;
% b=bone(number_of_frequency); % color depth from black  to white in number_width+2 stufen
% load('plot_cut_edge_model_power\cut_edge_model.mat')
hold on; 
xlabel( 'H in A/m', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$\Delta \mu_{c}(H) $', 'Interpreter', 'latex', 'fontsize',12 ); 
hold on
for j = 1: number_of_frequency
    f_Hz_tmp = eval(strcat(struct_cut_edge, '.', frequency_names{j})); 
    H_tmp = f_Hz_tmp.H_interp_pchip_final; 
    delta_mu_c = f_Hz_tmp.delta_mu_interp_pchip_final;
    plot(H_tmp,delta_mu_c )
end
hold off; 
txt_legend = string(frequency_names);
legend(txt_legend, 'Location', 'Northeast')
title('$\Delta \mu_{c}(H)$ under different frequency', 'Interpreter', 'latex', 'fontsize',12)
% saving figure
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('delta_mu_c_all_frequency');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_dir,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_dir,'\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_dir,'\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_dir,'\',convertStringsToChars(filename),'.fig'])    ;  

%% different eta(x) in one figure
figure;
b=bone(number_of_frequency+2); % color depth from black  to white in number_width+2 stufen
%load('plot_cut_edge_model\cut_edge_model.mat')
hold on; 
xlabel( 'distance in mm', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$\eta (x) $', 'Interpreter', 'latex', 'fontsize',12 ); 
xlim([0 120])
ylim([0 1])
% hold on
for j = 1: number_of_frequency
    f_Hz_tmp = eval(strcat(struct_cut_edge, '.', frequency_names{j})); 
    distance = [0:0.5:120]; 
    eta = f_Hz_tmp.eta.fitting_function(distance);
    plot(distance, eta, 'Color' ,b(j,:))
end
hold off; 
txt_legend = string(frequency_names);
legend(txt_legend, 'Location', 'Northeast')
title('$\eta (x)$ under different frequency', 'Interpreter', 'latex', 'fontsize',12)
% saving figure
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('eta_all_frequency');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_dir,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_dir,'\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_dir,'\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_dir,'\',convertStringsToChars(filename),'.fig']); 
close all
%% return to the main directory
eval(['cd ' main_dir])

