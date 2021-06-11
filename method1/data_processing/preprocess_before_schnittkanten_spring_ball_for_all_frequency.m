function [] = preprocess_before_schnittkanten_spring_ball_for_all_frequency(cur_dir, processed_data_dir)

% preprocessing of measurement data
allfrequency_data = fullfile(processed_data_dir, 'plot_original_data', 'allfrequency.mat');
load(allfrequency_data)
C = who; 
struct_name= C{:}; 
save_dir = fullfile(processed_data_dir, 'plot_corrected_data');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end
% frequency name vector extraction
frequency_names = fieldnames(eval(struct_name));
number_of_frequency = length(frequency_names); 
% figure number 
figure_num = 1;
for i = 1 : number_of_frequency
% width vector for each frequency determine
    frequency_name = frequency_names{i}; 
    var_name = strcat(struct_name, '.', frequency_name);                       % frequency struct, darunter verschiedene breite
    width_names = fieldnames(eval(var_name));                                  % width names under a frequency
    number_width = length(width_names);                                        % number of width, for each width there is a corresponding F(bi)
    width_vector = zeros(number_width,1);
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
            width_vector(j) = width_tmp ;                                      %  width vector are saved in width_vector
        end      

    end

    % spring-ball method to preprocess data
    spring_factor = 0.05; 
    for j = 1: number_width
        width_name_einzeln = width_names{j};
        if j== 1
            H_interp = eval(strcat(var_name,'.',width_name_einzeln,'.','H'));
            mu_relative_interp   =eval( strcat(var_name,'.',width_name_einzeln,'.','mu_relative'));
            eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp', '=H_interp')); 
            eval(strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp', '=mu_relative_interp')); 
        else
            H_origin = eval(strcat(var_name,'.',width_name_einzeln,'.','H'));
            mu_relative_origin=eval( strcat(var_name,'.',width_name_einzeln,'.','mu_relative'));
            H_interp = H_origin;
            mu_relative_interp = zeros(length(H_origin),1); 

            % corrected data from the width before
            H_origin_vector_one_width_later = eval(strcat(var_name,'.',width_names{j-1},'.','H_interp')); % important H_interp instead of H
            mu_relative_vector_one_width_later = eval(strcat(var_name,'.',width_names{j-1},'.','mu_relative_interp')); % important H_interp instead of H
            %
            number_of_H = length(H_origin);
            for k = 1: number_of_H
                mu_relative_interp_tmp = interp1(H_origin_vector_one_width_later,mu_relative_vector_one_width_later,H_interp(k),'linear' ,'extrap'); 
                if mu_relative_interp_tmp > mu_relative_origin(k)
                    mu_relative_interp(k) = mu_relative_interp_tmp + (mu_relative_interp_tmp - mu_relative_origin(k) ) * spring_factor; 
                else
                    mu_relative_interp(k) = mu_relative_origin(k);
                end
            end
            eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp', '=H_interp')); 
            eval(strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp', '=mu_relative_interp'));         
        end

    end

    % second time spring from higher to lower 
    spring_factor_second = 0.05; 
    for j =  number_width : -1:  1
        width_name_einzeln = width_names{j};
        if j== number_width
            H_interp = eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp'));
            mu_relative_interp   =eval( strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp'));
            eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp_2', '=H_interp')); 
            eval(strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp_2', '=mu_relative_interp')); 
        else
            H_origin = eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp'));
            mu_relative_origin=eval( strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp'));
            H_interp = H_origin;
            mu_relative_interp = zeros(length(H_origin),1); 

            % corrected data from the width before
            H_origin_vector_one_width_later = eval(strcat(var_name,'.',width_names{j+1},'.','H_interp')); % important H_interp instead of H
            mu_relative_vector_one_width_later = eval(strcat(var_name,'.',width_names{j+1},'.','mu_relative_interp')); % important H_interp instead of H
            %
            number_of_H = length(H_origin);
            for k = 1: number_of_H
                mu_relative_interp_tmp = interp1(H_origin_vector_one_width_later,mu_relative_vector_one_width_later,H_interp(k),'linear'); 
                if mu_relative_interp_tmp < mu_relative_origin(k)
                    mu_relative_interp(k) = mu_relative_interp_tmp - ( mu_relative_origin(k)-mu_relative_interp_tmp ) * spring_factor_second; 
                else
                    mu_relative_interp(k) = mu_relative_origin(k);
                end
            end
            eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp_2', '=H_interp')); 
            eval(strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp_2', '=mu_relative_interp'));         
        end

    end



    % plot with curves from interpolated data and original data

    figure(figure_num);
    hold on
    xlabel('H [A/m]', 'fontsize', 15)
    ylabel('\mu relative', 'fontsize', 15)
    title(strcat(frequency_name, ' first correction'))
    % color depth define for each witdh
    color_matrix = zeros(number_width,3); 
    rgb_min = 0.5;
    rgb_max = 0.9;
    for k = 1: number_width
        color_matrix(k,:) = [rgb_min rgb_min rgb_min] + (k-1) * 1/number_width * (rgb_max-rgb_min);
    end

    txt_legend = [];
    for j = 1 :  number_width
        % original data
        width_name_einzeln = width_names{j};
        H_vector_origin = eval(strcat(var_name,'.',width_name_einzeln,'.','H')); 
        mu_relative_origin = eval( strcat(var_name,'.',width_name_einzeln,'.','mu_relative')); 
        plot(H_vector_origin, mu_relative_origin,'.' ,'MarkerSize',12 ,'Color', color_matrix(j,:), 'linewidth', 1.3)
        txt_legend = [txt_legend, string( strcat(width_name_einzeln,'_ori'))]; 
        % interpolated data  
        H_vector_interp = eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp')); 
        mu_relative_interp = eval(strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp'));
        plot(H_vector_interp, mu_relative_interp,'.--' ,'MarkerSize',5 , 'linewidth', 1.3)
        txt_legend = [txt_legend, string(strcat(width_name_einzeln,'_interp'))]; 
    end

    txt_legend = replace(txt_legend,'_','.');
    legend(txt_legend)
    hold off
    set(gcf,'position',[0.2,0.2,720,405]); 
    % save figure
    filename = frequency_name;
    save_path = [save_dir,'\',convertStringsToChars(filename)];
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
     img =gcf; 
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_path,'\',convertStringsToChars(filename),'_first','.jpg'])
    print(img, '-depsc', '-r600', [save_path,'\',convertStringsToChars(filename),'_first','.eps'])
    print(img, '-dmeta', '-r600', [save_path,'\',convertStringsToChars(filename),'_first','.emf'])
    saveas(img,[save_path,'\',convertStringsToChars(filename),'_first','.fig'])    
    figure_num = figure_num +1; 
    % figure after second correction
    figure(figure_num);
    hold on
    xlabel('H [A/m]', 'fontsize', 15)
    ylabel('\mu relative', 'fontsize', 15)
    title(strcat(frequency_name, ' second correction'))
    % color depth define for each witdh
    color_matrix = zeros(number_width,3); 
    rgb_min = 0.5;
    rgb_max = 0.9;
    for k = 1: number_width
        color_matrix(k,:) = [rgb_min rgb_min rgb_min] + (k-1) * 1/number_width * (rgb_max-rgb_min);
    end
    % plot
    txt_legend = [];
    for j = 1 :  number_width
        % original data
        width_name_einzeln = width_names{j};
        H_vector_origin = eval(strcat(var_name,'.',width_name_einzeln,'.','H')); 
        mu_relative_origin = eval( strcat(var_name,'.',width_name_einzeln,'.','mu_relative')); 
        plot(H_vector_origin, mu_relative_origin,'.' ,'MarkerSize',12 ,'Color', color_matrix(j,:), 'linewidth', 1.3)
        txt_legend = [txt_legend, string( strcat(width_name_einzeln,'_ori'))]; 
        % interpolated data  
        H_vector_interp = eval(strcat(var_name,'.',width_name_einzeln,'.','H_interp_2')); 
        mu_relative_interp = eval(strcat(var_name,'.',width_name_einzeln,'.','mu_relative_interp_2'));
        plot(H_vector_interp, mu_relative_interp,'.--' ,'MarkerSize',5 , 'linewidth', 1.3)
        txt_legend = [txt_legend, string(strcat(width_name_einzeln,'_interp'))]; 
    end

    txt_legend = replace(txt_legend,'_','.');
    legend(txt_legend)
    hold off
    set(gcf,'position',[0.2,0.2,720,405]); 
    % save figure
    filename = frequency_name;
    save_path = [save_dir,'\',convertStringsToChars(filename)];
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
     img =gcf; 
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_path,'\',convertStringsToChars(filename) ,'.jpg'])
    print(img, '-depsc', '-r600', [save_path,'\',convertStringsToChars(filename),'.eps'])
    print(img, '-dmeta', '-r600', [save_path,'\',convertStringsToChars(filename),'.emf'])
    saveas(img,[save_path,'\',convertStringsToChars(filename),'.fig'])    
    figure_num = figure_num +1; 
    % save daten under each frequency
    save_matfile_frequency = strcat(save_path,'\',frequency_name,'.mat');
    eval(strcat(frequency_name,'_corrected',  '=',  'eval', '(' ,'var_name', ')'))
    
    save(save_matfile_frequency,strcat(frequency_name, '_corrected'));    
end

save_matfile_allfrequency =  strcat(save_dir,'\', struct_name,'.mat');    
save(save_matfile_allfrequency, struct_name);
eval(['cd ' cur_dir]) 
close all;