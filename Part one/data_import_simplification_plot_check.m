function [] = data_import_simplification_plot_check(original_data_dir, number_bad_lines, cur_dir)

%% daten importieren in matlab als struct
path_actual = original_data_dir;                                             % actual directory name
[match, noMatch] = regexp(path_actual,'\','match','split');                  % split through '\'
filename_actual= noMatch{end};                                              % the last one will be used to name the variable of measurement raw data

foldername = dir(original_data_dir);                                                      % subfolders under the actual directory will be listed
number_foldername = length(foldername);                                     % number of subfolders

for i = 1:number_foldername  
    name_folder = foldername(i).name; 
    if isempty(regexp(name_folder, '\d', 'once')) || ~endsWith(name_folder, 'mm')                  % if no number in the name of the folder, continue
        continue
    end
    path_name = strcat(foldername(i).folder, '\', foldername(i).name);     % fullpath for files under subfolder
    eval(['cd ' path_name])                                                 % enter subfolder where xlsx files are saved
    
    filename = dir('*.xlsx');
    number_files = length(filename);
    for j = 1 : number_files
        tmp = readtable(strcat(filename(j).folder,'\', filename(j).name) , 'Sheet','Results'); 
        tmp(1:number_bad_lines,:) = []; 
        % name adjusted to discard '.' '-' '+' symbols
        name_file = filename(j).name;                                       % name adjusted to replace '.' '-' '+' symbols with '_'
        filename_new = strrep(name_file,'-','_'); 
        filename_new = strrep(filename_new,'.','_'); 
        filename_new = strrep(filename_new,',','_');
        filename_new = extractBefore(filename_new,'_xlsx');                   
        v = genvarname(filename_new);         
        %
        eval([v '=tmp'])
        % name of folder postprocess so that no bad text symbols like '.'
        % and ',' exist
        name_folder = strrep(name_folder, '.','_');
        name_folder = strrep(name_folder, ',','_');
        %
        eval([strcat(filename_actual,'_messung_raw.','B',name_folder, '.') v '=' v])
        eval(['clear ' v])
    end
    
end
eval(['cd ' path_actual])
save('raw_data.mat', strcat(filename_actual,'_messung_raw'))

%% data extracted from original raw_data to reduce the matrix size, that will be behandled further
% vereinfachte Datensammlung mit wichtigen Daten speichert B, H. 
width_names = fieldnames(eval(strcat(filename_actual,'_messung_raw'))); 
number_of_width = length(width_names);
width_vector = zeros(number_of_width,1); 
for i =1 : number_of_width
    field_name_txt =  width_names{i};
    width_tmp = regexp(field_name_txt, '\d*\d*', 'match');
    if length(width_tmp) == 1
        width_tmp = width_tmp{:};
        width_tmp = str2double(width_tmp);
        width_vector(i) = width_tmp;
    else       
        width_tmp= strcat(width_tmp{1}, '.',width_tmp{end});
        width_tmp = str2double(width_tmp);
        width_vector(i) = width_tmp ;       %  width vector are saved in width_vector
    end   
end
width_vector = sort(width_vector);
width_vector_text = string(width_vector);
width_vector_text = strrep(width_vector_text,'.','_');

frequency_cell_array = cell(number_of_width,1);                             % frequencies for each width, 
for i = 1   :   number_of_width  
    frequency_subfield_names = fieldnames(eval(strcat(filename_actual,'_messung_raw','.', width_names{i} ))); % frequency cell array under each width 
    number_of_frequency = length(frequency_subfield_names); 
    frequency_vector = zeros(number_of_frequency,1);
    for j = 1: number_of_frequency
        j;
        frequency_tmp = eval(strcat(filename_actual,'_messung_raw','.', width_names{i}, '.', frequency_subfield_names{j}, '{:,7}')); % Frequency_Hz_ or Frequenz_Hz_ todo
        frequency_vector(j) = round(mean(frequency_tmp));       
    end
    frequency_vector = sort(frequency_vector);
    frequency_cell_array{i} = frequency_vector; 
end    
frequency_vector_text = string(frequency_vector);


for i = 1 : number_of_width
    width_name_text = strcat('B', width_vector_text(i), 'mm');
    frequency_names_under_width = fieldnames(eval(strcat(filename_actual,'_messung_raw','.', width_name_text )));
    for j = 1 : number_of_frequency
        frequency_name = frequency_vector_text(j);
        frequency_name_text = strcat(frequency_name,'Hz');
        for k = 1 : number_of_frequency
            if  contains(frequency_names_under_width{k,:}, frequency_name_text)
                table_tmp = eval(strcat(filename_actual,'_messung_raw','.', width_name_text, '.',frequency_names_under_width{k,:}  ));                
                break;
            else
                continue;
            end
        end
        try
            table_simplified = table_tmp(:,{'Sample', 'Jmax_T_', 'Hmax_A_m_', 'Bmax_T_', 'x_r'}); % todo sample oder probe
        catch
            table_simplified = table_tmp(:,{'Probe', 'Jmax_T_', 'Hmax_A_m_', 'Bmax_T_', 'x_r'}); % todo sample oder probe    
        end
        table_simplified.Properties.VariableNames = {'Sample', 'Jmax_T', 'Hmax_A_m', 'Bmax_T', 'mu_relative'}; % to do sample oder probe
        eval(strcat(filename_actual,'_messung_simplified','.', width_name_text, '.', 'f', frequency_name, 'Hz',  '=', 'table_simplified'  ) )
    end

end

save('simplified_data.mat', strcat(filename_actual,'_messung_simplified')) % sorted simplified data save


%% plot relative mu vs Feldstärke für verschieden Breite jeweils für verschiedene Frequenz
% data saving: mu = f(H) for different width under each frequency
% plot saving: mu = f(H) for different width under each frequency
save_dir = 'plot_original_data';
mkdir(save_dir)
for i = 1: number_of_frequency
    figure_num = i;
    figure(figure_num);
    hold on
    xlabel('H [A/m]', 'fontsize', 15)
    ylabel('\mu relative', 'fontsize', 15)
    title(strcat('f= ', frequency_vector_text(i), 'Hz'))
    for j = 1 :  number_of_width
        simplified_data_table_tmp = eval(strcat(filename_actual,'_messung_simplified', '.' , strcat('B', width_vector_text(j), 'mm'), strcat('.', 'f', frequency_vector_text(i), 'Hz')));
        H = simplified_data_table_tmp.Hmax_A_m; 
        mu_relative = simplified_data_table_tmp.mu_relative; 
        B = simplified_data_table_tmp.Bmax_T; 
        plot(H, mu_relative,'*-')
    % struct: B, H, mu for each witdh jeweils for different frequency        
        eval(strcat('f', frequency_vector_text(i), 'Hz', '.', strcat('B', width_vector_text(j), 'mm'), '.H', '= H' ) ); 
        eval(strcat('f', frequency_vector_text(i), 'Hz', '.', strcat('B', width_vector_text(j), 'mm'), '.B', '= B' ) );
        eval(strcat('f', frequency_vector_text(i), 'Hz', '.', strcat('B', width_vector_text(j), 'mm'), '.mu_relative', '= mu_relative' ) );
    end
    txt_legend = strcat('B', width_vector_text, 'mm'); 
    txt_legend = replace(txt_legend,'_','.');
    legend(txt_legend)
    hold off
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('f', frequency_vector_text(i), 'Hz');
%     curve_reasonable = input('Are theses curves of the current frequency reasonable?, input 1 for yes, 0 for no: ');
%     if curve_reasonable ~= 1
%         continue
%     end
%     curve_reasonable=1;
    save_path = [save_dir,'\',convertStringsToChars(filename)];
    if ~exist(save_path,'dir')
        mkdir(save_path);
    end
    img = gcf; 
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_path,'\',convertStringsToChars(filename),'.jpg'])
    print(img, '-depsc', '-r600', [save_path,'\',convertStringsToChars(filename),'.eps'])
    print(img, '-dmeta', '-r600', [save_path,'\',convertStringsToChars(filename),'.emf'])
    saveas(img,[save_path,'\',convertStringsToChars(filename),'.fig'])
    % save the data
    save_matfile_frequency = strcat(save_path,'\',strcat('f', frequency_vector_text(i), 'Hz'),'.mat');
    save(save_matfile_frequency,strcat('f', frequency_vector_text(i), 'Hz'));
    
    eval(strcat('allfrequency.',strcat('f', frequency_vector_text(i), 'Hz'), '=', strcat('f', frequency_vector_text(i), 'Hz')))
end
save_matfile_allfrequency =  strcat(save_dir,'\', 'allfrequency','.mat');    
save(save_matfile_allfrequency, 'allfrequency');
close all
%% check matrix 
save_checkmatrix_dir = 'check_measure_data';
mkdir(save_checkmatrix_dir);
frequency_field = fieldnames(allfrequency);
number_of_frequency = length(frequency_field);

check_vector = cell( number_of_frequency,2); 
check_vector(:,1) = frequency_field; 
check_vector(:,2) = num2cell(zeros(number_of_frequency,1)); 
for i = 1 : number_of_frequency
    frequency_name = frequency_field{i};
    width_field = fieldnames(eval(strcat('allfrequency.', frequency_name )));
    number_of_width = length(width_field);
    %eval(strcat(frequency_name,'_','checkmatrix', '=', 'zeros(number_of_width)'))
    mu_max_width = zeros(number_of_width,1);  % save maximum mu for each width
    for j = 1 : number_of_width       
        width_name = width_field{j}; 
        mu_vector = eval(strcat('allfrequency.', frequency_name, '.', width_name, '.', 'mu_relative' ));
        mu_max_width(j) = max(mu_vector);
    end
    [~,index_mu_max_width] = max(mu_max_width); % max mu within different width measurement   
    width_name_max_mu = width_field{index_mu_max_width} ; % width_name_max_mu represents the width name with max mu
    [~, index_mu_max] = max( eval(strcat('allfrequency.', frequency_name, '.', width_name_max_mu, '.', 'mu_relative' ))); % index_mu_max represent the index to position the max mu
    H_cor_mu_max = eval(strcat('allfrequency.', frequency_name, '.', width_name_max_mu, '.', 'H', strcat('(', int2str(index_mu_max), ')') )); % H_cor_mu_max to save the max H
    H_cor_mu_min = min(eval(strcat('allfrequency.', frequency_name, '.', width_name_max_mu, '.', 'H'))); 
    % between H_cor_mu_min and H_cor_mu_max, the mu relative will be interpolated for different width
    H_interpolate_vector = linspace(H_cor_mu_min, H_cor_mu_max, 20); 
    % for each width, the mu for different H_interpolate_vector will be interpolated
    figure;    
    xlabel('H [A/m]', 'fontsize', 15)
    ylabel('\mu relative interpolated', 'fontsize', 15)
    title(frequency_field{i}) 
    hold on; 
    for j = 1 : number_of_width       
        width_name = width_field{j}; 
        mu_vector = eval(strcat('allfrequency.', frequency_name, '.', width_name, '.', 'mu_relative' ));
        H_vector = eval(strcat('allfrequency.', frequency_name, '.', width_name, '.', 'H' ));
        mu_interpolate_vector = interp1(H_vector,mu_vector,H_interpolate_vector); 
        eval(strcat('BH_interpolate.', frequency_name, '.', width_name, '.', 'H_interpolate', '=', 'H_interpolate_vector' )  )
        eval(strcat('BH_interpolate.', frequency_name, '.', width_name, '.', 'mu_interpolate', '=', 'mu_interpolate_vector' )  )
        % plot
        plot(H_interpolate_vector, mu_interpolate_vector)
    end 
    hold off
    set(gcf,'position',[0.2,0.2,720,405]);
    txt_legend = string( width_field); 
    txt_legend = replace(txt_legend,'_','.');
    legend(txt_legend,'location', 'northwest');
    
    checkmatrix = zeros(number_of_width); 
    for k = 1 : number_of_width
        width_name_k = width_field{k};
        for j = k +1 : number_of_width
            width_name_j = width_field{j};
            logic_between_widths = eval( strcat('BH_interpolate.', frequency_name, '.', width_name_k,'.','mu_interpolate', '>', 'BH_interpolate.', frequency_name, '.', width_name_j,'.','mu_interpolate'));
            checkmatrix(k,j) =  nnz(logic_between_widths) > 0; 
        end       
    end
    [row_index, columb_index] = find(checkmatrix == 1);
    bad_width_cell_array = cell(0);
    if ~isempty(row_index)
         
        disp( strcat('for', ' ' , frequency_name, ' the messung for the following widths muss be renewed:'))
        for g = 1: length(row_index)
        disp(width_field{row_index(g)}) 
        disp(width_field{columb_index(g)})        
        bad_width_cell_array = [bad_width_cell_array, width_field{row_index(g)}, width_field{columb_index(g)}]; 
        disp(';')
        end
        
        check_vector{i,2} = 1; % check_vector adjustation, 1 for bad measurement for the corresponding frequency, 0 for good measurement 
    end
    eval( strcat('allcheckmatrix.', frequency_name, '.', 'checkmatrix', '=checkmatrix;' ))
    % table is besser to show checkmatrix
    checkmatrix_table = array2table(checkmatrix);
    checkmatrix_table.Properties.RowNames = width_field(:)'; 
    checkmatrix_table.Properties.VariableNames = width_field(:)';
    %
    save_check_path = [save_checkmatrix_dir,'\',convertStringsToChars(frequency_name)];
    if ~exist(save_check_path,'dir')
        mkdir(save_check_path);
    end  
    img =gcf; 
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_check_path,'\',convertStringsToChars(frequency_name),'.jpg'])
    print(img, '-depsc', '-r600', [save_check_path,'\',convertStringsToChars(frequency_name),'.eps'])
    print(img, '-dmeta', '-r600', [save_check_path,'\',convertStringsToChars(frequency_name),'.emf'])
    saveas(img,[save_check_path,'\',convertStringsToChars(frequency_name),'.fig']) 
    % save check matrix
    save_matfile_checkmatrix = strcat(save_check_path,'\', convertStringsToChars(frequency_name),'.mat');
    if ~isempty(bad_width_cell_array) 
        save(save_matfile_checkmatrix,'checkmatrix', 'bad_width_cell_array', 'checkmatrix_table');    
    end
end
save_matfile_all_checkmatrix = strcat(save_checkmatrix_dir,'\','all_checkmatrix' ,'.mat');
save(save_matfile_all_checkmatrix, 'allcheckmatrix'); 

save_matfile_checkvector = strcat(save_checkmatrix_dir,'\','checkvector' ,'.mat');
save(save_matfile_checkvector, 'check_vector'); 

save_matfile_BH_interpolate = strcat(save_checkmatrix_dir,'\','BH_interpolate' ,'.mat');
save(save_matfile_BH_interpolate, 'BH_interpolate'); 

eval(['cd ' cur_dir]) 
clear all;
close all;

