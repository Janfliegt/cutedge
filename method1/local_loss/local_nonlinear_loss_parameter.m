% different evoltuon is used to determine local iron loss parameters
function [] = local_nonlinear_loss_parameter(processed_data_dir, frequency_max, main_dir)
 
%% import measured original data related to iron loss
load(fullfile(processed_data_dir, 'simplified_data.mat'))
global struct_name
global width_names
global width_vector
struct_name = 'L_messung_simplified';
width_names = fieldnames(eval(struct_name)); % width names 

eval(['global ' struct_name]) % measured data will be changed to be globals
number_of_width = length(width_names); 
width_vector = zeros(number_of_width, 1);
save_dir =fullfile(processed_data_dir, strcat('local_loss_parameters'));
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end

%% determination of width vector
for j = 1 : number_of_width
    width_name_txt =  width_names{j} ;
    width_tmp = regexp(width_name_txt, '\d*\d*', 'match');
    if length(width_tmp) == 1
        width_tmp = width_tmp{:};
        width_tmp = str2double(width_tmp);   
        width_vector(j) = width_tmp;
    else       
        width_tmp= strcat(width_tmp{1}, '.', width_tmp{end});
        width_tmp = str2double(width_tmp);
        width_vector(j) = width_tmp;
    end      

end

%% determination of frequency vector
width_name = width_names{1}; 
width_struct = strcat(struct_name, '.', width_name); 
global frequency_names
frequency_names = fieldnames(eval(width_struct));       % measured frequencies under each width
number_frequency = length (frequency_names); 

global frequency_vector
frequency_vector = zeros(number_frequency,1);
%% determination of the frequency vector, it is assumed that for each width
% the same frequency vector is measured
num_frequency_used = 0;
for j = 1 : number_frequency
    frequency_name_txt =  frequency_names{j};
    frequency_tmp = regexp(frequency_name_txt, '\d*\d*', 'match');
    if length(frequency_tmp) == 1
        frequency_tmp = frequency_tmp{:};
        frequency_tmp = str2double(frequency_tmp);
        if frequency_tmp <= frequency_max
            frequency_vector(j) = frequency_tmp;
            num_frequency_used = num_frequency_used+1;
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
%% import of A-parameters for different width
global A_parameters     % measured data will be changed to be globals

path_of_A_parameters = fullfile(processed_data_dir, 'A_parameters', 'A_parameters.mat');
load(path_of_A_parameters)

% import cut edge modell
global cut_edge_model

path_of_cut_edge_model = fullfile(processed_data_dir, 'plot_cut_edge_model_cor_data','cut_edge_model.mat');
load(path_of_cut_edge_model);
%% initialization of frequency 
% frequency range and polarisation range to determine a_5
global f_min_a3_a4; 
global f_max_a3_a4; 
global J_min_a3_a4; 
global J_max_a3_a4; 
%% basis depth to determine the profile of local loss parameters, the depth is determined based on cut edge model, it do not need to be precise
global depth_basis_local_loss_parameter

for index_frequency = 1 : number_frequency
    f_min_a3_a4 = frequency_vector(index_frequency); % minimal frequency to extract a_3 a_4 parameter
    f_max_a3_a4 = frequency_vector(index_frequency); % maximal frequency to extract a_3 a_4 parameter
    J_min_a3_a4 = 0; % minimal polarisation to extract a_3 a_4 parameter
    J_max_a3_a4 = inf; % maximal frequency to extract a_3 a_4 parameter    
    % depth in the normalization
    depth_basis_local_loss_parameter = 5;   
    
    if f_min_a3_a4 > frequency_max
       break 
    end
    
    %% evolutionary algorithm for local hysteresis parameters a1,1 a1,2 alpha1, alpha2
    number_designvariable = 4;             % a3_1 a3_2, a4_1 a4_2
    range = [[2 4 2 4];[0 0 -2 0]];         % a3 ~ [1 1+2], a4~[0,2]
    % initiation
    P_Size = number_designvariable*4;         % Population Size
    VAR = number_designvariable;              % Number of design Variables
    ITE = 16;              % Anzahl der Iterationen
    Population = zeros(P_Size,VAR);        % Individuen in der Population
    Fitness_Population = zeros(P_Size,3);  % erste Spalte für Fitnessbewertung (feasible und unfeasible beide Fälle).
                                           % zweite spalte für feasible oder unfeasible,
                                           % dritte Spalte  für Verletzungsgrad der Nebenbedingungen
    grenze = range;             % untere und obere Grenze für jede Dimension
    % erste Zeile von grenze für Obergrenze, zweite Zeile für Untergrenze
    CR = 0.9;                  % Rekombinationsrate
    F_alpha_0 =   0.6;           % eine Skalierungsfaktor: groß für bessere lokale Suche, klein für bessere Exploration
    F_alpha_end = 0.3;
    F_beta = 0.1;              % Skalierungsfaktor
    S_r0 = 0.55;               %  Sr beschreibt  Wahrscheinlichkeit, wie viel der Vergleich aufgrund der Zielgröße anstatt der  Nebenbedingung stattfindet
    S_rGmax = 0.025;           % Endwert für Sr
    delta_Sr = 2*(S_r0 - S_rGmax ) / ITE;     % Abnahme der Sr nach einer Iteration
    child = zeros(1,VAR);                  % Speichert die Designvariablen des neusten Kinds
    Fitness_child = zeros(1,3);            % erste Spalte für Fitnessbewertung (feasible und unfeasible beide Fälle).
                                           % zweite Spalte gültig nur wenn unfeasible, um die Violation grad zu zeigen,
                                           % dritte spalte für feasible oder unfeasible 
    u_G_1 = zeros(1,VAR);                  % Speicherung der Designvariablen des bis dahin besten Kindes
    Fitness_u_G_1 = zeros(1,3);            % wie  der Vektor Fitness des bis dahin besten Kindes
    print_index= 0;
    n0 = 1;                                % 2 Kinder für ein Elter
    x_best = zeros(1,VAR);                 % Speicherung der Designvariablen der besten Lösung in der Population in x_best 
    Fitness_x_best = zeros(1);             % Speicherung der Eigenschaft der besten Lösung, inklusive Zielgröße, Verletzunggrad der Nebenbedingung, feasibility
    RES = 5; % restarten mal              % 10-mal Optimierungsversuche
    Fitness_xbest_evolution = zeros(RES,ITE); % Speicherung die Verläufe der Bewertung der besten Lösung in jedem Optimierungsversuch

    Best_Result = zeros(RES,1);            % Speicherung der letztendlichen Fitnesserwerte der besten Lösung nach jedem Optimierungsversuch
    Best_Solution = zeros(RES,VAR);        % Speicherung der letztendlichen Designvariabllen der besten Lösung nach jedem Optimierungsversuch



    % RES restarting of evolutionary algorithm
    for ggg = 1:RES
        % initiation der Population innerhalb der Grenze
        for j = 1: P_Size
            j
                   while 1
                       for nn = 1 : VAR
                           Population(j,nn) = (range(1,nn)-range(2,nn))*rand + range(2,nn);
                       end

                        [Temp1,Temp2,Temp3]= model_nonlinear_loss_HP(Population(j, :), print_index);  
                        if Temp1 ~= -1                        % gleich -1, bedeutet Recycling, bis ungelich -1
                           Fitness_Population(j,:) = [Temp1,Temp2,Temp3]; 
                            break;
                        end
                   end
        end

        % bestimmen x_best innerhalb der anfangs Population  aufgrund der
        % Constraintsbehandlung: feasible besser als unfeasible; wenn beide
        % feasible, dann Individuum mit niedrigeren Fitness-Funktion besser; 
        % wenn beide unfeasible, dann Indivium mit weniger Verletzungsgrad besser
        term = 1;  % term speichert die Numerierung der besten Lösung
        for jj = 2:P_Size
            if Fitness_Population(term,2) == 0 &&  Fitness_Population(jj,2) == 0  % beide gültig
               if  Fitness_Population(term,1) > Fitness_Population(jj,1)
                   term = jj;
                   continue;
               end
            elseif Fitness_Population(term,2) == 1 &&  Fitness_Population(jj,2) == 0 % eine gültig, die andere ungültig
                   term = jj;
                   continue;
            elseif Fitness_Population(term,2) == 1 &&  Fitness_Population(jj,2) == 1 % beide ungültig
                   if Fitness_Population(term,3) >  Fitness_Population(jj,3) 
                       term = jj;
                       continue;
                   end
            end
        end

        % Iteration beginnt
        for ii = 1: ITE
            ii
            F_alpha = F_alpha_0+ (F_alpha_end - F_alpha_0)* (ii/ITE)^4;  % Anpassung des Skalierungsfaktor mit der Iteration
            for i = 1 : P_Size
                   for k = 1: n0              
                       while 1
                           j_rand = unidrnd(VAR);           % ganze Zufallszahl, damit mindenstens eine Dimension sich verändert
                           while 1
                              r1 =  unidrnd(P_Size);        % r1, r2, r3 sind ganze Zufallszahl, zur Mutation
                              r2 =  unidrnd(P_Size);
                              r3 =  unidrnd(P_Size);
                              if r1 ~= r2 && r1 ~= r3 && r2 ~= r3 && r1~= i && r2~= i && r3~= i  % unterschiedlich und ungleich i  
                                 break; 
                              end
                           end

                           % binäre Rekombination
                           for j = 1:VAR
                               if rand < CR || j == j_rand      
                                        child(j) = Population(r3,j) + F_alpha* (Population(term,j)-Population(r2,j)) + F_beta * (Population(i,j) - Population(r1,j) );
                               else
                                        child(j) = Population(i,j); 
                               end
                           end

                           % hier beginnt boundary einschänkung
                           for mm = 1:VAR

                               if child (mm) >  grenze(1,mm)
                                  child(mm) = 2* grenze(1,mm) - child (mm);
                               elseif child (mm) <  grenze(2,mm)
                                  child(mm) = 2* grenze(2,mm) - child (mm); 
                               end

                           end
                           % hier boundary einschränkung fertig    
    %                        diff_child = diff(child); % todo
    %                        if nnz(diff_child > 0) > 0
    %                            continue;
    %                        end


                            [Temp1,Temp2,Temp3]= model_nonlinear_loss_HP(child, print_index);  
                            if Temp1 ~= -1 
                               Fitness_child = [Temp1,Temp2,Temp3];           
                                break;
                            end
                        end

                       if k > 1
                           % Vergleich der Kinder basierend auf Nebenbedingung
                                if Fitness_child(1,2) == 0 && Fitness_u_G_1(1,2) == 0  % beide gültig
                                    if Fitness_child(1,1) < Fitness_u_G_1(1,1)
                                        u_G_1 = child; 
                                        Fitness_u_G_1 = Fitness_child;
                                    end                            
                                elseif Fitness_child(1,2) == 0 &&  Fitness_u_G_1(1,2) == 1 % eine gültig, die ander ungültig
                                        u_G_1 = child; 
                                        Fitness_u_G_1 = Fitness_child;
                                elseif Fitness_child(1,2) == 1 &&  Fitness_u_G_1(1,2) == 1 % beide ungültig
                                    if  Fitness_child(1,3) < Fitness_u_G_1(1,3)  % spalte 2 steht für , wie viel Violation gegen die Nebenbedingung es ist
                                        u_G_1 = child; 
                                        Fitness_u_G_1 = Fitness_child;
                                    end                            
                                end


                       else 
                                u_G_1 = child;
                                Fitness_u_G_1 = Fitness_child;
                       end
                   end

                 if ii == 1
                    Sr =  S_r0 - delta_Sr;  % Sr steht für die Möglichkeit, wie möglich der Verlgeich nur aufgrund der Fitnessbewertung statt der Feasiblität stattfindet.
                 else 
                     Sr = Sr - delta_Sr;    % die Wahrscheinlichkeit Sr verkleinert sich mit der Iterationen
                 end
                            if rand < Sr    % Vergleich nur aufgrund der Zeilfunktion
                                     if Fitness_Population(i,1) >=   Fitness_u_G_1(1,1)   % bei bessere Zielgröße, die Ersetzung findet statt
                                        Population(i,:) = u_G_1;
                                        Fitness_Population(i,:) = Fitness_u_G_1;
                                     else                      
                                        Population(i,:) =  Population(i,:);               
                                     end
                            else
                                     % Vergleich  aufgrund der Nebenbedingung
                                     if Fitness_Population(i,2) == 0 && Fitness_u_G_1(1,2) == 0
                                            if Fitness_Population(i,1) > Fitness_u_G_1(1,1)
                                                    Population(i,:) = u_G_1;
                                                    Fitness_Population(i,:) = Fitness_u_G_1;                                       
                                            end
                                     elseif Fitness_Population(i,2) == 1 && Fitness_u_G_1(1,2) == 0
                                                    Population(i,:) = u_G_1;
                                                    Fitness_Population(i,:) = Fitness_u_G_1; 
                                     elseif  Fitness_Population(i,2) == 1 &&  Fitness_u_G_1(1,2) == 1 
                                            if Fitness_Population(i,3) > Fitness_u_G_1(1,3)
                                                    Population(i,:) = u_G_1;
                                                    Fitness_Population(i,:) = Fitness_u_G_1;
                                            end
                                     end
                            end
            end
            % erneuen x_best bzw. Population_best, Vergleich 
            term = 1; % term speicher die Numerierung der Besten in der Popolation
        for jj = 2:P_Size
            if Fitness_Population(term,2) == 0 &&  Fitness_Population(jj,2) == 0
               if  Fitness_Population(term,1) > Fitness_Population(jj,1)
                   term = jj;
                   continue;
               end
            elseif Fitness_Population(term,2) == 1 &&  Fitness_Population(jj,2) == 0
                   term = jj;
                   continue;
            elseif Fitness_Population(term,2) == 1 &&  Fitness_Population(jj,2) == 1
                   if Fitness_Population(term,3) >  Fitness_Population(jj,3) 
                       term = jj;
                       continue;
                   end
            end
        end
           Fitness_xbest_evolution(ggg,ii) = Fitness_Population(term,1) ;    % 
        end
           Best_Result(ggg) = Fitness_Population(term,1);   % Speicherung der Zielgröße der besten Lösung
           Best_Solution(ggg,:) = Population(term,:);       % Speicherung der Designvariablen der besten Lösung
    end
    %% convergence bahavior
    print_index = 1;
    [z,y] = min(Best_Result);               % Auswahl der besten Lösung innherhalb 10-mal Optimierungsversuche
    %x = Population(term,:);
    figure;
    hold on;
    % Label axes
    xlabel( 'number of restart of algorithm', 'Interpreter', 'none' );
    ylabel( 'min. fitness function at end of each restart', 'Interpreter', 'none' ); 
    %
    plot(Best_Result,'blue');              % Darstellung der Zielgrößen aller Optimierungsversuche
    hold off;
    xx = Best_Solution(y,:);
    model_nonlinear_loss_HP(xx, print_index);
    figure;
    xlabel( 'number of iterations', 'Interpreter', 'none' );
    ylabel( 'cost function value after each iteration', 'Interpreter', 'none' ); 
    title(strcat('convergence behavior of DE to determine local parameters for ',frequency_names{index_frequency}))   
    hold on;
    for n = 1:RES
        x = 1:1:ITE;     
        plot(x,Fitness_xbest_evolution(n,:),'blue');   % Verläufe der Zielgröße bei jeder Optimeirungsversuche
    end
    hold off
    frequency_name = frequency_names{index_frequency};
    save_frequency = [save_dir,'\','nonlinear','\',convertStringsToChars(frequency_name)];
    if ~exist(save_frequency,'dir')
        mkdir(save_frequency);
    end
     
    save_dir_DE_results = [save_frequency,'\' ,'DE_results'];
    if ~exist(save_dir_DE_results,'dir')
        mkdir(save_dir_DE_results);
    end
    set(gcf,'position',[0.2,0.2,720,405]);
    img = gcf;
    filename = strcat('convergence_of_DE');
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_dir_DE_results,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_dir_DE_results,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_dir_DE_results,'\',convertStringsToChars(filename),'.emf']); 
    saveas(img,[save_dir_DE_results,'\',convertStringsToChars(filename),'.fig']) 
    
    %% xx to a3_1, a3_2, a4_1, a4_2
    a3_1 = xx(1);               % 
    a3_2 = xx(2);               %    
    a4_1 = xx(3);               % 
    a4_2 = xx(4);               % 
    
    %% definition of local loss parameter functions  
    % a3(x)    
    a3_u = A_parameters.(width_names{end}).a3   ;
    a3_x = @ (x) a3_u.*(1+a3_1.*exp(-a3_2.*x/depth_basis_local_loss_parameter));
    % data strutc for later saving
    a3_distance.a3_u = a3_u; 
    a3_distance.a3_1 = a3_1;  
    a3_distance.delta_a3 = depth_basis_local_loss_parameter/ a3_2; 
    % a1_distance.a1_x = a1_x;
    delta_a3 = depth_basis_local_loss_parameter/ a3_2;  
    a3_distance.a3_x = @ (x) a3_u.*(1+a3_1.*exp(-x/delta_a3)); 
    
    % eval(strcat(frequency_name,'_local_a1',  '=', 'a1_distance'))   
    eval(strcat('allfrequency', '.', frequency_name, '.', 'a3_distance', '=', 'a3_distance'))
    % plot a3_x 
    figure;
    hold on; 
    xlabel( 'distance in mm', 'Interpreter', 'none','Fontsize', 16  );
    ylabel( '$a_3(x)$', 'Interpreter', 'latex', 'Fontsize', 16  );
    title(strcat('a_3(x) ', ' for frequency ', frequency_name), 'Fontsize', 16 )
    distance = [0:0.1:20]; 
    plot(distance, a3_x(distance))
    hold off; 
    % saving a3_x
    save_a3_x = [save_frequency,'\' ,'a3_x'];

    if ~exist(save_a3_x,'dir')
        mkdir(save_a3_x);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('a3_x');
    img =gcf;  
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_a3_x,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_a3_x,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_a3_x,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_a3_x,'\',convertStringsToChars(filename),'.fig'])    ;      
    % save a3_x data
    data_name = 'a3_distance' ; 
    save_matfile_frequency = strcat(save_a3_x,'\',data_name,'.mat');
    save(save_matfile_frequency,'a3_distance');     
    %% 
    % plot a4_x
    a4_u = A_parameters.(width_names{end}).a4   ;
    a4_x = @ (x) a4_u.*(1-a4_1.*exp(-a4_2.*x/depth_basis_local_loss_parameter));      
    % data strutc for later saving
    a4_distance.a4_u = a4_u; 
    a4_distance.a4_1 = a4_1; 
    a4_distance.delta_a4 = depth_basis_local_loss_parameter/ a4_2;
    % alpha_distance.alpha_x = alpha_x; 
    delta_a4 = depth_basis_local_loss_parameter/ a4_2;  
    a4_distance.a4_x = @ (x) a4_u .* (1-a4_1.* exp(-x/delta_a4));       
    eval(strcat('allfrequency', '.', frequency_name, '.', 'a4_distance', '=', 'a4_distance'))
    figure;
    distance = [0:0.1:20]; 
    plot(distance, a4_x(distance))
    xlabel( 'distance in mm', 'Interpreter', 'none','Fontsize', 20  );
    ylabel( '$a_4(x)$', 'Interpreter', 'latex', 'Fontsize', 20  );
    title(strcat('a4(x)', ' for frequency ', frequency_name), 'Fontsize', 20 )
    ylim([0 a4_u*2])
    % saving alpha_x
    save_a4_x = [save_frequency,'\' ,'a4_x'];
    if ~exist(save_a4_x,'dir')
        mkdir(save_a4_x);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('alpha_x');
    img =gcf;  
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_a4_x,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_a4_x,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_a4_x,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_a4_x,'\',convertStringsToChars(filename),'.fig'])    ;  

    % save alpha_distance data
    data_name = 'a4_distance' ; 
    save_matfile_frequency = strcat(save_a4_x,'\',data_name,'.mat');
    save(save_matfile_frequency,'a4_distance'); 
end
save_matfile_allfrequency =  fullfile(save_dir,'nonlinear', 'allfrequency.mat');    
save(save_matfile_allfrequency, 'allfrequency');
close all
%% different a3(x) in one figure
save_allfrequency = fullfile(save_dir,'nonlinear');
save_a3_x_frequency = fullfile(save_allfrequency, 'a3_x_frequency'); 
if ~exist(save_a3_x_frequency,'dir')
    mkdir(save_a3_x_frequency);
end

figure;
hold on; 
xlabel( 'distance in mm', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$a_3(x)$', 'Interpreter', 'latex', 'fontsize',12 );
depth_vs_frequency = []; 
frequency_vs_depth = [];
for j = 1: number_frequency
    frequency_name = frequency_names{j};

    if frequency_vector(j) > frequency_max
       break 
    end
    
    a3_x = allfrequency.(frequency_name).a3_distance.a3_x   ;  
    distance = [0:0.1:20];
    plot(distance, a3_x(distance));
    depth_vs_frequency = [depth_vs_frequency allfrequency.(frequency_name).a3_distance.delta_a3 ];
    frequency_vs_depth = [frequency_vs_depth frequency_vector(j)]; 
    
end
hold off; 
txt_legend = string(frequency_names(frequency_vector<=frequency_max));
legend(txt_legend, 'Location', 'Northeast')
title('$a_3(x)$ under different frequency', 'Interpreter', 'latex', 'fontsize',12)
% saving figure
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('a3_distance_all_frequency');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_a3_x_frequency ,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_a3_x_frequency, '\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_a3_x_frequency , '\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_a3_x_frequency,'\',convertStringsToChars(filename),'.fig'])    ; 
close all; 

figure;
hold on
xlabel( 'frequency in Hz', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$\delta_{a3}$ in mm', 'Interpreter', 'latex', 'fontsize',12 );
plot(frequency_vs_depth, depth_vs_frequency, 'o--r', 'Linewidth', 2)
title('Einflusstiefe von a3 in dependency of frequency')
hold off
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('a3_Einflusstiefe');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_a3_x_frequency ,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_a3_x_frequency, '\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_a3_x_frequency , '\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_a3_x_frequency,'\',convertStringsToChars(filename),'.fig'])    ; 
close all; 
%% different a4(x) in one figure
save_allfrequency = fullfile(save_dir,'nonlinear');
save_a4_x_frequency = fullfile(save_allfrequency, 'a4_x_frequency'); 
if ~exist(save_a4_x_frequency,'dir')
    mkdir(save_a4_x_frequency);
end

figure;
hold on; 
xlabel( 'distance in mm', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$a_4(x)$', 'Interpreter', 'latex', 'fontsize',12 );
depth_vs_frequency = []; 
frequency_vs_depth = [];
for j = 1: number_frequency
    frequency_name = frequency_names{j};
    
    if frequency_vector(j) > frequency_max
       break 
    end
    
    a4_x = allfrequency.(frequency_name).a4_distance.a4_x   ;  
    distance = [0:0.1:20];
    plot(distance, a4_x(distance));
    depth_vs_frequency = [depth_vs_frequency allfrequency.(frequency_name).a4_distance.delta_a4 ]; 
    frequency_vs_depth = [frequency_vs_depth frequency_vector(j)]; 
end
hold off; 
txt_legend = string(frequency_names(frequency_vector<=frequency_max));
legend(txt_legend, 'Location', 'Northeast')
title('$a_4(x)$ under different frequency', 'Interpreter', 'latex', 'fontsize',12)
% saving figure
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('a4_distance_all_frequency');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_a4_x_frequency ,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_a4_x_frequency, '\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_a4_x_frequency , '\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_a4_x_frequency,'\',convertStringsToChars(filename),'.fig'])    ; 
close all; 

figure;
hold on
xlabel( 'frequency in Hz', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$\delta_{a4}$ in mm', 'Interpreter', 'latex', 'fontsize',12 );
plot(frequency_vs_depth, depth_vs_frequency, 'o--r', 'Linewidth', 2)
title('Einflusstiefe von a4 in dependency of frequency')

hold off
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('a4_Einflusstiefe');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_a4_x_frequency ,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_a4_x_frequency, '\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_a4_x_frequency , '\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_a4_x_frequency,'\',convertStringsToChars(filename),'.fig'])    ; 
close all; 
%% 
eval(['cd ' main_dir])