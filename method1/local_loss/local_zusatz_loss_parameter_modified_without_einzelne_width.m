% different evoltuon is used to determine local iron loss parameters
function [] = local_zusatz_loss_parameter_modified_without_einzelne_width(processed_data_dir, validate_width, frequency_max, main_dir)
 
%% import measured original data related to iron loss
load(fullfile(processed_data_dir, 'simplified_data.mat'))
global struct_name
global width_names
global width_vector
struct_name = 'L_messung_simplified';
width_names = fieldnames(eval(struct_name)); % width names 

eval(['global ' struct_name]) % measured data will be changed to be globals
number_of_width = length(width_names); 
width_vector = zeros(number_of_width-1,1);
save_dir =fullfile(processed_data_dir, strcat('local_loss_parameters','_without_', num2str(validate_width),'mm'));
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end
%% determination of width vector
index_width = 1;
for j = 1 : number_of_width
    width_name_txt =  width_names{j} ;
    width_tmp = regexp(width_name_txt, '\d*\d*', 'match');
    if length(width_tmp) == 1
        width_tmp = width_tmp{:};
        width_tmp = str2double(width_tmp);
        if width_tmp ~= validate_width     
            width_vector(index_width) = width_tmp;
            index_width = 1 + index_width;
        end
    else       
        width_tmp= strcat(width_tmp{1}, '.', width_tmp{end});
        width_tmp = str2double(width_tmp);
        if width_tmp ~= validate_width     
            width_vector(index_width) = width_tmp;
            index_width = 1 + index_width;
        end
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
% determination of the frequency vector, it is assumed that for each width
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
%% 
global f_min_a5; 
global f_max_a5; 
global J_min_a5; 
global J_max_a5;    
global depth_basis_local_loss_parameter

for index_frequency = 1 : number_frequency
    % initialization of frequency 
    % frequency range and polarisation range to determine a_5

    f_min_a5 = frequency_vector(index_frequency); % minimal frequency to extract a_5 parameter
    f_max_a5 = frequency_vector(index_frequency); % maximal frequency to extract a_5 parameter
    J_min_a5 = 0;   % minimal Polarisation to extract a_5 parameter
    J_max_a5 = inf; % maximal Polarisation to extract a_5 parameter

    if f_min_a5 > frequency_max
       break 
    end

    %% basis depth to determine the profile of local loss parameters, the depth is determined based on cut edge model, it do not need to be precise
    depth_basis_local_loss_parameter = 5; 

    %% evolutionary algorithm for local zusatzloss parameters a5_1 a5_2
    number_designvariable = 2;                 % a5_1 a5,2
    range = [[8 10 ];[0 0.1]];         % a5_1 ~ [0 8], a5_2~[0.1 , 10]
    % initiation
    P_Size = number_designvariable*4;         % Population Size
    VAR = number_designvariable;              % Number of design Variables
    ITE = 16;              % Anzahl der Iterationen
    Population = zeros(P_Size,VAR);        % Individuen in der Population
    Fitness_Population = zeros(P_Size,3);  % erste Spalte f??r Fitnessbewertung (feasible und unfeasible beide F??lle).
                                           % zweite spalte f??r feasible oder unfeasible,
                                           % dritte Spalte  f??r Verletzungsgrad der Nebenbedingungen
    grenze = range;             % untere und obere Grenze f??r jede Dimension
    % erste Zeile von grenze f??r Obergrenze, zweite Zeile f??r Untergrenze
    CR = 0.9;                  % Rekombinationsrate
    F_alpha_0 =   0.6;           % eine Skalierungsfaktor: gro?? f??r bessere lokale Suche, klein f??r bessere Exploration
    F_alpha_end = 0.3;
    F_beta = 0.1;              % Skalierungsfaktor
    S_r0 = 0.55;               %  Sr beschreibt  Wahrscheinlichkeit, wie viel der Vergleich aufgrund der Zielgr????e anstatt der  Nebenbedingung stattfindet
    S_rGmax = 0.025;           % Endwert f??r Sr
    delta_Sr = 2*(S_r0 - S_rGmax ) / ITE;     % Abnahme der Sr nach einer Iteration
    child = zeros(1,VAR);                  % Speichert die Designvariablen des neusten Kinds
    Fitness_child = zeros(1,3);            % erste Spalte f??r Fitnessbewertung (feasible und unfeasible beide F??lle).
                                           % zweite Spalte g??ltig nur wenn unfeasible, um die Violation grad zu zeigen,
                                           % dritte spalte f??r feasible oder unfeasible 
    u_G_1 = zeros(1,VAR);                  % Speicherung der Designvariablen des bis dahin besten Kindes
    Fitness_u_G_1 = zeros(1,3);            % wie  der Vektor Fitness des bis dahin besten Kindes
    print_index= 0;
    n0 = 2;                                % 2 Kinder f??r ein Elter
    x_best = zeros(1,VAR);                 % Speicherung der Designvariablen der besten L??sung in der Population in x_best 
    Fitness_x_best = zeros(1);             % Speicherung der Eigenschaft der besten L??sung, inklusive Zielgr????e, Verletzunggrad der Nebenbedingung, feasibility
    RES = 5; % restarten mal              % 10-mal Optimierungsversuche
    Fitness_xbest_evolution = zeros(RES,ITE); % Speicherung die Verl??ufe der Bewertung der besten L??sung in jedem Optimierungsversuch

    Best_Result = zeros(RES,1);            % Speicherung der letztendlichen Fitnesserwerte der besten L??sung nach jedem Optimierungsversuch
    Best_Solution = zeros(RES,VAR);        % Speicherung der letztendlichen Designvariabllen der besten L??sung nach jedem Optimierungsversuch



    % RES restarting of evolutionary algorithm
    for ggg = 1:RES
        % initiation der Population innerhalb der Grenze
        for j = 1: P_Size
            j
                   while 1
                       for nn = 1 : VAR
                           Population(j,nn) = (range(1,nn)-range(2,nn))*rand + range(2,nn);
                       end

                        [Temp1,Temp2,Temp3]= model_zusatz_loss_HP(Population(j, :), print_index);  
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
        term = 1;  % term speichert die Numerierung der besten L??sung
        for jj = 2:P_Size
            if Fitness_Population(term,2) == 0 &&  Fitness_Population(jj,2) == 0  % beide g??ltig
               if  Fitness_Population(term,1) > Fitness_Population(jj,1)
                   term = jj;
                   continue;
               end
            elseif Fitness_Population(term,2) == 1 &&  Fitness_Population(jj,2) == 0 % eine g??ltig, die andere ung??ltig
                   term = jj;
                   continue;
            elseif Fitness_Population(term,2) == 1 &&  Fitness_Population(jj,2) == 1 % beide ung??ltig
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
                           j_rand = unidrnd(VAR);           % ganze Zufallszahl, damit mindenstens eine Dimension sich ver??ndert
                           while 1
                              r1 =  unidrnd(P_Size);        % r1, r2, r3 sind ganze Zufallszahl, zur Mutation
                              r2 =  unidrnd(P_Size);
                              r3 =  unidrnd(P_Size);
                              if r1 ~= r2 && r1 ~= r3 && r2 ~= r3 && r1~= i && r2~= i && r3~= i  % unterschiedlich und ungleich i  
                                 break; 
                              end
                           end

                           % bin??re Rekombination
                           for j = 1:VAR
                               if rand < CR || j == j_rand      
                                        child(j) = Population(r3,j) + F_alpha* (Population(term,j)-Population(r2,j)) + F_beta * (Population(i,j) - Population(r1,j) );
                               else
                                        child(j) = Population(i,j); 
                               end
                           end

                           % hier beginnt boundary einsch??nkung
                           for mm = 1:VAR

                               if child (mm) >  grenze(1,mm)
                                  child(mm) = 2* grenze(1,mm) - child (mm);
                               elseif child (mm) <  grenze(2,mm)
                                  child(mm) = 2* grenze(2,mm) - child (mm); 
                               end

                           end
                           % hier boundary einschr??nkung fertig    
    %                        diff_child = diff(child); % todo
    %                        if nnz(diff_child > 0) > 0
    %                            continue;
    %                        end


                            [Temp1,Temp2,Temp3]= model_zusatz_loss_HP(child, print_index);  
                            if Temp1 ~= -1 
                               Fitness_child = [Temp1,Temp2,Temp3];           
                                break;
                            end
                        end

                       if k > 1
                           % Vergleich der Kinder basierend auf Nebenbedingung
                                if Fitness_child(1,2) == 0 && Fitness_u_G_1(1,2) == 0  % beide g??ltig
                                    if Fitness_child(1,1) < Fitness_u_G_1(1,1)
                                        u_G_1 = child; 
                                        Fitness_u_G_1 = Fitness_child;
                                    end                            
                                elseif Fitness_child(1,2) == 0 &&  Fitness_u_G_1(1,2) == 1 % eine g??ltig, die ander ung??ltig
                                        u_G_1 = child; 
                                        Fitness_u_G_1 = Fitness_child;
                                elseif Fitness_child(1,2) == 1 &&  Fitness_u_G_1(1,2) == 1 % beide ung??ltig
                                    if  Fitness_child(1,3) < Fitness_u_G_1(1,3)  % spalte 2 steht f??r , wie viel Violation gegen die Nebenbedingung es ist
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
                    Sr =  S_r0 - delta_Sr;  % Sr steht f??r die M??glichkeit, wie m??glich der Verlgeich nur aufgrund der Fitnessbewertung statt der Feasiblit??t stattfindet.
                 else 
                     Sr = Sr - delta_Sr;    % die Wahrscheinlichkeit Sr verkleinert sich mit der Iterationen
                 end
                            if rand < Sr    % Vergleich nur aufgrund der Zeilfunktion
                                     if Fitness_Population(i,1) >=   Fitness_u_G_1(1,1)   % bei bessere Zielgr????e, die Ersetzung findet statt
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
           Best_Result(ggg) = Fitness_Population(term,1);   % Speicherung der Zielgr????e der besten L??sung
           Best_Solution(ggg,:) = Population(term,:);       % Speicherung der Designvariablen der besten L??sung
    end
   %% best results for F(b1)....F(bn) bestimmen, und die entsprechende delta_mu_c_H
    print_index = 1;
    [z,y] = min(Best_Result);               % Auswahl der besten L??sung innherhalb 10-mal Optimierungsversuche
    %x = Population(term,:);
    figure;
    hold on;
    % Label axes
    xlabel( 'number of restart of algorithm', 'Interpreter', 'none' );
    ylabel( 'min. fitness function at end of each restart', 'Interpreter', 'none' ); 
    %
    plot(Best_Result,'blue');              % Darstellung der Zielgr????en aller Optimierungsversuche
    hold off;
    xx = Best_Solution(y,:);
    model_zusatz_loss_HP(xx, print_index);
    figure;
    for n = 1:RES
        x = 1:1:ITE;
        % x = x*P_Size;
        % Label axes
        xlabel( 'number of iterations', 'Interpreter', 'none' );
        ylabel( 'cost function value after each iteration', 'Interpreter', 'none' ); 
        hold on;
        plot(x,Fitness_xbest_evolution(n,:),'blue');   % Verl??ufe der Zielgr????e bei jeder Optimeirungsversuche
    end
    hold off
    frequency_name = frequency_names{index_frequency};
    save_frequency = fullfile(save_dir,'zusatz',convertStringsToChars(frequency_name));
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
    saveas(img,[save_dir_DE_results,'\',convertStringsToChars(filename),'.fig'])    ;
    
     %% xx to a5_1, a5_2
    a5_1 = xx(1);               % 
    a5_2 = xx(2);               % 
    
    %% definition of local loss parameter functions  
    % a5(x) 
    a5_u = A_parameters.(width_names{end}).a5   ;
    a5_x = @ (x) a5_u.*(1+a5_1.*exp(-a5_2.*x/depth_basis_local_loss_parameter));      
    
    % data strutc for later saving
    a5_distance.a5_u = a5_u; 
    a5_distance.a5_1 = a5_1; 
    
    a5_distance.delta_a5 = depth_basis_local_loss_parameter/ a5_2; 
    % 
    delta_a5 = depth_basis_local_loss_parameter/ a5_2;  
    a5_distance.a5_x = @ (x) a5_u.*(1+a5_1.*exp(-x/delta_a5)); 
    eval(strcat('allfrequency', '.', frequency_name, '.', 'a5_distance', '=', 'a5_distance'))
    % plot a5_x 
    figure;
    hold on;
    distance = [0:0.1:20]; 
    plot(distance, a5_x(distance)) 
    xlabel( 'distance in mm', 'Interpreter', 'none','Fontsize', 16  );
    ylabel( '$a_5(x)$', 'Interpreter', 'latex', 'Fontsize', 16  );
    title(strcat('a_5(x) ', 'for frequency ', frequency_name), 'Fontsize', 16 )
    hold off;
    % saving a5_x
    save_a5_x = [save_frequency,'\' ,'a5_x'];
    if ~exist(save_a5_x,'dir')
        mkdir(save_a5_x);
    end      
    set(gcf,'position',[0.2,0.2,720,405]); 
    filename = strcat('a5_x');
    img =gcf;  
    print(img, '-dpng', '-r600', [save_a5_x,'\',convertStringsToChars(filename),'.jpg']);
    print(img, '-depsc', '-r600', [save_a5_x,'\',convertStringsToChars(filename),'.eps']);
    print(img, '-dmeta', '-r600', [save_a5_x,'\',convertStringsToChars(filename),'.emf']);
    saveas(img,[save_a5_x,'\',convertStringsToChars(filename),'.fig'])    ;      
    % save a5_x data
    data_name = 'a5_distance' ; 
    save_matfile_frequency = strcat(save_a5_x,'\',data_name,'.mat');
    save(save_matfile_frequency,'a5_distance');     
    
    
end

save_matfile_allfrequency =  fullfile(save_dir,'zusatz', 'allfrequency.mat');    
save(save_matfile_allfrequency, 'allfrequency');
close all
%% different a5(x) in one figure
save_allfrequency = fullfile(save_dir,'zusatz');
save_a5_x_frequency = fullfile(save_allfrequency, 'a5_x_frequency'); 
if ~exist(save_a5_x_frequency,'dir')
    mkdir(save_a5_x_frequency);
end

figure;
hold on; 
xlabel( 'distance in mm', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$a_5(x)$', 'Interpreter', 'latex', 'fontsize',12 );
depth_vs_frequency = []; 
frequency_vs_depth = [];
for j = 1: number_frequency
 
    if frequency_vector(j) > frequency_max
       break 
    end
    
    frequency_name = frequency_names{j};
    a5_x = allfrequency.(frequency_name).a5_distance.a5_x   ;  
    distance = [0:0.1:20];
    plot(distance, a5_x(distance));
    depth_vs_frequency = [depth_vs_frequency allfrequency.(frequency_name).a5_distance.delta_a5 ];    
    frequency_vs_depth = [frequency_vs_depth frequency_vector(j)]; 

end
hold off; 
txt_legend = string(frequency_names(frequency_vector<=frequency_max));
legend(txt_legend, 'Location', 'Northeast')
title('$a_5(x)$ under different frequency', 'Interpreter', 'latex', 'fontsize',12)
% saving figure
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('a5_distance_all_frequency');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_a5_x_frequency ,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_a5_x_frequency, '\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_a5_x_frequency , '\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_a5_x_frequency,'\',convertStringsToChars(filename),'.fig'])    ; 
close all; 

figure;
hold on
xlabel( 'frequency in Hz', 'Interpreter', 'none', 'fontsize',12 );
ylabel( '$\delta_{a5}$ in mm', 'Interpreter', 'latex', 'fontsize',12 );
plot(frequency_vs_depth, depth_vs_frequency, 'o--r', 'Linewidth', 2)
title('Einflusstiefe von a5 in dependency of frequency')
ylim([0 5])
hold off
set(gcf,'position',[0.2,0.2,720,405]); 
filename = strcat('a5_Einflusstiefe');
img =gcf;  
% print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
print(img, '-dpng', '-r600', [save_a5_x_frequency ,'\',convertStringsToChars(filename),'.jpg']);
print(img, '-depsc', '-r600', [save_a5_x_frequency, '\',convertStringsToChars(filename),'.eps']);
print(img, '-dmeta', '-r600', [save_a5_x_frequency , '\',convertStringsToChars(filename),'.emf']);
saveas(img,[save_a5_x_frequency,'\',convertStringsToChars(filename),'.fig'])    ; 
close all; 
%% return to the main dierectory
eval(['cd ' main_dir])