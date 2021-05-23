function [] = optimierung_main_find_Fb_corrected_data_for_all_frequency(original_data_dir, used_frequency, cur_dir)
eval(['cd ' original_data_dir])
load('plot_corrected_data\allfrequency.mat');
frequency_fields = fieldnames(allfrequency);
number_of_frequency = length(frequency_fields);
for i = 1 : number_of_frequency
    frequency_name = frequency_fields{i};
    cur_frequency = str2num(frequency_name(2:end-2));
    if ~ismember(cur_frequency, used_frequency)
        allfrequency = rmfield(allfrequency, frequency_name);
    end
end
    
C = who; 
struct_name= C{:}; 
struct_cut_edge = 'cut_edge_model'; 
save_dir = 'plot_cut_edge_model_cor_data';
mkdir(save_dir)

% frequency name vector extraction
frequency_names = fieldnames(eval(struct_name));
number_of_frequency = length(frequency_names) ;

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
for i_frequency = 1 : number_of_frequency
    frequency_name = frequency_names{i_frequency} 
    
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
            width_vector(j) = width_tmp ;                                      %  width vector are saved in width_vector
        end      

    end
    % fitting method for F(b) after discrete F(b) are determined
    Fb_fitting_method = 2; % 1 for smoothing spline; 2 for a*exp(x*b); 3 for a*x^b+c
    %% evolutionary algorithm
    number_designvariable = number_width-1;                 % number of F(b) 4, 5, 7.5, 10, 15, 30, 60, without F(120) 
    range = [ones(1,number_designvariable);zeros(1,number_designvariable)];
    % initiation
    P_Size = number_designvariable*5;         % Population Size
    VAR = number_designvariable;              % Number of design Variables
    ITE = 12;              % Anzahl der Iterationen
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
    print_DE= 0;
    n0 = 2;                                % 2 Kinder für ein Elter
    x_best = zeros(1,VAR);                 % Speicherung der Designvariablen der besten Lösung in der Population in x_best 
    Fitness_x_best = zeros(1);             % Speicherung der Eigenschaft der besten Lösung, inklusive Zielgröße, Verletzunggrad der Nebenbedingung, feasibility
    RES = 4; % restarten mal              % 10-mal Optimierungsversuche
    Fitness_xbest_evolution = zeros(RES,ITE); % Speicherung die Verläufe der Bewertung der besten Lösung in jedem Optimierungsversuch

    Best_Result = zeros(RES,1);            % Speicherung der letztendlichen Fitnesserwerte der besten Lösung nach jedem Optimierungsversuch
    Best_Solution = zeros(RES,VAR);        % Speicherung der letztendlichen Designvariabllen der besten Lösung nach jedem Optimierungsversuch



    % RES restarting of evolutionary algorithm
    for ggg = 1:RES
        % initiation der Population innerhalb der Grenze
        for j = 1: P_Size
            j;
                   while 1
                       Population(j,1) = (range(1,1)-range(2,1))*rand + range(2,1); 
                       for nn = 2 : VAR
                           Population(j,nn) = (Population(j,nn-1)-range(2,nn))*rand + range(2,nn);          % F(b) decreases monoton
                       end

                        [Temp1,Temp2,Temp3]= modeling_corrected_data(Population(j, :), print_DE);  
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
                           diff_child = diff(child);
                           if nnz(diff_child > 0) > 0
                               continue;
                           end


                            [Temp1,Temp2,Temp3]= modeling_corrected_data(child, print_DE);  
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
    %% best results for F(b1)....F(bn-1) bestimmen, und die entsprechende delta_mu_c_H
    print_DE = 1;
    [z,y] = min(Best_Result);               % Auswahl der besten Lösung innherhalb 10-mal Optimierungsversuche
    xx = Best_Solution(y,:);
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

save_matfile_allfrequency =  strcat(save_dir,'\', struct_cut_edge,'.mat');    
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
saveas(img,[save_dir,'\',convertStringsToChars(filename),'.fig'])    ; 
close all
eval(['cd ' cur_dir])
end
function [cost,y3,y2] = modeling_corrected_data(xx,print_DE)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y3 = 0;     % no constraints exist
    y2 = 0;     % no constraints exist
    F_b = xx;               % from 4,5, ..., 60mm
    % declare the global variables
    global f_in_Hz         
    global width_names
    global width_vector
    global H_interpolation
    global mu_no_degradation
    
    % define mu without degradation, corresponding to 120mm
    % measurement points in H are used
    mu_no_degradation = eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'mu_relative_interp_2'));
            
    % definition of interpolation points in H
    H_interpolation = eval(strcat('f_in_Hz', '.', width_names{end}, '.', 'H_interp_2'));  % here can be changed, if more interpolated H are desired to approximate the F(b)

    % length 
    number_width = length(width_names);
    number_interpolated_H = length(H_interpolation);
    % define matrix for f(H,bi)
    global f_H_b_matrix
    f_H_b_matrix = zeros(number_interpolated_H, number_width-1);  % each row for a H, each column corresponds to a width
    % filling f_H_b_matrix
    for k = 1 : number_interpolated_H   % iteration along H
        for g = 1: number_width-1       % iteration along width
            F_b_tmp = F_b(g);
            mu_without_degradation_tmp = mu_no_degradation(k);
            % mu under other width will be interplolated
            mu_tmp = eval(strcat('f_in_Hz', '.', width_names{g}, '.', 'mu_relative_interp_2'));
            H_tmp = eval(strcat('f_in_Hz', '.', width_names{g}, '.', 'H_interp_2')); 
            mu_with_degradation_vector = interp1(H_tmp, mu_tmp, H_interpolation, 'pchip');
            % mu under the same H is saved in mu_measure_with_degradation_tmp
            mu_measure_with_degradation_tmp = mu_with_degradation_vector(k); 
            % 
            f_H_b_tmp = (mu_without_degradation_tmp - mu_measure_with_degradation_tmp)/ F_b_tmp / mu_without_degradation_tmp; % normalized, different from silas formula
            %  f_H_b_tmp = (mu_without_degradation_tmp - mu_measure_with_degradation_tmp)/ F_b_tmp ; % non-normalized, the silas formula
            %  f_H_b_tmp = min(f_H_b_tmp,1);                                  % f_H_b_tmp can not over 1
            %  matrix filling
            f_H_b_matrix(k,g) = f_H_b_tmp;
        end
    end
    cost = mean(var(f_H_b_matrix,0,2));
    if max(max(f_H_b_matrix)) > 1  % f_H_b_tmp should not be greater than 1, which is a constraint to be considered
        y3 = 1; 
        y2 =  max(max(f_H_b_matrix)) -1; 
    end
    
    if print_DE == 0
        return;
    end
    str = sprintf('RESULTS:\n');
    disp(str);
    str = sprintf('F_b:\n');
    disp(str);
    F_b;
    str = sprintf('cost function:\n');
    disp(str);
    cost;
    f_H_b_matrix;
    global variance_for_different_H
    variance_for_different_H = var(f_H_b_matrix,0,2);
    % delta_mu_c_H calculate 
    global delta_mu_c_H
    delta_mu_c_H =  zeros(length(H_interpolation),1); 
    % preprocess f_H_b_matrix so that no element bigger than 1 due to numeric error
    f_H_b_matrix(f_H_b_matrix > 1) = 1; 
    %
    for k = 1:length(H_interpolation)
        delta_mu_c_H_tmp = mean(f_H_b_matrix(k,:)) .* mu_no_degradation(k); % average f_H_b used to calculate delta_mu_c_H 
        delta_mu_c_H(k) = delta_mu_c_H_tmp; 
    end
end
