function [Best_Result, Best_Solution, Fitness_xbest_evolution] = evolutionary(upper, number_width, ITE, RES)
%% evolutionary algorithm
    
    number_designvariable = number_width-1;                 % number of F(b) 4, 5, 7.5, 10, 15, 30, 60, without F(120) 
    range = [upper* ones(1,number_designvariable);zeros(1,number_designvariable)];
    % initiation
    P_Size = number_designvariable*5;         % Population Size
    VAR = number_designvariable;              % Number of design Variables
%     ITE = 12;              % Anzahl der Iterationen
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
%     RES = 4; % restarten mal              % 10-mal Optimierungsversuche
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
%                            diff_child = diff(child);
%                            if nnz(diff_child > 0) > 0
%                                continue;
%                            end


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