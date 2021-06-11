% A_parameter determination

function [] = A_parameter_determination_without_interface(processed_data_dir, hyst_art, cur_dir)
%%
data_simplified = fullfile(processed_data_dir, 'simplified_data.mat');
load(data_simplified)
width_names = fieldnames(eval('L_messung_simplified')); % width names 
struct_name = 'L_messung_simplified';
number_of_width = length(width_names); 
width_vector = zeros(number_of_width,1);

save_dir = fullfile(processed_data_dir, 'A_parameters');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

%% initialisation
%hyst_art = 1;     % hyst_art = 1 means the formula of Silas, a1.*J.^alpha.*f % hyst_art = 2 means the formula of Silas, a1.*J.^(alpha+beta.*J).*f
% characteristic paramters related to the a_2 parameter todo
d = 0.35*10^(-3); % dicke in m
Si = 3.5; % % Siliziumgehalt in [%]
rho_e = (12+11*Si)*(10^-8);     % spezifische Widerstand
rho = 7500;    % massendichte kg/m^3
% frequency range and polarisation range to determine a_5
f_min_a5 = 10; % minimal frequency to extract a_5 parameter
f_max_a5 = 50; % maximal frequency to extract a_5 parameter
J_min_a5 = 0.3;   % minimal Polarisation to extract a_5 parameter
J_max_a5 = 1.3; % maximal Polarisation to extract a_5 parameter

% frequency range and polarisation range to determine a_3 and a_4
% here all the measure data are used to fit a_3 a_4 parameters
f_min_a3_a4 = 100; %0; % minimal frequency to extract a_3 a_4 parameter
f_max_a3_a4 = inf; % maximal frequency to extract a_3 a_4 parameter
J_min_a3_a4 = 1.3; %0; % minimal polarisation to extract a_3 a_4 parameter
J_max_a3_a4 = inf; % maximal frequency to extract a_3 a_4 parameter

%% determination of width vector
for j = 1 : number_of_width
    width_name_txt =  width_names{j};
    width_tmp = regexp(width_name_txt, '\d*\d*', 'match');
    if length(width_tmp) == 1
        width_tmp = width_tmp{:};
        width_tmp = str2double(width_tmp);
        width_vector(j) = width_tmp;
    else       
        width_tmp= strcat(width_tmp{1}, '.', width_tmp{end});
        width_tmp = str2double(width_tmp);
        width_vector(j) = width_tmp;                                      %  width vector are saved in width_vector
    end      

end

%% loops to determine the A-parameters for different width
for i = 1 : number_of_width
    width_name = width_names{i};
    var_name = strcat(struct_name, '.', width_name); 
    frequency_names = fieldnames(eval(var_name));       % measured frequencies under each width
    number_frequency = length(frequency_names); 
    frequency_vector = zeros(number_frequency,1);
    
    % directory creation
    
    filename_width = width_name;
    save_width_path = fullfile(save_dir,filename_width);
    if ~exist(save_width_path,'dir')
        mkdir(save_width_path);
    end    
    
    % determination of the frequency vector under each width 
    for j = 1 : number_frequency
        frequency_name_txt =  frequency_names{j} ;
        frequency_tmp = regexp(frequency_name_txt, '\d*\d*', 'match');
        if length(frequency_tmp) == 1
            frequency_tmp = frequency_tmp{:};
            frequency_tmp = str2double(frequency_tmp);
            frequency_vector(j) = frequency_tmp;
        else       
            frequency_tmp= strcat(frequency_tmp{1}, '.', frequency_tmp{end});
            frequency_tmp = str2double(frequency_tmp);
            frequency_vector(j) = frequency_tmp ;                                      %  width vector are saved in width_vector
        end      

    end
 
%% a_2 parameter determination    
% Blechdaten = inputdlg({'Blechdicke eingeben [mm]:', 'Siliziumgehalt eingeben [%]:', 'Dichte des Materials eingeben [kg/m^3]:'}, 'Blechdaten eingeben');
% d = str2double(Blechdaten(1))*10^(-3); % dicke in m
% Si = str2double(Blechdaten(2)); % % Siliziumgehalt in [%]
% rho_e = (12+11*Si)*(10^-8);     % spezifische Widerstand
% rho = str2double(Blechdaten(3));    % massendichte kg/m^3
% a2 = ((pi^2)*(d^2))/(6*rho*rho_e);

a2 = ((pi^2)*(d^2))/(6*rho*rho_e);

%% a_1 und alpha determination
frequency_min = min(frequency_vector);      % the measurement with the minimal frequency is used to determine hysterese-parameters 
frequency_min_tex = ['f', num2str(frequency_min), 'Hz'];
data_tmp = eval(strcat(struct_name,  '.', width_name, '.',frequency_min_tex)); % data read
J = data_tmp.Jmax_T; 
f =   data_tmp.Frequenz_Hz; 
P_fe = data_tmp.Ps_W_kg; 

Peddy = a2.*J.^2.*f.^2; 
P_hyst = P_fe - Peddy;             % extraction of the eddy current loss, and the rest is the hysteresis loss
P_hyst_normalized = P_hyst./ f;    % divided by frequency for prepraring for the fitting in the following

% hyst_art = inputdlg({'Art der Hystereseverlust [1 or 2]:',  }, 'Hystereseverlustart eingeben');
% hyst_art = str2double(hyst_art(1)); 

if hyst_art == 1
    %Fitting
    ini = [0.01,1.5];
    lowbound = zeros(1,2);
    HystFunktion = @(HystOpt, J) (HystOpt(1) .* J .^(HystOpt(2)));
    HystOpt = lsqcurvefit(HystFunktion, ini, J, P_hyst_normalized, lowbound);

    %Verlustvektor mit den ermittelten Parametern fuer das plotten berechnen
    Physt_berechnet = HystOpt(1) .* J .^(HystOpt(2)) .*f;
    figure; 
    plot(J, P_fe, 'r-')
    hold on; 
    plot(J, Physt_berechnet,'b')
    hold on; 
    txt_legend = ["measured P_{hyst}", "calculated P_{hyst}"]; 
    xlabel( 'J in T', 'Interpreter', 'none', 'fontsize', 15 );
    ylabel( '$P_{\rm hyst}$ in W/kg', 'Interpreter', 'latex', 'fontsize', 15 );
    legend(txt_legend, 'Location', 'northwest', 'fontsize', 15)
    a1 = HystOpt(1); 
    alpha = HystOpt(2);
    set(gca,'FontSize',15)
    hold off; 
else            
    %Fitting
    ini = [0.01,1.5, 0.1];
    lowbound = zeros(1,3);
    HystFunktion = @(HystOpt, J) (HystOpt(1) .* J .^(HystOpt(2)+HystOpt(3).*J));
    HystOpt = lsqcurvefit(HystFunktion, ini, J, P_hyst_normalized, lowbound);

    %Verlustvektor mit den ermittelten Parametern fuer das plotten berechnen
    Physt_berechnet = HystOpt(1) .* J .^(HystOpt(2)+HystOpt(3).*J) .*f;
    figure; 
    plot(J, P_fe, 'r-')
    hold on; 
    plot(J, Physt_berechnet,'b')
    hold on; 
    txt_legend = ["measured P_{hyst}", "calculated P_{hyst}"]; 
    xlabel( 'J in T', 'Interpreter', 'none', 'fontsize', 15 );
    ylabel( '$P_{\rm hyst}$ in W/kg', 'Interpreter', 'latex', 'fontsize', 15 );
    legend(txt_legend, 'Location', 'northwest', 'fontsize', 15)
    a1 = HystOpt(1); 
    alpha = HystOpt(2);
    beta = HystOpt(3) ; 
    set(gca,'FontSize',15)
    hold off;     
end


%% a_5 parameter determination
% a5_measured_daten = inputdlg({'frequency min [Hz]:', 'frequency max [Hz]:', 'Polarisation min [T]:', 'Polarisation max [T]:'}, 'Range in f and J to cal. a5');
% f_min_a5 = str2double(a5_measured_daten(1)); % dicke in m
% f_max_a5 = str2double(a5_measured_daten(2)); % dicke in m
% J_min_a5 = str2double(a5_measured_daten(3)); % dicke in m
% J_max_a5 = str2double(a5_measured_daten(4)); % dicke in m

index_frequency = (frequency_vector >= f_min_a5 & frequency_vector <= f_max_a5);  % filter by frequency
frequency_range_a5 = frequency_vector(index_frequency); % only these frequencies are used to fit a_5 parameters
% initialization with empty array
J_a5 = []; 
f_a5 = []; 
P_fe_a5 = []; 
J = []
f = []
P_fe = []
% data extraction 
for k = 1 : length(frequency_range_a5)
    frequency_text = ['f', num2str(frequency_range_a5(k)), 'Hz']; 
    data_tmp = eval(strcat(struct_name,  '.', width_name, '.',frequency_text)); 
    J = data_tmp.Jmax_T; 
    f =   data_tmp.Frequenz_Hz; 
    P_fe = data_tmp.Ps_W_kg; 
    % filter by the range of J 
    index_tmp = (J >= J_min_a5 & J <= J_max_a5); % filter by the range of polarisation
    J_filterd = J(index_tmp); 
    f_filterd = f(index_tmp); 
    P_fe_filterd = P_fe(index_tmp); 
    % saved in J_a5
    J_a5 = [J_a5; J_filterd]; 
    f_a5 = [f_a5; f_filterd];
    P_fe_a5 = [P_fe_a5; P_fe_filterd];
end

if hyst_art == 1
    Physt = a1.*J_a5.^(alpha).*f_a5; 
else
    Physt = a1.*J_a5.^(alpha+beta.*J_a5).*f_a5 ; 
end

Peddy = a2*(J_a5.^2).*(f_a5.^2);
Pexc = P_fe_a5 - Physt - Peddy;         % the excess loss is the difference between the measured loss and the sum of calculated hysteris and eddy current
Pexc_normalized = Pexc ./ (f_a5.^1.5);  % normalization before fitting      
%Fitting
lowbound = 0;
ini = 0.0001;
ExzessFunktion = @(ExzessOpt, J) (ExzessOpt .* J.^1.5 );
ExzessOpt = lsqcurvefit(ExzessFunktion, ini, J_a5, Pexc_normalized, lowbound);
a5 = ExzessOpt; 
%   ExzessFunktion = @(ExzessOpt, J) (ExzessOpt .* J.^1.5 .* f_a5.^1.5  );
%   ExzessOpt = lsqcurvefit(ExzessFunktion, ini, J_a5, Pexc, lowbound);
%   a5 = ExzessOpt
% 
% Pexc_berechnet = a5 .* J_a5.^1.5 .* f_a5.^1.5 ;
Pexc_berechnet_normalized = a5 .* J_a5.^1.5;
figure; 
plot(J_a5, Pexc_berechnet_normalized, 'r-')
hold on; 
plot(J_a5, Pexc_normalized, 'b')
txt_legend = [ "calculated P_{exc, normalized}" , "measured P_{exc, normalized}"]; 
xlabel( 'J in T', 'Interpreter', 'none', 'fontsize', 15 );
ylabel( '$P_{\rm exc}/f^{1.5}$ in W/(kg$\cdot$Hz)', 'Interpreter', 'latex', 'fontsize', 15 );
legend(txt_legend, 'Location', 'northwest', 'fontsize', 15)
set(gca,'FontSize',15)
hold off;    

%% a3 und a4 parameters determination
% a3_a4_measured_daten = inputdlg({'frequency min [Hz]:', 'frequency max [Hz]:', 'Polarisation min [T]:', 'Polarisation max [T]:'}, 'Range in f and J to cal. a5');
% f_min_a3_a4 = str2double(a3_a4_measured_daten(1)); % dicke in m
% f_max_a3_a4 = str2double(a3_a4_measured_daten(2)); % dicke in m
% J_min_a3_a4 = str2double(a3_a4_measured_daten(3)); % dicke in m
% J_max_a3_a4 = str2double(a3_a4_measured_daten(4)); % dicke in m

index_frequency = (frequency_vector >= f_min_a3_a4 & frequency_vector <= f_max_a3_a4); % filter
frequency_range_a3_a4 = frequency_vector(index_frequency); 
% initialization with empty array
J_a3_a4 = []; 
f_a3_a4 = []; 
P_fe_a3_a4 = []; 
J = []
f = []
P_fe = []
% data extraction 
for k = 1 : length(frequency_range_a3_a4)
    frequency_text = ['f', num2str(frequency_range_a3_a4(k)), 'Hz']; 
    data_tmp = eval(strcat(struct_name,  '.', width_name, '.',frequency_text)); 
    J = data_tmp.Jmax_T; 
    f =   data_tmp.Frequenz_Hz; 
    P_fe = data_tmp.Ps_W_kg; 
    % filter of the range of J 
    index_tmp = (J >= J_min_a3_a4 & J <= J_max_a3_a4); 
    J_filterd = J(index_tmp); 
    f_filterd = f(index_tmp); 
    P_fe_filterd =P_fe(index_tmp); 
    % saved in J_a5
    J_a3_a4 = [J_a3_a4; J_filterd]; 
    f_a3_a4 = [f_a3_a4; f_filterd];
    P_fe_a3_a4 = [P_fe_a3_a4; P_fe_filterd];
end

if hyst_art == 1
    Physt = a1.*J_a3_a4.^(alpha).*f_a3_a4; 
else
    Physt = a1.*J_a3_a4.^(alpha+beta.*J_a3_a4).*f_a3_a4 ; 
end    

Peddy = a2.*(J_a3_a4.^2).*(f_a3_a4.^2);
Pexc = a5.*(J_a3_a4.^1.5).*(f_a3_a4.^1.5);
Psat = P_fe_a3_a4 - Physt - Peddy - Pexc;   % the saturation loss is the difference between the measured loss and the sum of calculated hys, eddy, and excess losses
Psat_normalized = Psat./ (f_a3_a4.^2); 

%Fitting
lowbound = [0,0];
ini = [0.1,0.1];

if i== number_of_width
    lowbound = [0,6];
    ini = [0.1,6];    
end

% SatFunktion = @(SatOpt, J) a2*SatOpt(1).* (J.^(2+SatOpt(2)));
% SatOpt = lsqcurvefit(SatFunktion, ini, J_a3_a4, Psat_normalized, lowbound);
% a3 = SatOpt(1)
% a4 = SatOpt(2)
SatFunktion = @(SatOpt, J) a2*SatOpt(1).* (J.^(2+SatOpt(2))) .* f_a3_a4.^2;
SatOpt = lsqcurvefit(SatFunktion, ini, J_a3_a4, Psat, lowbound);
a3 = SatOpt(1)
a4 = SatOpt(2)
% plot normalized saturation loss 
Psat_berechnet_normalized = a2 .* a3 .* J_a3_a4.^(2+a4)  ;
figure; 
plot(J_a3_a4, Psat_berechnet_normalized, 'r')
hold on; 
plot(J_a3_a4, Psat_normalized,'b')
txt_legend = [ "calculated P_{sat, normalized}" , "measured P_{sat, normalized}"]; 
xlabel( 'J in T', 'Interpreter', 'none', 'fontsize', 15 );
ylabel( '$P_{\rm sat}/f^{2}$ in W/(kg$\cdot$Hz)', 'Interpreter', 'latex', 'fontsize', 15 );
legend(txt_legend, 'Location', 'northwest', 'fontsize', 15)
set(gca,'FontSize',15)
hold off;  

%% save a-parameters for each width
if hyst_art == 1
    A_parameters.(width_name).a1 = a1; 
    A_parameters.(width_name).alpha = alpha;
    A_parameters.(width_name).a2 = a2;
    A_parameters.(width_name).a3 = a3;
    A_parameters.(width_name).a4 = a4;
    A_parameters.(width_name).a5 = a5;
else 
    A_parameters.(width_name).a1 = a1; 
    A_parameters.(width_name).alpha = alpha;
    A_parameters.(width_name).beta = beta;
    A_parameters.(width_name).a2 = a2;
    A_parameters.(width_name).a3 = a3;
    A_parameters.(width_name).a4 = a4;
    A_parameters.(width_name).a5 = a5;    
end

%% save plots calculated eisenverluste vs measured eisenverluste
for j = 1 : number_frequency
    frequency_text =  frequency_names{j} ;
    data_tmp = eval(strcat(struct_name,  '.', width_name, '.',frequency_text));  % data extraction for each frequency
    J = data_tmp.Jmax_T; 
    f =   data_tmp.Frequenz_Hz; 
    P_fe_measured = data_tmp.Ps_W_kg;    

    if hyst_art == 1
        Physt = a1.*J.^(alpha).*f; 
    else
        Physt = a1.*J.^(alpha+beta.*J).*f; 
    end  

 %  Physt = a1.*J.^(alpha+beta.*J)*f;
    Peddy = a2.*(J.^2).*(f.^2);
    Pexc = a5.*(J.^1.5).*(f.^1.5);
    Psat = a2.*a3.*(J.^(2+a4)).*(f.^2);
    Pges = Physt + Peddy + Pexc + Psat;    
    load('D:\Materials\thesis\schnittkanteeffekt\A-Parameter\02_12_2020_A_parameter_determination\RWTHColors.mat') 
    figure; 
    hold on
    grid on
    AreaYMatrix = [Peddy,Physt,Pexc,Psat];
    VerlustArea = area(J,AreaYMatrix);
    VerlustArea(1).FaceColor = RWTHblue1;
    VerlustArea(2).FaceColor = RWTHblue2;
    VerlustArea(3).FaceColor = RWTHblue3;
    VerlustArea(4).FaceColor = RWTHpetrol2;   
    xlabel('Magnetische Polarisation in T', 'FontSize', 12, 'FontName', 'Arial')
    ylabel('Eisenverluste in W/kg', 'FontSize', 12, 'FontName', 'Arial')
    title(sprintf('Verlustzusammensetzung bei %dHz', round(mean(f))), 'FontSize', 12, 'FontName', 'Arial', 'FontWeight', 'Normal')
% measured eisenverluste
    plot(J, P_fe_measured, 'Color', RWTHorange1, 'LineWidth', 1.5, 'LineStyle', '--')
    legend({'Wirbelstromverluste', 'Hystereseverluste', 'Exzessverluste', 'Sättigungsverluste', 'Messung'}, 'location', 'northwest')
    hold off;
% save figure
    filename = frequency_text;
    save_frequency_path = [save_width_path,'\',convertStringsToChars(filename)];
    if ~exist(save_frequency_path,'dir')
        mkdir(save_frequency_path);
    end
    img =gcf; 
    % print(img, '-dpng', '-r600', [save_path,'\',filename,'.png'])
    print(img, '-dpng', '-r600', [save_frequency_path,'\',convertStringsToChars(filename),'.jpg'])
    print(img, '-depsc', '-r600', [save_frequency_path,'\',convertStringsToChars(filename),'.eps'])
    print(img, '-dmeta', '-r600', [save_frequency_path,'\',convertStringsToChars(filename),'.emf'])
    saveas(img,[save_frequency_path,'\',convertStringsToChars(filename),'.fig'])     
    close all; 
end

% save A_parameters for different widths 
eval(strcat('A_parameters_', width_name, '=' ,'A_parameters.(width_name)')); 
save_matfile_width = strcat(save_width_path,'\','A_parameters','.mat');
save(save_matfile_width,strcat('A_parameters_', width_name)); 

end

save_matfile_width_all = strcat(save_dir,'\','A_parameters','.mat');
save(save_matfile_width_all,'A_parameters'); 

%% calculate the A-parameters and draw figures
widths = fieldnames(A_parameters);
num_width = length(widths);
A_parameters_cell = fieldnames(A_parameters.(widths{1}));
num_A_parameters = length(A_parameters_cell);
for i = 1:num_A_parameters
    A_parameter_cur = A_parameters_cell{i};
    A_parameter_cur_array = zeros(1, num_width);
    for j = 1:num_width
        A_parameter_cur_array(j) = A_parameters.(widths{j}).(A_parameter_cur);
    end
    
    filename = A_parameter_cur;
    plot(width_vector, A_parameter_cur_array, 'o--r', 'Linewidth', 2)  
    xlabel('width [mm]','fontsize',14,'interpreter','latex')
    ylabel(A_parameters_cell(i),'fontsize',14,'interpreter','latex')
    set(gca,'TickLabelInterpreter','latex')  % f¨¹r die Zahl an der Achse
    set(gca,'FontSize', 14)  % Achse Fontsize

    save_path = fullfile(save_dir, filename);
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    % sava the figure
    img =gcf;  
    print(img, '-dpng', '-r600', [save_path,'/',filename,'.jpg'])
    print(img, '-dpng', '-r600', [save_path,'/',filename])
    print(img, '-deps', '-r600', [save_path,'/',filename])
    saveas(gcf, [save_path,'/',filename, '.fig'])
    close all
end
close all
eval(['cd ' cur_dir])

