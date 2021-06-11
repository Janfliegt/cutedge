function [] = validation_cut_edge_model(processed_data_dir, validate_width, cur_dir)

%% 
numbins = 1000;  % how many intervals used to calculate the integration of eta
%%
% To run the script, cut_edge_model and allfrequency.mat should first be loaded to workspace
load(fullfile(processed_data_dir, 'plot_cut_edge_model_cor_data', 'cut_edge_model.mat'))
load(fullfile(processed_data_dir, 'plot_corrected_data', 'allfrequency.mat'))
test_breite_str = strcat('B', num2str(validate_width), 'mm');
test_breite_value = validate_width;
save_dir = fullfile(processed_data_dir,strcat('validate_cut_edge_', num2str(validate_width), 'mm'));
%%
hzs = fieldnames(cut_edge_model);   % cell array to save different frequency
num_hz = length(hzs);               % number of frequency to be validated

cut_edge_model_hz = cell([1, num_hz]);
allfrequency_hz = cell([1, num_hz]);
for i = 1:num_hz
    cut_edge_model_hz = cut_edge_model.(hzs{i});  % generate Field names from variables by using parentheses https://www.mathworks.com/help/matlab/matlab_prog/generate-field-names-from-variables.html
    allfrequency_hz = allfrequency.(hzs{i});
    % Messdaten
    mu_mess = allfrequency_hz.(test_breite_str).mu_relative; % original measurement, without any correction
    H_mess = allfrequency_hz.(test_breite_str).H;            % original measurement points
    
    eta_fun = cut_edge_model_hz.eta.fitting_function;        % eta(x) is introduced to later calculation
    
    H_predict = cut_edge_model_hz.H_interp_pchip_final;      % interpolation points in H in the cut edge model
    delta_mu_c_H = cut_edge_model_hz.delta_mu_interp_pchip_final;   % the corresponding delta_mu_c_H in the cut edge model
    mu_no_degradation = cut_edge_model_hz.mu_no_degreation_pchip;   % the corresponding relative mu of undamaged Eisenmaterialien
    num_H = length(H_predict);                                      % number of interpolation points
    
    test_vec = linspace(0, test_breite_value/2, numbins);           % numbins points are used to numericly calculate the F(b)from the integral of eta(x)
    eta = eta_fun(test_vec);                                        % interpolated values in eta
    Fb_fit = trapz(test_vec, eta)/(test_breite_value/2);            % integral by using trapez
    
    mu_predict = mu_no_degradation - delta_mu_c_H .* Fb_fit;
    
    mu_predict_interp = interp1(H_predict, mu_predict, H_mess, 'linear', 'extrap');
    mu_diff = mu_predict_interp - mu_mess;
    mu_diff_rela = abs((mu_predict_interp - mu_mess)./mu_mess);
    res = table(H_mess, mu_mess, mu_predict_interp, mu_diff, mu_diff_rela);
    filename = hzs{i};
    save_path = [save_dir, '/', filename, '/data/'];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    save([save_path, '/data.mat'], 'res')
    
    % plot(H_predict, mu_predict, H_mess, mu_mess)
    plot(H_mess, mu_predict_interp, 'b--', 'LineWidth', 2)
    hold on;
    plot( H_mess, mu_mess,'ro' ,'LineWidth', 2)
    hold off;
    xlabel("$H$ in A/m",'fontsize',16,'interpreter','latex')
    ylabel('relative permeability','fontsize',16,'interpreter','latex')
    h = legend('$\mu_{\rm predicted}$', '$\mu_{\rm measured}$');
    % title(hzs{i}, 'Interpreter', 'latex', 'fontsize',12)
  
    set(h, 'fontsize',16, 'Interpreter','latex');  
    set(gca,'TickLabelInterpreter','latex')  % f¨¹r die Zahl an der Achse
    set(gca,'FontSize',16) % Achse Fontsize

    filename = hzs{i};
    save_path = [save_dir, '/', filename, '/figures/'];
    if ~exist(save_path, 'dir')
        mkdir(save_path);
    end
    % sava the figure
    img =gcf;  
    set(gcf,'position',[0.2,0.2,720,405])

    print(img, '-dpng', '-r600', [save_path,'/',filename,'.jpg'])
    print(img, '-dpng', '-r600', [save_path,'/',filename, '.png'])
    print(img, '-depsc', '-r600', [save_path,'/',filename, '.eps'])
    print(img, '-dmeta', '-r600', [save_path,'/',filename, '.emf']);    
    saveas(gcf, [save_path,'/',filename, '.fig'])
end
close all
eval(['cd ' cur_dir])
    
