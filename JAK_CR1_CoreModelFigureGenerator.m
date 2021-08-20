% This script generates figures for the core erosion model.
clear all
close all

% addpath('/Users/alexandrabalter/MATLAB/projects/bedrock_core/Data')
addpath('/SteadyStateProfiles')
addpath('/ModelOutputs')

% do you want to plot all modeled Be-10 results on Figure 2?
Fig2Flag = 0; % 0 for no, 1 for yes.


%%
none = load('JAK_CoreModel_3.3-4.0permil_GlErate_0.0-0.8_LIAErate_0.0-0.6_zeroconcstart.mat');
fifty = load('JAK_CoreModel_3.3-4.0permil_GlErate_0.0-0.8_LIAErate_0.0-0.6_50SS.mat');
ten = load('JAK_CoreModel_3.3-4.0permil_GlErate_0.0-0.8_LIAErate_0.0-0.6_10SS.mat');
five = load('JAK_CoreModel_3.3-4.0permil_GlErate_0.0-0.8_LIAErate_0.0-0.6_5SS.mat');
%% Sample Info

filename = 'Users/alexandrabalter/MATLAB/projects/bedrock_core/Data/18JAK-CR1_dataanalysis.xlsx';

% sample info will read into table from excel or .txt file. Columns should
% be in order: Sample ID, Nuclide, top depth, bottom depth, quartz weight,
% atoms/g, and error (atoms/g). See example spreadsheet for formatting. 

samples = readtable(filename);

% convert table to structural array

samples = table2struct(samples);

%% Get sample concentration vector

for a = 1:length(samples)
    sample_conc(a, 1) = samples(a).Nconc;
    sample_error(a, 1) = samples(a).ErrorConc;
    avg_depth(a, 1) = mean([samples(a).TopDepth samples(a).BottomDepth]);
end

%% Loop through scenarios

for a = 1:4
    
 % call correct steady-state erosion scenario   
    if a == 1
        model = five;
    elseif a == 2
        model = ten;
    elseif a == 3
        model = fifty;
    elseif a == 4
        model = none;
    end
    
 % loop through glaciation thresholds
 
    for b = 1:length(model.glaciation_threshold)  
        
        N_model = model.N_model(:, :, :, b);

% Calculate chi-2 statistic
        
    mismatch = squeeze(sum(((sample_conc - N_model)./sample_error).^2, 1));

    reduced_chi2 = mismatch./(length(samples)-2);

    % find best-fitting combination of erosion rates for each d18O threshold.
    minimum = min(min(reduced_chi2));

    [minGL, minLIA] = find(reduced_chi2 == minimum);
    
% plot chi-squared values as colormap
    
    plot_index = sub2ind([numel(model.glaciation_threshold), 4], b, a);
    figure(1)
    subplot(5, length(model.glaciation_threshold), plot_index+length(model.glaciation_threshold))
    
    % plot glacial erosion on y-axis and LIA erosion on x-axis (NEED TO CHECK). Color by
    % mismatch
        rho = 2.65;

        [X, Y] = meshgrid(model.erosion.historical_erosion, model.erosion.gl_erosion);

        pcolor(model.erosion.historical_erosion, model.erosion.gl_erosion, reduced_chi2);
        colormap(parula);
        c = colorbar;
        caxis([0 3]);
        c.Label.String = 'Red-chi^{2}';
%         c.Label.FontSize = 14;
        set(gca, 'xlim', [0 0.6], 'ylim', [0 0.8])
        xlabel('Hist. e-rate (mm/yr)');
        ylabel('Pleisto. e-rate (mm/yr)')
        
        shading flat
        
        txt = '$r\\chi^{2}$ = %0.2f\nPl. = %0.2f\nHist. = %0.2f';
        bestfit = sprintf(txt, minimum, model.erosion.gl_erosion(minGL),model.erosion.historical_erosion(minLIA));
        hold on
        text(0.1, 0.6, bestfit, 'interpreter', 'latex', 'FontName', 'Arial');
        
     %% Get exposure history from d18O threshold
     
        history.deglac_t = 7508;
        
        [model_times, glacial_lengths, interglacial_lengths, time, switch_times, oxygen_isotopes] = ...
            get_oxygen_times(history, model.glaciation_threshold(b), 15); 


 %% plot d18o threshold and resulting exposure histories
        
        figure(1)
        hold on
        subplot(5, length(model.glaciation_threshold), b)
            plot(time/1000, oxygen_isotopes, 'k');
            hold on
            plot([time(1)/1000 time(end)/1000], [model.glaciation_threshold(b) model.glaciation_threshold(b)], 'r');
            hold on

            for i = 1:length(switch_times)-1
                if mod(i, 2) == 1
                    plot([switch_times(i)/1000 switch_times(i+1)/1000], [2.8 2.8], '-', 'Color', [135/256, 206/256, 235/256], 'LineWidth', 8)
                else
                    plot([switch_times(i)/1000 switch_times(i+1)/1000], [2.8 2.8], '-', 'Color', [128/256,0,0],  'LineWidth', 8)
                end
                hold on
            end
            
            txt2 = '$\\delta^{18}O$ = %0.1f';
            gl_thresh = sprintf(txt2, model.glaciation_threshold(b));
            
            hold on
            text(1e3, 4.7, gl_thresh, 'Interpreter', 'Latex', 'FontName', 'Arial');
            
            set(gca, 'YDir', 'Reverse', 'YLim', [2.5 5], 'XLim', [0 2.7e3])
            xlabel('Time (kyr)')
            ylabel('\delta^{18}O (permil)')
    
 
 %% Plot modeled Be-10 concentrations for best-fitting erosion rates
 
 if Fig2Flag == 1
 model_conc = N_model(:, minGL, minLIA);
 
  figure(2) 
  subplot(4, length(model.glaciation_threshold), plot_index+length(model.glaciation_threshold))
  errorbar(sample_conc, avg_depth, sample_error, 'ko', 'horizontal');
    hold on
    plot(model_conc, avg_depth, 'ro');

    set(gca, 'YDir', 'reverse', 'XScale', 'log')
    xlabel('^{10}Be conc. (at. g^{-1})')
    ylabel('Depth (cm)')
%     legend('Measured conc.', 'Modeled Concentration', 'Location', 'Southeast')
    grid on
    
    hold on
    bestfit = sprintf(txt, minimum, model.erosion.gl_erosion(minGL),model.erosion.historical_erosion(minLIA));
        hold on
        text(2.5e3, 300, bestfit, 'Interpreter', 'Latex', 'FontName', 'Arial');
    
 %% plot d18o on fig 2
        
        figure(2)
        hold on
        subplot(4, length(model.glaciation_threshold), b)
            plot(time/1000, oxygen_isotopes, 'k');
            hold on
            plot([time(1)/1000 time(end)/1000], [model.glaciation_threshold(b) model.glaciation_threshold], 'r');
            hold on

            for i = 1:length(switch_times)-1
                if mod(i, 2) == 1
                    plot([switch_times(i)/1000 switch_times(i+1)/1000], [2.8 2.8], '-', 'Color', [135/256, 206/256, 235/256], 'LineWidth', 8)
                else
                    plot([switch_times(i)/1000 switch_times(i+1)/1000], [2.8 2.8], '-', 'Color', [128/256,0,0],  'LineWidth', 8)
                end
                hold on
            end
            hold on 
            text(1e3, 4.7, gl_thresh, 'interpreter', 'latex');
            set(gca, 'YDir', 'Reverse', 'YLim', [2.5 5], 'XLim', [0 2.7e3])
            xlabel('Time (kyr)')
            ylabel('\delta^{18}O (permil)')
 else 
 end
    end
        
end
% 
