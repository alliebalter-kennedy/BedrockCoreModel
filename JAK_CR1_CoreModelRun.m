%% Run bedrock core model for 18JAK-CR1

clear all
close all

%% Location information

% Latitude and elevation of site
loc.lat = 69.2308;
loc.long = -49.8089;
loc.elv = 93;
corename = '18JAK-CR1';

%% Define Model Parameters

% add constraints from Holocene. This is important for calculating LIA
% erosion.

history.deglac_t = 7508; % local deglaciation age
history.historical_cover = 213; % duration of historical cover
history.historical_deglac = 10; % recent deglaciation (years before 2018)

% steady state erosion rate

SSerate = 5; % m/myr; steady state erosion rate to initialize model.

%% set length of glacials and interglacials using benthic d18O

glaciation_threshold = [3.3:0.1:4.0]; % d18O values [permil] for creating exposure histories

for ii = 1:length(glaciation_threshold)

smoothing_time = 15; % kyr, 1/2 time of smoothing window (i.e., value of 15 will be 30 kyr smoothing window)

[history.model_times, history.glacial_lengths, history.interglacial_lengths] = get_oxygen_times(history, glaciation_threshold(ii), smoothing_time); % create pre-deglaciation exposure history
%% Set model space for erosion

erosion.gl_erosion = [0:0.01:0.8]; % [mm/yr]; erosion rates to use during ice cover prior to deglaciation (constant for all glacial periods)
erosion.historical_erosion = [0:0.01:0.6];  % [mm/yr]; erosion rates to use during historical cover

% Define a rock core depth -- make it longer than your actual core. 

erosion.dz = 0.1; % [cm] needs to be same precision as your sample thickness measurements.
erosion.profile_vec = [0:erosion.dz:500]; % 5 meters - actual core is 4 m. 


%% steady erosion start

% find correct pre-made SS erosion profiel for starting cosmogenic nuclide
% concentrations

formatfileSS = 'SteadyStateProfiles/SteadyStateErosionRate_%d.mat';
filenameSS = sprintf(formatfileSS, SSerate);
steady_state = load(filenameSS); 
steady_state.N_start(size(steady_state.N_start)+1:4000001) = linspace(steady_state.N_start(end), 0, 4000001-(size(steady_state.N_start,1)));

steady_state.flag = 0; % 1 for use steady state profile at start, 0 for startinc concentrations of zero


%% Run model

tic
[N_final, depths_final] = CoreModel(loc, history, erosion, steady_state);
toc

%% Check that all final top depths are zero

if sum(sum(sum(depths_final(1, :, :)))) == 0
else 
    display('WARNING: not all depth profiles start at 0 cm');
end

%% Sample Info

dataname = '18JAK-CR1_dataanalysis.xlsx';

% sample info will read into table from excel or .txt file. Columns should
% be in order: Sample ID, Nuclide, top depth, bottom depth, quartz weight,
% atoms/g, and error (atoms/g). See example spreadsheet for formatting. 

samples = readtable(dataname);

% convert table to structural array

samples = table2struct(samples);

% extract sample depths

for a = 1:length(samples)
    samples_td(a) = samples(a).TopDepth.*2.65;
    samples_bd(a) = samples(a).BottomDepth.*2.65;
    N_samples(a) = samples(a).Nconc;
    sample_error(a) = samples(a).ErrorConc;
    avg_depth(a) = mean([samples_td(a) samples_bd(a)]);
end

%% Calculate Concentration at Sample Depths

depth_profile = depths_final(:, 1, 1); % depths from model output -- should be same as input, but if not unit test above will break. 

for a = 1:length(samples)
    [~, min_td] = min(abs(depth_profile - samples_td(a))); % maps sample depths to depth profile in case precision is different
    [~, min_bd] = min(abs(depth_profile - samples_bd(a)));
    for b = 1:length(erosion.gl_erosion)
        for c = 1:length(erosion.historical_erosion)
            N_temp = N_final(min_td:min_bd, b, c);
            depth_temp = depth_profile(min_td:min_bd);
            N_model_temp(a, b, c) = trapz(depth_temp, N_temp)./(samples_bd(a)-samples_td(a)); % sum modeled nuclide concentrations over sample depth
        end
    end
end

N_model(:, :, :, ii) = N_model_temp;
end

%% Save output

if steady_state.flag == 1
    FormatFile = '/ModelOutputs/JAK_CoreModel_%0.1f-%0.1fpermil_GlErate_%0.01f-%0.01f_LIAErate_%0.01f-%0.01f_%dSS.mat';
    filename = sprintf(FormatFile, round(glaciation_threshold(1), 2), round(glaciation_threshold(end), 2), round(erosion.gl_erosion(1), 2), ...
    round(erosion.gl_erosion(end), 2), round(erosion.historical_erosion(1), 2), round(erosion.historical_erosion(end), 2), SSerate);
elseif steady_state.flag == 0
    FormatFile = '/ModelOutputs/JAK_CoreModel_%0.1f-%0.1fpermil_GlErate_%0.01f-%0.01f_LIAErate_%0.01f-%0.01f_zeroconcstart.mat';
    filename = sprintf(FormatFile, round(glaciation_threshold(1),2), round(glaciation_threshold(end),2), round(erosion.gl_erosion(1), 2), ...
    round(erosion.gl_erosion(end), 2), round(erosion.historical_erosion(1), 2), round(erosion.historical_erosion(end), 2));
end

save(filename, 'N_model', 'erosion', 'glaciation_threshold');

beep on