function [model_times, glacial_lengths, interglacial_lengths, time, switch_times, oxygen_isotopes] = get_oxygen_times(history, glaciation_threshold, smoothing_time); 

% derives exposure history from benthic oxygen stack of Lisiecki & Raymo
% (2005).

% Lisiecki, L. E., & Raymo, M. E. (2005). A Pliocene?Pleistocene stack of 
% 57 globally distributed benthic ?18O records. Paleoceanography, 20(1)
% https://doi.org/10.1029/2004pa001071

%% Load and parse dataset

load('LR04.txt'); % load LRO4 stack

LR04(LR04(:, 1) > 2700, :) = []; % cap at 2.7 Myr

oxygen_raw = LR04(:, 2); % vectorize d18O data

time = LR04(:, 1).*1e3; % vectorize time data

%% clean data

% interpolate data to 1 kyr timesteps

time_interp = 0:1e3:2.7e6; 
oxygen_interp = interp1(time, oxygen_raw, time_interp);


%% filter

% smooth data using resolution specified in input
b = (1/smoothing_time)*ones(1, smoothing_time);
a = 1; 

oxygen_isotopes = filter(b, a, oxygen_interp);

oxygen_isotopes(1:2*smoothing_time) = oxygen_interp(1:2*smoothing_time);

%% get all in correct shape and direction
oxygen_isotopes = oxygen_isotopes';
time_interp = time_interp';

% start with earliest time
oxygen_isotopes = flipud(oxygen_isotopes);
time = flipud(time_interp);

%% create mask usign glaciation threshold

glaciation_mask = oxygen_isotopes >= glaciation_threshold; % 1 = ice-covered, 0 = ice-free

% find index where each ice cover period starts and ends

[switches, ~] = find(or(diff(glaciation_mask) == 1, diff(glaciation_mask) == -1));

% if ice-covered at the beginning of the pleistocene, as set by the d18O
% threshold, then set first index as a switch into a glacial period so that 
% model begins at 2.7 Myr

if oxygen_isotopes(1) > glaciation_threshold
    switches = [1; switches]; 
end

%% find times when ice-cover begins or ends

switch_times = time(switches);
switch_times(end) = history.deglac_t;


%% determine model time

% find duration of each time step. odd indices are glacial periods; 
% even indices are interglacial periods. 

model_times = abs(diff(switch_times));

% pull out vector with durations of glacials and interglacials

glacial_lengths = model_times(1:2:end);
interglacial_lengths = model_times(2:2:end);

end
