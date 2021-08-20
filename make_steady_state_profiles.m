clear all
close all

% 
%% Location information

% Latitude and elevation of site
loc.lat = 69.161154;
loc.long = -49.797913;
loc.elv = 93;
corename = '18JAK-CR1';


%% Load Constants

consts = bedrock_constants();

% note that if I want to change the production rate, I have to do so in
% that function. I'd like to come up with a better way to do this. 

%% Inputs to calculate production 

dz = 0.1;

rho = 2.65;            % [g cm^-3]; rock density
z_cm = [0:dz:40000]; % [cm]; length of starting production profile.
                       % needs to be longer than total erosion over model.
z_gcm2 = z_cm.*rho;    % [g cm^-2]; depth


% Atmospheric pressure at site
site_p = ERA40atm(loc.lat,loc.long,loc.elv); % site air pressure


% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later.
m = build_muon_profile_w14c(site_p,consts,0);

% Define production rate info
SFsp = stone2000(loc.lat,site_p,1); % scaling factor

% Build a data structure with production rate information
p.P10sp = consts.P10q_St .* SFsp; % Be-10 spallation production rate at surface 
p.P26sp = p.P10sp.*consts.R2610q; % Al-26 spallation production rate at surface
p.P14sp = consts.P14q_St.*SFsp; % C-14 spallation production rate at surface

% Attenuation
p.Lsp = 160; % g/cm2, appropriate for high latitude.


%% Define erosion rate

ee = 50; % m/Ma -- need to get to g/cm^2/yr
eee = (ee .* 100 .* rho) ./ 10^6; 

%% Integrate through time

tic
for i = 1:length(z_gcm2);

    pfun = @(t) PofZ((z_gcm2(i) + eee.*t), m, p, 10) .* exp(-consts.l10 .* t);

    N_start(i, 1) = integral(pfun, 0, 15e6);
end
toc

%% save
FormatFile = '/SteadyStateProfiles/SteadyStateErosionRate_%d.mat';

filename = sprintf(FormatFile, ee);

save(filename, 'N_start');

beep on