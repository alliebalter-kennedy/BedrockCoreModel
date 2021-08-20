function [N_final, depths_final] = LIA_CoreModel(loc, history, erosion, steady_state)

%% Load Constants

consts = bedrock_constants();

% should update this so it is a stand-alone script. 

%% Inputs to calculate production 

dz = erosion.dz; 

rho = 2.65;            % [g cm^-3]; rock density
z_cm = [0:dz:4000000]; % [cm]; length of starting production profile.
                       % needs to be longer than total erosion over model.
z_gcm2 = z_cm.*rho;    % [g cm^-2]; depth


profile_vec = erosion.profile_vec'.*rho; % [g cm^-2];

% Atmospheric pressure at site
site_p = ERA40atm(loc.lat,loc.long,loc.elv); % site air pressure


% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later after Balco, 2017. 
m = build_muon_profile_w14c(site_p,consts,0);

% Define production rate info
SFsp = stone2000(loc.lat,site_p,1); % scaling factor

% Build a data structure with production rate information
p.P10sp = consts.P10q_St .* SFsp; % Be-10 spallation production rate at surface
p.P26sp = p.P10sp.*consts.R2610q; % Al-26 spallation production rate at surface
p.P14sp = consts.P14q_St.*SFsp; % C-14 spallation production rate at surface

% Attenuation
p.Lsp = 160; % g/cm2.

% Define total production

P10z = PofZ(z_gcm2, m, p, 10); % sum of production by spallation and muons; Balco (2017)

%% Unwrap variables 

deglac_t = history.deglac_t; % local deglaciation age
historical_cover = history.historical_cover; % duration of historical cover
historical_deglac = history.historical_deglac; % known recent deglaciation from satellite imagery.
model_times = history.model_times'; % pre-deglaciation exposure history as determined by the d18O threshold
length_gl = history.glacial_lengths; % length of each ice-covered period
length_int = history.interglacial_lengths; % length of each ice-free period
 
%% Model Set Up

% Holocene History information

EarlyHolo_t = deglac_t - historical_cover - historical_deglac; % Exposure duration from deglaciation to the start of historical cover

%% Create gl_int vectors with d18O data

% add Holocene history information
timesteps = [model_times EarlyHolo_t historical_cover historical_deglac]; % concatenate all exposure history information

starttime = max(fliplr(cumsum(timesteps))); % define model start time

N_final = zeros(length(profile_vec), length(erosion.gl_erosion), length(erosion.historical_erosion)); % create matrix for final nuclide concentrations

for b = 1:length(erosion.gl_erosion) % loop over orbital-scale erosion rates
    for c = 1:length(erosion.historical_erosion) % loop over historical erosion rates

    % define erosion

    gl_erosion_rate = (erosion.gl_erosion(b)./10) .* rho; % erosion rate for each period of ice cover [cm/yr]
    gl_erosion = gl_erosion_rate .* length_gl; % total erosion depth for each period of ice cover [cm]

    historical_erosion_rate = (erosion.historical_erosion(c)./10) .* rho; % historical erosion rates [cm/yr]
    historical_erosion = historical_erosion_rate .* historical_cover; % total erosion depth for historical ice cover [cm]

    addint = 0; % force model to start during a period of ice cover to start eroding steady state profile. inititlize binary variable. 
    if length(length_int) < length(length_gl) % if fewer interglacials than glacials
        addint = 1; 
        length_int = [0; length_int]; % add interglcial of length zero at first interglacial. 
    end

    gl_erosion = gl_erosion.*ones(length(length_gl), size(gl_erosion, 2)); % create vector with erosion depths during ice cover
    int_erosion = 0.*ones(length(length_int), 1); % create vector with erosion depths during ice free periods. right now just 0, but could add e/L to nuclide calcs.        

    % interleave vectors for erosion during glacials and interglacials.

    erode = zeros(size(gl_erosion, 2), size(length_int, 1).*2);
    erode = (reshape([int_erosion(:) gl_erosion(:)]', 2*size(length_int, 1), [])');


    if addint == 1;
        length_int = length_int(2:end); % if needed to add interglacial at the beginning of model time, select erosion starting with first glacial period. 
        int_erosion = int_erosion(2:end);
        erode = erode(2:end);
    else
    end


    % add Holocene erosion information to "erode" vector
    erode = [erode 0 historical_erosion 0];

    erode = round(erode, 2); % round to nearest 0.1 cm.
    startdepth = max(fliplr(cumsum(erode))); % find starting depth of modern rock surface

    start_index = ceil(((startdepth./rho)./dz)+1); % find index of the starting depth

    %% Run Model

    % intitialize

    topdepth_new = startdepth; 
    topindex_new = start_index;
    bottomindex_new = ceil(start_index + (length(profile_vec)-1));

    if steady_state.flag == 0 % use pre-made steady state profile if indicated by flag. this will be starting nuclide concentrations. 
        N_old = zeros(length(profile_vec), 1);
    elseif steady_state.flag == 1
        N_old = steady_state.N_start(topindex_new:bottomindex_new);
    else  display('WARNING: Pick steady state flag 1 or 0');
    end

    topdepth_old = topdepth_new;
    topindex_old = topindex_new;

    % loop through model time

    for a = 1:length(timesteps)
    
        % expose for even a (interglacials)
    
        if mod(a, 2) == 0
            % calculate cosmogenic nuclide accumulation during exposure
            N_new = N_old.*exp(-consts.l10.*timesteps(a)) + P10z(topindex_new:bottomindex_new)'./consts.l10 .* (1-exp(-consts.l10.*timesteps(a)));
            % no erosion takes place
            topindex_new = topindex_old;
        else
        % bury and erode for odd a (glacials)
            % cosmogenic nucldie decay
            N_new = N_old.*exp(-consts.l10.*timesteps(a));
            
            % shift depth profile up so that next exposure period begins at
            % new production rate
            topdepth_new = round(topdepth_old - erode(a), 2); % depth at the end of the glacial; round to nearest 0.1 cm
            topindex_new = ceil(((topdepth_new./rho)./dz)+1); % find index for top of depth profile. will be used to find production rate during next exposure period. 
            bottomindex_new = ceil(topindex_new + (length(profile_vec)-1)); % find index for bottom of depth profile.
        end
    
        topdepth_old = topdepth_new;
        topindex_old = topindex_new;
        N_old = N_new;
   
    end
    
    % save results
    
    N_final(:, b, c) = N_new; % save final depth profiles at 0.1 cm spacing for each combo of erosion rates.
    depths_final(:, b, c) = z_gcm2(topindex_new:bottomindex_new)'; % all final depths should be zero. this passes through test in model run script.


end 
end

out = [N_final, depths_final];
end

