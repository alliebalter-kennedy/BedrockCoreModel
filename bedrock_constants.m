%% Constants for bedrock profiles 

% Code in this scipt is adapted from code written by
% Greg Balco, script titled GISP2_data_analysis.m (lines 124-191), which
% accompanies Schaefer et al. (2016), andcan be found at 
% http://hess.ess.washington.edu/repository/GISP2/. 

function consts = bedrock_constants()

consts.P10q_St = 4.043; % arctic PR; Young et al., 2013
consts.delP10q_St = 0.153; % Corresponding uncert.

consts.P14q_St = 13.590; % Svalbard Saturated Peaks 
consts.delP14q_St = 0.443; % Svalbard Saturated Peaks

consts.R2610q = 6.75;

% Muon interaction cross-sections

% Be-10 from Balco et al., 2017
consts.mc10q.Natoms = 2.006e22;
consts.mc10q.k_neg = 0.00191 .* 0.704 .* 0.1828; % Balco, 2017
consts.mc10q.sigma0 = 0.280e-30; % Balco, 2017

% Al-26 from Balco et al., 2017
consts.mc26q.Natoms = 1.003e22;
consts.mc26q.k_neg = 0.0133 .* 0.296 .* 0.6559; % from Balco, 2017
consts.mc26q.sigma0 = 3.89e-30; % Balco, 2017

% C-14 from Balco et al., 2017
consts.mc14q.Natoms = 2.006e22;
consts.mc14q.k_neg = 0.116 .* 0.704 .* 0.1828; % from Balco, 2017
consts.mc14q.sigma0 = 0.45e-27./190; % Heisinger et al., 2002

% Also define decay constants

consts.l10 = -log(0.5)./1.387e6; % Chmeleff/Korschinek value
consts.l26 = 9.83e-7; % Compatible with Nishiizumi (2004)
consts.l14 = -log(0.5)./5730; 

consts = consts;
end