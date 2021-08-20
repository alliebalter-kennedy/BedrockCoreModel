function out = build_muon_profile_w14c(pressure,consts,plotFlag)

% This code builds muon production profile for SLHL for interpretation of 
% GISP2 bedrock data. Basically, this is building a lookup table to
% facilitate numerical integration later.
%
% Greg Balco -- Berkeley Geochronology Center -- May 2016 
% 
% Syntax: out = build_muon_profile(pressure,consts,plotFlag)
% pressure is site air pressure in hPa
% consts is consts struct with muon params
% plotFlag is 0 or 1 to disable/enable check plots

% published in Schaefer et al. (2016)

% Updated by Allie Balter-Kennedy: January 2020
% Added constants for 14C production from Balco (2017)

% Schaefer, J. M., Finkel, R. C., Balco, G., Alley, R. B., Caffee, M. W., 
% Briner, J. P., et al. (2016). Greenland was nearly ice-free for extended 
% periods during the Pleistocene. Nature, 540(7632), 252?255. 
% https://doi.org/10.1038/nature20146

% Define lower boundary of precise calculation
zbot = 20000;

% Plotflag default
if nargin < 3; plotFlag = 0; end;

% Compute muon production rates on log mesh
m.zz = logspace(0,log10(zbot),200);
for a = 1:length(m.zz);
    m.P10mu(a) = P_mu_total_alpha1(m.zz(a),pressure,consts.mc10q);
    m.P26mu(a) = P_mu_total_alpha1(m.zz(a),pressure,consts.mc26q);
    m.P14mu(a) = P_mu_total_alpha1(m.zz(a),pressure,consts.mc14q);
end;

% Add zero depth to the log mesh
m.zz = [0 m.zz];
m.P10mu = [P_mu_total_alpha1(0,pressure,consts.mc10q) m.P10mu];
m.P26mu = [P_mu_total_alpha1(0,pressure,consts.mc26q) m.P26mu];
m.P14mu = [P_mu_total_alpha1(0,pressure,consts.mc14q) m.P14mu];

% Plotting 

if plotFlag == 1;
    figure;
    loglog(m.zz,m.P10mu,'b',m.zz,m.P10mu,'bo');
    hold on; 
    plot(m.zz,m.P26mu,'g',m.zz,m.P26mu,'go');
    hold on; 
    plot(m.zz,m.P14mu,'g',m.zz,m.P14mu,'go');
end;


%% This part does does an exponential fit to bottom 2 data. Used for 
% extrapolating beyond lowermost calculated depth from above. 

L10 = 1/(-(1/(m.zz(end)-m.zz(end-1))).*log(m.P10mu(end)./m.P10mu(end-1)));
L26 = 1/(-(1/(m.zz(end)-m.zz(end-1))).*log(m.P26mu(end)./m.P26mu(end-1)));
L14 = 1/(-(1/(m.zz(end)-m.zz(end-1))).*log(m.P14mu(end)./m.P14mu(end-1)));

m.L10bot = L10;
m.L26bot = L26;
m.L14bot = L14;
% More plotting

if plotFlag == 1;
    figure;
    semilogy(m.zz((end-5):end),m.P10mu((end-5):end),'go');
    hold on; plot(m.zz((end-5):end),m.P26mu((end-5):end),'bo');
    hold on; plot(m.zz((end-5):end),m.P14mu((end-5):end),'ko');

    px = (m.zz(end):100:(m.zz(end)+10000));
    py10 = m.P10mu(end).*exp(-(px-m.zz(end))./L10);
    py26 = m.P26mu(end).*exp(-(px-m.zz(end))./L26);
    plot(px,py10,'g',px,py26,'b');

end;

% Save results

out = m;



