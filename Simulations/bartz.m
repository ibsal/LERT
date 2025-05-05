function h = bartz(Pstation, Tstation, M, Astation, Athroat, Tfillet, of, Thw, cstar, tbl)
%% Bartz Heat Transfer Coeffecient Function
% Pstation = pressure at station
% Tstation = temperature at station 
% M = mach number at station
% Astation = cross sectional area at station 
% Athroat = throat area
% Tfillet = "fillet" or radius of curvature at throat 
% of = o/f ratio
% Thw = Temperature of combustion chamber wall "hot wall". This value will
% cstar = cstar, from interpolation of rocket sims with o/f ratios
% be iterated, so be smart about how you use it.

% In final implementation, tbl will be an additional input to avoid parsing
% large txt files many times. 
hgas = zeros(length(Pstation),1);
for i = 1:length(Pstation)

    %% Reference Temperature (T*) Calculation   
% Used for calculating some transport properties in this calculation. 

    Tstar =  Tstation(i) .* (1 + 0.032.*M(i).^2 + 0.58.*(Thw./Tstation(i) - 1));
    
    gamma = interpCEAtransport(tbl, Pstation(i), Tstar, of, 'gamma');
    Tstag = Tstation(i) * (1 + 0.5 * (gamma-1) * M(i)^2);
    Pstag = Pstation(i) * (1+ 0.5 * (gamma-1) * M(i)^2)^(gamma/(gamma-1));

    %% Bartz Calculation
    % Uses frozen solution 
    % Use Tstar for all transport properties: Viscocity and conductivity
    % Use temperature at station for all other properties 
    
    Dt = 2 * sqrt(Athroat/pi); % Throat Diameter
    mu = interpCEAtransport(tbl, Pstation(i), Tstar, of, 'mu', 'fr', 'nearest');  %combustion gas viscociyt 
    Cp = interpCEAtransport(tbl, Pstation(i), Tstation(i), of, 'Cp', 'fr', 'nearest');
    lambda = interpCEAtransport(tbl, Pstation(i), Tstar, of, 'k', 'fr', 'nearest'); %combustion gas thermal conductivity 
    Pr = (Cp * mu)/lambda; %combustion gas prandtl number
    
    
    % actual bartz calc
    
    sigma = ((0.5 * (Thw/Tstag) * (1 + 0.5*(gamma-1) * M(i)^2) + 0.5)^0.68 * (1 + 0.5*(gamma-1) * M(i)^2)^0.12)^-1;

    hgas(i) = 0.026 * Dt^-2 * mu^0.2 * Cp * Pr^-0.6 * Pstag^0.8 * cstar^-0.8 * Dt^0.1 * Tfillet^-0.1 * Athroat^0.9 * Astation(i)^-0.9 * sigma;
end
h = hgas;
end