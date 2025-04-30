function qdot = bartzHeatFlux(mdot, At, Dt, rc, ...
                               cp0, T0, gamma, Pr, mu0, mu_w, ...
                               Mvec, Twall)
% BARTZHEATFLUX   Compute convective heat-flux [W/m^2] along a nozzle
%
%   qdot = bartzHeatFlux(mdot, At, Dt, rc, cp0, T0, gamma, Pr, mu0, mu_w,
%                        Mvec, Twall)
%
%   mdot   – total mass-flow   [kg/s]
%   At     – throat area       [m^2]
%   Dt     – throat diameter   [m]
%   rc     – radius of curvature at throat (hot-gas side) [m]
%   cp0    – stagnation cp (CEA)  [J/(kg·K)]
%   T0     – stagnation temperature (CEA Tc) [K]
%   gamma  – ratio of specific heats at exit
%   Pr     – Prandtl number (use throat value or ~0.7-0.8)
%   mu0    – viscosity at T0   [Pa·s]
%   mu_w   – viscosity at wall Twall [Pa·s]
%   Mvec   – Mach number vector at station positions
%   Twall  – assumed wall temperature [K] (scalar)
%
%   Returns q'' as a vector the same length as Mvec.

    % throat mass flux
    Gt = mdot / At;                        % kg m^-2 s^-1

    % constant front factor
    C  = 0.026 * ...
          cp0 * Gt^0.8 * T0^0.73 / ...
          ( Dt^0.2 * rc^0.1 * Pr^0.6 ) * ...
          (mu0/mu_w)^0.2 * ...
          (Twall/T0)^0.68;

    % axial variation term
    theta = (1 + (gamma-1)/2 .* Mvec.^2).^(-0.34);  % –0.68/2 exponent

    qdot = C * theta;                     % W/m^2
end
