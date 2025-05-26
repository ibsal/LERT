%% Engine Solver
% Isolate engine 




function out = engineSolver(in)
%% INPUT (struct 'in')
%   .nozzle      – geometry sub‑struct
%   [ChamberLength, ChamberRadius, ExitRadius, ThroatRadius, NozzleLength, BlendRadius, ConvergingAngle, DivergingAngle, ExitAngle]
%   .prop      – propellant & CEA sub‑struct
%   [TransportCEA]
%   .cool      – coolant & channel sub‑struct
%   [Ncc, WallThickness, Height, Width, FinWidth]
%   .material
%   [YoungsModulus, CTE, Poisson, SurfaceRoughness, Strength, Conductivity, MaxTemperature]
%   .environment    - Actual "inputs" into the engine
%   [AmbientPressure, FuelPressure, FuelTemperature, OxPressure]
%   .numerics   - Useful for simulation refinement 
%   [simDx, dt]
%   .injector     - Injector information
%   [OxOrifices, FuelOrifices, OxDiameter, FuelDiameter, OxCd, FuelCd]
%% OUTPUT (struct 'out')
%   .station   - x values at which values are computed 
%   [x, Rad, Area]
%   .coolant
%   [VaporState, Temperature, Pressure, Velocity]
%   .nozzle
%   [Temperature, Pressure, MachNumber]
%   .chamber
%   [HotWallTemperature, ColdWallTemperature, Stress]
%   .mass
%   [MdotO, MdotF, Mdot, OF]
%   .meta    - remember input file 
%   .meta = in

%% Stage 0 Setup

% Store input metadata
out.meta = in;

% Output station positions used for calculations
out.station.x = (0:in.numerics.simDx:(in.nozzle.ChamberLength + in.nozzle.NozzleLength));
%% NOTE: IT COULD BE SMART TO RESOLVE TO HAVE GREATER RESOLUTION NEAR THE THROAT? IDK???

% Generate geometry data for all other processes

x = out.station.x;

[StationRad, StationArea] = combustionChamberProfile(x, in.nozzle.ChamberRadius, ...
        in.nozzle.BlendRadius, in.nozzle.ConvergingAngle, in.nozzle.BlendRadius, in.nozzle.ThroatRadius, in.nozzle.ChamberLength, ...
        in.nozzle.DivergingAngle, in.nozzle.ExitRadius, in.nozzle.ExitAngle,in.nozzle.NozzleLength);

out.station.Rad = StationRad;
out.station.Area = StationArea;
OxArea = in.injector.OxDiameter^2 * 0.25 * pi * in.injector.OxOrifices;
FuelArea = in.injector.FuelDiameter^2 * 0.25 * pi * in.injector.FuelOrifices;

%% Stage 1 Sim
converged = 0;
PC = min(in.environment.OxPressure, in.environment.FuelPressure) * 0.5; % First guess chamber pressure
Pdrop = in.environment.FuelPressure * 0;  % First guess pressure drop across channels
Pb = in.environment.AmbientPressure;
while ~converged
    % Calculate oxidizer and fuel densities
    rhoO = 1220;
    rhoF = 786;  %% OVERIDE WITH ACTUAL LOOKUP TABLES PLS
    % Calculate mass flow of oxidiser based on first guess chamber
    % pressure and external feed pressure
    MdotO = sign(in.environment.OxPressure-PC) *OxArea * in.injector.OxCd * sqrt(2 * rhoO * abs(in.environment.OxPressure-PC));
    % Calculate mass flow of fuel based on first guess chamber pressure,
    % first guess pressure drop over regen channels, and feed pressure
    MdotF = sign(in.environment.FuelPressure - Pdrop - PC) * FuelArea * in.injector.FuelCd * sqrt(2 * rhoF * abs(in.environment.FuelPressure - Pdrop - PC));
    % Calculate first o/f ratio
    of = MdotO/MdotF;
    % Calculate thermo stuff based on o/f and PC  %% THIS IS MISSING PC
    % INPUT ON THE SOLVER!!!! NEED NEW CEA DATA UGHHHH
    Tcomb   = ceaGridInterp('rocketCEAgrid.csv','T_c', PC, of);
    GammaC = ceaGridInterp('rocketCEAgrid.csv','gamma_c',  PC, of);
    if GammaC < 1
        error("GAMMAC < 1")
    end
    Rhoc = ceaGridInterp('rocketCEAgrid.csv', 'rho_c', PC, of);
    Cp = ceaGridInterp('rocketCEAgrid.csv', 'Cp_c', PC, of);
    Cv = Cp/GammaC;
    Rgas = Cp - Cv;
    % Calculate required mass flow from nozzle based on thermo properties
    Pratio = in.environment.AmbientPressure/PC;
    PratioCritical = (2/(GammaC + 1))^(GammaC/(GammaC-1));
    if(PC<Pb)
        %Backflow
        error("BACKFLOW HAS OCCURED");
        Backflow = 1;
    elseif(Pratio < PratioCritical)
        % Choked flow
        RealMdot = PC * min(StationArea) * (Tcomb)^-0.5 * sqrt(GammaC/Rgas) * (0.5*(GammaC + 1))^((-1 - GammaC)/(2 * (GammaC -1)));
        Sonic = 1;
        Backflow = 0;
    else 
        %% THIS DOESN't WORK 
        func = @(M) (1 + 0.5*(GammaC - 1)*M^2)^((-1*GammaC)/(GammaC-1)) - Pb/PC;
        Mb = fzero(func, [1e-6, 0.999999999]); %Mach number at exit area 
        RealMdot = StationArea(end) * PC * sqrt(Tcomb) * sqrt(GammaC/Rgas) * Mb * (1 + 0.5*(GammaC-1) * Mb^2) ^((-1 - GammaC)/(2* (GammaC-1)));
        Sonic = 0;
        Backflow = 0;
    end
    % Compute mass flow error
    Merror = RealMdot - (MdotF + MdotO);
    out.mass.OF = of;
    out.mass.MdotO = MdotO;
    out.mass.MdotF = MdotF;
    out.mass.Mdot = RealMdot;
    [M, T, P] = nozzleMach(PC,Pb,GammaC, StationArea, Tcomb);
    out.nozzle.Temperature = T;
    out.nozzle.Pressure = P;
    out.nozzle.MachNumber = M;
    %% Pressure drop section 
    % Calculate heating equilibrium for hot wall, cold wall, coolant
    % temperature, and coolant pressure along the channel. 
    % Establish new pressure drop 
    % Compare new pressure drop with old pressure drop guess for error
    % value 
    % Adjust pressure drop based on error values
    % Adjust chamber pressure based on mass flow error
    PC = (PC - (PC * 0.1 * Merror));
    if abs(Merror) < 1e-3
        converged = 1;
    else
        converged = 0;
    end

    % Check for convergance where errors are both within bounds
end

if out.mass.MdotF<0 || out.mass.MdotO<0
    error("FLOW THROUGH INJECTOR")
end