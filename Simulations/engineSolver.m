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
%   .coolant
%   [VaporState, Temperature, Pressure, Velocity]
%   .nozzle
%   [Temperature, Pressure, MachNumber]
%   .chamber
%   [HotWallTemperature, ColdWallTemperature, Stress]
%   .mass
%   [MdotO, MdotF]
%   .meta    - remember input file 
%   .meta = in

%% Stage 0 Setup

% Store input metadata
out.meta = in;

% Output station positions used for calculations
out.station = (0:in.numerics.simDx:(in.nozzle.ChamberLength + in.nozzle.NozzleLength));

% Generate geometry data for all other processes

x = out.station;

[StationRad, StationArea] = combustionChamberProfile(x, in.geom.ChamberRadius, ...
        in.geom.BlendRadius, in.geom.ConvergingAngle, in.geom.BlendRadius, in.geom.ThroatRadius, in.geom.ChamberLength, ...
        in.geom.DivergingAngle, in.geom.ExitRadius, in.geom.ExitAngle,in.geom.NozzleLength);


%% Stage 1 Sim
converged = false;
PC = 300 * 6894.76; % First guess chamber pressure
Pdrop = 100 * 6894.76;  % First guess pressure drop across channels
while ~converged
    % Calculate mass flow of oxidiser based on first guess chamber
    % pressure and external feed pressure
    % Calculate mass flow of fuel based on first guess chamber pressure,
    % first guess pressure drop over regen channels, and feed pressure
    % Calculate first o/f ratio 
    % Calculate thermo stuff based on o/f and PC
    % Calculate required mass flow from nozzle based on thermo properties
    % Compute mass flow error, use to scale new chamber pressure
    % Calculate mach number along nozzle length, temperature, and pressure,
    % including logic for subsonic/supersonic/shock wave in expanding
    % section solutions 
    % Calculate heating equilibrium for hot wall, cold wall, coolant
    % temperature, and coolant pressure along the channel. 
    % Establish new pressure drop 
    % Compare new pressure drop with old pressure drop guess for error
    % value 
    % Adjust pressure drop based on error values

    % Check for convergance where errors are both within bounds
end



end