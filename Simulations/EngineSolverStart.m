addpath 'Interpolation WIP'

IN2M = 0.0254;
PSI2PA = 6894.76;
%% ─────────────────────────── Engine Input Struct ───────────────────────────
% All units are SI unless otherwise noted.  Adjust any placeholder values
% (marked with “…”) to match your design or test‑case data.

in = struct();  % root struct container
% ────────────────────────────────────────────────────────────────────────────
% 1) Nozzle geometry
% ────────────────────────────────────────────────────────────────────────────
in.nozzle.ChamberLength   = 3*IN2M ;   % [m]
in.nozzle.ChamberRadius   = 1*IN2M ;   % [m]
in.nozzle.ExitRadius      = 2*IN2M ;   % [m]
in.nozzle.ThroatRadius    = 0.436*IN2M ;   % [m]
in.nozzle.NozzleLength    = 1.5*IN2M ;   % [m]
in.nozzle.BlendRadius     = 0.5*IN2M ;   % [m]  throat/contour fillet
in.nozzle.ConvergingAngle = 45 ;   % [deg]
in.nozzle.DivergingAngle  = 15 ;   % [deg]
in.nozzle.ExitAngle       = 2 ;   % [deg] half‑angle at nozzle exit

% ────────────────────────────────────────────────────────────────────────────
% 2) Propellant / CEA data
% ────────────────────────────────────────────────────────────────────────────
% Path or identifier for the CEA transport‑property file (set once—rest of
% the engine picks up all γ, Cp, μ, k, etc. from this reference):
in.prop.TransportCEA      = parseCEAtransport("2transportN2OxIPA.txt",'true');

% ────────────────────────────────────────────────────────────────────────────
% 3) Coolant‑channel geometry
% ────────────────────────────────────────────────────────────────────────────
in.cool.Ncc          = 20 ;   % [‑]   number of cooling channels
in.cool.WallThickness= 0.02*IN2M ;   % [m]   hot‑wall thickness
in.cool.Height       = 0.02*IN2M ;   % [m]   channel height
in.cool.Width        = 0.02*IN2M ;   % [m]   channel width
in.cool.FinWidth     = 0.02*IN2M ;   % [m]   fin (land) thickness between channels

% ────────────────────────────────────────────────────────────────────────────
% 4) Material properties (hot‑wall alloy or composite)
% ────────────────────────────────────────────────────────────────────────────
in.material.YoungsModulus    = 180e9 ;   % [Pa]
in.material.CTE              = 17e-6 ;   % [1/K]  coefficient of thermal expansion
in.material.Poisson          = 0.27 ;   % [‑]
in.material.SurfaceRoughness = 6.3e-6 ;   % [m]    arithmetic mean roughness (Ra)
in.material.Strength         = 400e6 ;   % [Pa]   ultimate or yield strength
in.material.Conductivity     = 15 ;   % [W/(m·K)]
in.material.MaxTemperature   = 1100 ;   % [K]    allowable wall temp

% ────────────────────────────────────────────────────────────────────────────
% 5) External/operating conditions
% ────────────────────────────────────────────────────────────────────────────
in.environment.AmbientPressure = 10135;   % [Pa]
in.environment.FuelPressure    = 900 *PSI2PA;   % [Pa]  injector inlet
in.environment.FuelTemperature = 300;   % [K]
in.environment.OxPressure      = 1000*PSI2PA ;   % [Pa]

% ────────────────────────────────────────────────────────────────────────────
% 6) Numerical‑grid controls
% ────────────────────────────────────────────────────────────────────────────
in.numerics.simDx = 0.01*IN2M ;   % [m] axial cell size
in.numerics.dt    = 1e-5 ;   % [s] time step

% ────────────────────────────────────────────────────────────────────────────
% 7) Injector parameters
% ────────────────────────────────────────────────────────────────────────────
in.injector.OxOrifices   = 12 ;     % [‑] count of oxidizer orifices
in.injector.FuelOrifices = 12 ;     % [‑] count of fuel  orifices
in.injector.OxDiameter   = 0.063*IN2M ;     % [m] orifice diameter
in.injector.FuelDiameter = 0.063*IN2M ;     % [m]
in.injector.OxCd         = 0.95 ;     % [‑] discharge coefficient
in.injector.FuelCd       = 0.95 ;     % [‑]

%% ─────────────────────── End of Input Definition ──────────────────────────
% Save or pass "in" to your solver/analysis script:
results = engineSolver(in);

plot(results.station.x, results.nozzle.Pressure)
hold on