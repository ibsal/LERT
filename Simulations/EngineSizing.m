%--------------------------------------------------------------
%  REGENERATIVE-COOLED NOZZLE 1-D THERMAL / STRESS MODEL
%  -- Stand-alone version  --
%--------------------------------------------------------------
%  * User-adjustable inputs are grouped in SECTION 1 below *
%--------------------------------------------------------------
%  * Code by Isaac Sal. Nice comment blocks by ChatGPT
%--------------------------------------------------------------
clear; clc;

%% ========= 1) USER-ADJUSTABLE PARAMETERS =============================
% -- GLOBAL DESIGN TARGETS ---------------------------------------------
Thrust         = 250 * 4.44822;        % [N]  design thrust
OFRatio        = 3;                  % [-]  oxidizer / fuel mass ratio

% -- COMBUSTION CONDITIONS ---------------------------------------------
CombustionPressure = 300 * 6894.76;    % [Pa]
AmbientPressure    = 14.6959 * 6894.76;% [Pa]

% -- CHAMBER / NOZZLE GEOMETRY (inches, gets converted internally) -----
ChamberLengthIN = 3;                   % chamber length  [in]
ChamberRadiusIN = 1;                   % chamber radius  [in]
NozzleLengthIN  = 1.5;                 % nozzle length   [in]
filletsIN       = 1;                   % blend radius    [in]

% -- COOLANT-CHANNEL & WALL GEOMETRY -----------------------------------
Ncc            = 40;       % number of coolant channels
wallthickness  = 0.03;     % [in] wall separating hot & cold sides
height         = 0.03;     % [in] channel height  (open gap)
width          = 0.05;     % [in] channel width
finwidth       = 0.04;     % [in] fin (land) width between channels

% -- THERMAL / MATERIAL PROPERTIES -------------------------------------
hcoolant  = 25e3;          % [W/m^2-K] coolant-side h-coefficient
k          = 15;           % [W/m-K]   wall conductivity
Tcoolantinit = 300;        % [K]       coolant inlet temperature
Pcoolantinit = 500 * 6894.76; % [Pa]   coolant inlet pressure
cpIPA      = 2600;         % [J/kg-K]  IPA specific heat (assumed const.)

Echamber   = 180e9;        % [Pa] Young's modulus of wall
CTEchamber = 17e-6;        % [1/K] thermal expansion coeff.
Poissons   = 0.27;         % [-]
Strength   = 400e6;        % [Pa] yield (for margin calc)
roughness  = 6.3e-6;       % [m]

StressFOS  = 2;            % [-]
TempFOS    = 1.2;          % [-]

RecoveryFactor = 0.91;     % 1 is worst case: assumes stagnation heating

% -- NUMERICAL SETTINGS -------------------------------------------------
simDx = 0.0025;            % [in] axial marching step
Tmax  = 1100;              % [K]   hot-wall over-temp shading limit
sfos  = StressFOS;                 % [-]   stress factor-of-safety requirement
tfos  = TempFOS;               % [-]   temp  factor-of-safety requirement

%% ========= 2) DERIVED CONSTANTS & PRECALC ============================
IN2M   = 0.0254;            % in → m
PSI2PA = 6894.76;           % psi → Pa  (for convenience later)

offset = 0;  % dummy for plotting coolant at x(end)+dx (keep orig logic)

%% ========= 3) CEA LOOK-UPS & INITIAL ENGINE GEOMETRY =================
cea     = readCEA('rocketN2OxIPA.txt');
Cstar   = ceaInterp(cea,'Cstar', 'O_F', OFRatio);
Tcomb   = ceaInterp(cea,'T_c',   'O_F', OFRatio);
Aexit   = ceaInterp(cea,'a_e',   'O_F', OFRatio);
Mexit   = ceaInterp(cea,'Mach_e','O_F', OFRatio);
CEAreaRatio = ceaInterp(cea,'Ae_At','O_F', OFRatio);
Vexit   = Mexit * Aexit;

Mflow = Thrust / Vexit;     % total propellant mass-flow [kg/s]
mdotF = Mflow / (1 + OFRatio);      % fuel  mdot [kg/s]
mdotO = Mflow - mdotF;              % oxid. mdot [kg/s]

ThroatArea   = Mflow * Cstar / CombustionPressure;  % [m^2]
ThroatRadius = sqrt(ThroatArea/pi);                 % [m]
ExitArea     = CEAreaRatio * ThroatArea;            % [m^2]

ThroatRadiusIN = ThroatRadius / IN2M;               % [in]
ExitRadiusIN   = sqrt(ExitArea/pi) / IN2M;          % [in]

Lstar = (ChamberLengthIN * ChamberRadiusIN^2 * pi) / (ThroatRadiusIN^2 * pi);

%% ========= 4) AXIAL GRID & COMBUSTION-SIDE PROFILES ==================
% axial stations from nozzle exit (x=0) → chamber rear-wall
x = 0:simDx:(ChamberLengthIN + NozzleLengthIN);
[stationrad, stationarea] = combustionChamberProfile(x, ChamberRadiusIN, ...
        filletsIN, 45, filletsIN, ThroatRadiusIN, ChamberLengthIN, ...
        15, ExitRadiusIN, 3, NozzleLengthIN);

gamma_e    = ceaInterp(cea,'gamma_e','O_F',OFRatio);
Machcurve  = machFromArea(stationrad, gamma_e, ThroatRadiusIN); % supersonic
T          = Tcomb ./ (1 + (gamma_e-1)/2 .* Machcurve.^2);
P          = CombustionPressure .* (1 + (gamma_e-1).*0.5 .* Machcurve.^2) ...
                                          .^(-gamma_e./(gamma_e-1));
Taw        = T .* (1 + RecoveryFactor*0.5*(gamma_e - 1) .* Machcurve.^2);                          % simple 0.9 factor

%% ========= 5) COOLANT & WALL SOLUTION LOOP ===========================
N = numel(stationrad);        % number of axial cells

% ----- Pre-allocate result arrays -------------------------------------
Tcw  = zeros(1,N);    Thw  = Tcw;    q = Tcw;    he = Tcw;
Tbulkcoolant      = zeros(1,N+1);    Tbulkcoolant(end) = Tcoolantinit;
Pcoolant          = zeros(1,N+1);    Pcoolant(end)     = Pcoolantinit;
Rhocoolant        = zeros(1,N+1);
Vcoolant          = zeros(1,N+1);
IPAkinematicVisc  = zeros(1,N+1);
ReIPA             = zeros(1,N+1);
FrictionFactor    = zeros(1,N+1);
Vapor             = false(1,N);
Stress            = zeros(1,N);

% ---- Initial coolant properties at nozzle exit (node N+1) ------------
Ripa                = 8.3144598/60.0950 * 1e3;            % [J/kg-K]
Rhocoolant(end)     = Pcoolant(end)/(Ripa*Tbulkcoolant(end));
HydraulicDiameter_m = 4*(width*IN2M*height*IN2M)/ ...
                       (2*IN2M*width + 2*IN2M*height);
Vcoolant(end)       = (mdotF/Ncc) / (Rhocoolant(end) * width*IN2M * height*IN2M);
IPAkinematicVisc(end) = ipaKinV(Tbulkcoolant(end));
ReIPA(end)          = Vcoolant(end) * HydraulicDiameter_m / IPAkinematicVisc(end);
FrictionFactor(end) = darcy(roughness,HydraulicDiameter_m,ReIPA(end));

% ---- Helper area vectors ---------------------------------------------
Agas = 2*pi*stationrad*IN2M .* (simDx*IN2M) / Ncc;    % per channel [m^2]
Lc   = height*IN2M + 0.5*finwidth*IN2M/2;             % fin length [m]
mFin = sqrt(2*hcoolant/(k*0.5*finwidth*IN2M));
nfn  = tanh(mFin*Lc)/(mFin*Lc);                       % fin efficiency
Acoolant = (2*nfn*height*IN2M + width*IN2M) * simDx*IN2M; % [m^2/cell]

% ---- Anon funcs -------------------------------------------------------
combustioncea = parseCEAtransport("2transportN2OxIPA.txt",'true');

h  = @(Thw_,i) bartz(P(i),T(i),Machcurve(i),stationarea(i)*0.00064516, ...
                     ThroatArea,filletsIN*IN2M,OFRatio,Thw_,Cstar,combustioncea);
qCoolant = @(Tcw_,Tbulk) hcoolant*Acoolant.*(Tcw_-Tbulk);
qWall    = @(Thw_,Tcw_,i) k*Agas(i).*(Thw_-Tcw_)/(wallthickness*IN2M);
qGas     = @(Thw_,i) h(Thw_,i).*Agas(i).*(Taw(i)-Thw_);
optsFS   = optimset('Display','off');

% ---- March from NOZZLE exit → chamber (reverse x) --------------------
for i = N:-1:1
    bulkIn = Tbulkcoolant(i+1);
    pin    = Pcoolant(i+1);

    f = @(x) [ qCoolant(x(1),bulkIn) - qGas(x(2),i) ; ...
               qCoolant(x(1),bulkIn) - qWall(x(2),x(1),i) ];
    var = fsolve(f,[bulkIn; T(i)],optsFS);   % x(1)=Tcw , x(2)=Thw

    Tcw(i) = var(1);   Thw(i) = var(2);   he(i) = h(Thw(i),i);
    q(i)   = qGas(Thw(i),i);

    % --- coolant property update --------------------------------------
    dT             = q(i)/( (mdotF/Ncc)*cpIPA );
    Tbulkcoolant(i)= bulkIn + dT;

    dP = (FrictionFactor(i+1)*Rhocoolant(i+1)*Vcoolant(i+1)^2* ...
           simDx*IN2M)/(2*HydraulicDiameter_m);
    Pcoolant(i) = pin - dP;

    Rhocoolant(i) = Pcoolant(i)/(Ripa*Tbulkcoolant(i));
    Vcoolant(i)   = (mdotF/Ncc)/(Rhocoolant(i)*width*IN2M*height*IN2M);
    IPAkinematicVisc(i)= ipaKinV(Tbulkcoolant(i));
    ReIPA(i)      = Vcoolant(i)*HydraulicDiameter_m/IPAkinematicVisc(i);
    FrictionFactor(i)= darcy(roughness,HydraulicDiameter_m,ReIPA(i));

    Vapor(i) = strcmpi(ipa(Pcoolant(i),Tbulkcoolant(i)),"VAPOR");

    % --- simple thin-wall hoop+thermal stress -------------------------
    Stress(i) = ((Pcoolant(i)-P(i))*stationrad(i)*IN2M)/(wallthickness*IN2M) + ...
                 (Echamber*CTEchamber*q(end)*wallthickness*IN2M)/ ...
                 (2*(1-Poissons)*k);
end

%% ========= 6) PLOTS ==================================================
figure('Name','Regen Nozzle Results','Color','w');

% 6a) Wall / coolant / gas temperatures
axT = subplot(2,2,1); hold(axT,'on'); grid(axT,'on');
plot(axT,x,Thw,'DisplayName','Hot Wall');
plot(axT,x,Tcw,'DisplayName','Cold Wall');
plot(axT,[x, x(end)+simDx],Tbulkcoolant,'DisplayName','Coolant');
plot(axT,x,T,'DisplayName','Combustion');
plot(axT,x,Taw,'DisplayName','Adiabatic');

% --- shading vapour and over-temp regions -----------------------------
yl = ylim(axT);
shade = @(idx,color,alpha) patch('XData',[x(idx(1)) x(idx(end)) x(idx(end)) x(idx(1))], ...
          'YData',[yl(1) yl(1) yl(2) yl(2)], ...
          'FaceColor',color,'FaceAlpha',alpha,'EdgeColor','none', ...
          'Parent',axT);

if any(Vapor),   shade(find(Vapor),[1 0.6 0.6],0.30); end
if any(Thw>Tmax),shade(find(Thw>Tmax),[1 0.8 0.4],0.25);end

legend(axT,'Location','best');
xlabel(axT,'Nozzle Station, in'); ylabel(axT,'Temperature [K]');

% 6b) Stress & temperature margins
subplot(2,2,2);
plot(x,Strength./Stress - sfos,'DisplayName','Stress Margin'); hold on; grid on;
plot(x,Tmax./Thw - tfos,'DisplayName','Temp  Margin');
legend; xlabel('Nozzle Station, in'); ylabel('Margin (-)');

% 6c) Coolant pressure profile
subplot(2,2,3);
plot([x, x(end)+simDx],Pcoolant*0.000145038); grid on;
xlabel('Nozzle Station, in'); ylabel('Coolant Pressure [psi]');

% 6d) Nozzle profile
subplot(2,2,4);
plot(x,stationrad,'k',x,-stationrad,'k'); axis equal;
xlabel('Axial Station, in'); ylabel('Radius, in'); title('Nozzle Profile');

%% ========= 7) SUMMARY OUTPUT =========================================
summary = struct( ...
    'MassFlow_kg_s',Mflow, ...
    'mdot_fuel',mdotF, 'mdot_oxidizer',mdotO, ...
    'ThroatRadius_in',ThroatRadiusIN, 'ExitRadius_in',ExitRadiusIN, ...
    'Lstar_in',Lstar, ...
    'MaxHotWallTemp_K',max(Thw), 'MaxStress_Pa',max(Stress), ...
    'MinStressMargin',min(Strength./Stress) );

fprintf('\n===== RUN SUMMARY =====\n'); disp(summary);
% Uncomment next line to save for post-processing:
% save('regen_summary.mat','summary');
