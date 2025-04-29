CombustionPressure = 300 * 6894.76; %PA
AmbientPressure = 14.6959 * 6894.76; %PA Standard atmosphere
CombustionGamma = 1.2;

ExitMach = sqrt( 2/(CombustionGamma-1) * ( (CombustionPressure./AmbientPressure).^((CombustionGamma-1)/CombustionGamma) - 1 ) );
AreaRatio = (1./ExitMach) .* ((2./(CombustionGamma+1)) .* (1 + (CombustionGamma-1)./2 .* ExitMach.^2)) .^ ((CombustionGamma+1)./(2*(CombustionGamma-1)));

Thrust = 250 * 4.44822; %Thrust, Newtons
OFRatio = 5;

cea = readCEA('ceaoutput.txt');
Mexit = ceaInterp(cea, 'Mach_e', 'O_F', OFRatio);
Aexit = ceaInterp(cea, 'a_e', 'O_F', OFRatio);
Cstar = ceaInterp(cea, 'Cstar', 'O_F', OFRatio);
Vexit = Mexit * Aexit;
CEAreaRatio= ceaInterp(cea, 'Ae_At', 'O_F', OFRatio);
Mflow = Thrust/Vexit; %kg/s

ThroatArea = Mflow * Cstar / CombustionPressure; %m^2
ThroatRadius = sqrt(ThroatArea/pi); %m
ExitArea = CEAreaRatio * (ThroatArea);
ThroatRadiusIN = ThroatRadius * 39.3701;
ExitRadius = sqrt(ExitArea/pi);
ExitRadiusIN = ExitRadius * 39.3701;

NozzleLengthIN = ExitRadiusIN * 1.8 * 2;

ChamberLengthIN = 3;
ChamberRadiusIN = .65;
Lstar = (ChamberLengthIN * ChamberRadiusIN^2 * pi)/(ThroatRadiusIN^2 * pi);

x = 0:0.01:(ChamberLengthIN + NozzleLengthIN);
[stationrad, stationarea] = combustionChamberProfile(x, ChamberRadiusIN, 0.25, 45, 0.25, ThroatRadiusIN, ChamberLengthIN, 15, ExitRadiusIN,2, NozzleLengthIN);

Tcombustion = ceaInterp(cea, 'T_c', 'O_F', OFRatio);
gamma_e = ceaInterp(cea, 'gamma_e', 'O_F', OFRatio);
Machcurve = machFromArea(stationrad, gamma_e, ThroatRadiusIN);          % supersonic

T    = Tcombustion ./ (1 + (gamma_e-1)/2 .* Machcurve.^2);

% Target mass flow for oxidizer and fuel
mdotF = Mflow/(1+OFRatio) 
mdotO = Mflow - mdotF

%% --- Bartz convective heat-flux -----------------------------------------
% throat geometry in metres
Dt   = 2*ThroatRadius;              % throat diameter
rc   = ThroatRadius;                % blend radius ≈ 1·rt

% bulk-gas properties at this O/F (CEA exit columns are close enough)
cp0  = ceaInterp(cea,'Cp','O_F',OFRatio);      % J/(kg·K)
Pr0  = ceaInterp(cea,'Pr','O_F',OFRatio);      % –
T0   = Tcombustion;                            % K (stagnation)
gamma= gamma_e;                                % exit γ

% simple Sutherland μ(T) [Pa·s]  (≈10 % error for IPA/N2O products)
muFcn = @(T) 1.458e-6 .* T.^1.5 ./ (T + 110.4);
mu0   = muFcn(T0);
Twall = 700;                                   % assume Cu wall limit
mu_w  = muFcn(Twall);

% heat-flux vector (W m⁻²) along x
qconv = bartzHeatFlux(Mflow, ThroatArea, Dt, rc, ...
                      cp0, T0, gamma, Pr0, mu0, mu_w, ...
                      Machcurve, Twall);


Gt = Mflow / ThroatArea;          % kg m^-2 s^-1   ≈ 1.3e3 ?
qEst = 0.026*cp0*Gt^0.8*T0^0.73  ... % plug values
        /(Dt^0.2*rc^0.1*Pr0^0.6)...
        *(mu0/mu_w)^0.2*(Twall/T0)^0.68;