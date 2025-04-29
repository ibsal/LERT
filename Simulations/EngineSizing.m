

CombustionPressure = 300 * 6894.76; %PA
AmbientPressure = 14.6959 * 6894.76; %PA Standard atmosphere
CombustionGamma = 1.2;

ExitMach = sqrt( 2/(CombustionGamma-1) * ( (CombustionPressure./AmbientPressure).^((CombustionGamma-1)/CombustionGamma) - 1 ) );
AreaRatio = (1./ExitMach) .* ((2./(CombustionGamma+1)) .* (1 + (CombustionGamma-1)./2 .* ExitMach.^2)) .^ ((CombustionGamma+1)./(2*(CombustionGamma-1)));

Thrust = 250 * 4.44822; %Thrust, Newtons
OFRatio = 7;

cea = readCEA('ceaoutput.txt');
Mexit = ceaInterp(cea, 'Mach_e', 'O_F', OFRatio);
Aexit = ceaInterp(cea, 'a_e', 'O_F', OFRatio);
Cstar = ceaInterp(cea, 'Cstar', 'O_F', OFRatio);
Vexit = Mexit * Aexit;
CEAreaRatio= ceaInterp(cea, 'Ae_At', 'O_F', OFRatio);
Mflow = Thrust/Vexit; %kg/s

ThroatArea = Mflow * Cstar / CombustionPressure; %m^2
ThroatRadius = sqrt(ThroatArea)/pi; %m
ExitArea = CEAreaRatio * (ThroatArea);
ThroatRadiusIN = ThroatRadius * 39.3701
ExitRadius = sqrt(ExitArea)/pi;
ExitRadiusIN = ExitRadius * 39.3701

NozzleLengthIN = ExitRadiusIN * 1.8 * 2

ChamberLengthIN = 3;
ChamberRadiusIN = .65;
Lstar = (ChamberLengthIN * ChamberRadiusIN^2 * pi)/(ThroatRadiusIN^2 * pi)


x = 0:0.01:(ChamberLengthIN + NozzleLengthIN);
[stationrad, stationarea] = combustionChamberProfile(x, ChamberRadiusIN, 0.25, 45, 0.25, ThroatRadiusIN, ChamberLengthIN, 15, ExitRadiusIN,2, NozzleLengthIN);

Tcombustion = ceaInterp(cea, 'T_c', 'O_F', OFRatio);
gamma_e = ceaInterp(cea, 'gamma_e', 'O_F', OFRatio);
Machcurve = machFromArea(stationrad, gamma_e, ThroatRadiusIN);          % supersonic

T    = Tcombustion ./ (1 + (gamma_e-1)/2 .* Machcurve.^2);