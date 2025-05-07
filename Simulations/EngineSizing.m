CombustionPressure = 300 * 6894.76; %PA
AmbientPressure = 14.6959 * 6894.76; %PA Standard atmosphere
CombustionGamma = 1.2;

ExitMach = sqrt( 2/(CombustionGamma-1) * ( (CombustionPressure./AmbientPressure).^((CombustionGamma-1)/CombustionGamma) - 1 ) );
AreaRatio = (1./ExitMach) .* ((2./(CombustionGamma+1)) .* (1 + (CombustionGamma-1)./2 .* ExitMach.^2)) .^ ((CombustionGamma+1)./(2*(CombustionGamma-1)));

Thrust = 250 * 4.44822; %Thrust, Newtons
OFRatio = 3.5;

cea = readCEA('rocketN2OxIPA.txt');
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

NozzleLengthIN = 1.5;

ChamberLengthIN = 3;
ChamberRadiusIN = 3;
Lstar = (ChamberLengthIN * ChamberRadiusIN^2 * pi)/(ThroatRadiusIN^2 * pi);

simDx = 0.0025;

x = 0:simDx:(ChamberLengthIN + NozzleLengthIN);
fillets = 1; %in
[stationrad, stationarea] = combustionChamberProfile(x, ChamberRadiusIN, fillets, 45, fillets, ThroatRadiusIN, ChamberLengthIN, 15, ExitRadiusIN,3, NozzleLengthIN);

Tcombustion = ceaInterp(cea, 'T_c', 'O_F', OFRatio);
gamma_e = ceaInterp(cea, 'gamma_e', 'O_F', OFRatio);
Machcurve = machFromArea(stationrad, gamma_e, ThroatRadiusIN);          % supersonic

T    = Tcombustion ./ (1 + (gamma_e-1)/2 .* Machcurve.^2);
P = CombustionPressure .* (1 + (gamma_e-1).*0.5 .* Machcurve.^2).^(-1*gamma_e./(gamma_e-1));

% Target mass flow for oxidizer and fuel
mdotF = Mflow/(1+OFRatio) 
mdotO = Mflow - mdotF


% bulk-gas properties at this O/F (CEA exit columns are close enough)
cp0  = ceaInterp(cea,'Cp','O_F',OFRatio);      % J/(kg·K)
Pr0  = ceaInterp(cea,'Pr','O_F',OFRatio);      % –
T0   = Tcombustion;                            % K (stagnation)

combustioncea = parseCEAtransport("2transportN2OxIPA.txt", 'true');

h = @(Thw, i) bartz(P(i), T(i), Machcurve(i), stationarea(i).*0.00064516, ThroatArea, fillets*0.0254, OFRatio, Thw, Cstar, combustioncea);

%hg(Taw - Twg) = q = (k/t) * (Twg - Twc) = hc * (Twc - Tco)


Ncc = 100; %number of coolant channels
Agas = (2*pi*stationrad.*0.0254.*simDx.*0.0254)/Ncc;

Taw = T.*0.9;


wallthickness = 0.03; %in
hcoolant = 25*10^3 %W/M^2K
height = 0.03; %in
width = 0.03; %in
k = 15;
finwidth = 0.04; %in
Lc = height*0.0254 + (0.5 * finwidth*0.0254)/2;
m = sqrt(2*(hcoolant)/(k * 0.5 * finwidth*0.0254));
nfn = tanh(m*Lc)/(m*Lc);

Acoolant = (2*nfn*height*0.0254 + width*0.0254)*simDx*0.0254;

qCoolant = @(Tcw, Tbulkcoolant) hcoolant.*Acoolant.*(Tcw - Tbulkcoolant);

qWall = @(Thw, Tcw, i) k.*Agas(i).*(Thw -Tcw)./(wallthickness.*0.0254);

qGas = @(Thw, i) h(Thw, i).*Agas(i).*(Taw(i)-Thw);

Tcoolantinit = 300; %Initial Coolant temperature, K
Pcoolantinit = 400 * 6894.76; % Pa

N   = length(stationrad);          % total axial cells

% --- pre‑allocate result vectors (faster & avoids confusing zeros) ------
Tcw          = zeros(1,N);         % cold‑wall temperature
Thw          = zeros(1,N);         % hot‑wall temperature
q            = zeros(1,N);         % gas‑side heat flux
he           = zeros(1,N);         % Bartz hg
deltaT       = zeros(1,N);         % coolant ΔT per cell

Tbulkcoolant = zeros(1,N+1);       % bulk coolant temp at cell nodes
Tbulkcoolant(end) = Tcoolantinit;  % *** inlet is NOW at the NOZZLE exit ***
Pcoolant = zeros(1,N+1);
Pcoolant(end) = Pcoolantinit;
opts = optimset('Display','off');
cpIPA = 2600;
Vapor = zeros(1, N);

Ripa = 8.3144598/60.0950 * 10^3; %J/kg K
Rhocoolant = zeros(1,N+1);
Rhocoolant(end) = Pcoolant(end)/(Ripa * Tbulkcoolant(end));
Vcoolant = zeros(1,N+1);
Vcoolant(end) = (mdotF/Ncc)/(Rhocoolant(end) * width * 0.0254 * height * 0.0254);

IPAkinematicViscocity = ipaKinV(Tbulkcoolant(end)); %m^2/s
HydraulicDiameter = 4*(width *0.0254 *  height * 0.0254)/ (2*0.0254*width + 2*0.0254*height);
ReIPA = zeros(1, N+1);
ReIPA(end) = (Vcoolant(end) * HydraulicDiameter) / IPAkinematicViscocity;

roughness = 6.3*10^-6;
FrictionFactor = zeros(1, N+1);
FrictionFactor(end) = darcy(roughness, HydraulicDiameter,ReIPA(end));

for i = N:-1:1                         % *** march from nozzle → chamber ***
    
    bulkIn = Tbulkcoolant(i+1);       % coolant entering this cell (downstream)
    pin = Pcoolant(i+1);
    % --- FSOLVE for wall / coolant interface temps ---------------------
    f  = @(x) [ ...
        qCoolant(x(1), bulkIn) - qGas(x(2), i) ; ...
        qCoolant(x(1), bulkIn) - qWall(x(2), x(1), i) ];

    x0 = [bulkIn ;  T(i)];            % initial guess
    var= fsolve(f, x0, opts);         % var(1)=Tcw , var(2)=Thw

    % --- store results --------------------------------------------------
    Tcw(i) = var(1);
    Thw(i) = var(2);
    he(i)  = h(var(2), i);

    q(i)       = qGas(Thw(i), i);
    deltaT(i)  = q(i) / ((mdotF/Ncc) * cpIPA);
    Tbulkcoolant(i) = bulkIn + deltaT(i);
    deltaP = (FrictionFactor(i+1) * Rhocoolant(i+1) * Vcoolant(i+1)^2 * simDx * 0.0254)/(2*HydraulicDiameter);
    Pcoolant(i) = pin - deltaP;
    Rhocoolant(i) = Pcoolant(i)/(Ripa * Tbulkcoolant(i));
    Vcoolant(i) = (mdotF/Ncc)/(Rhocoolant(i) * width * 0.0254 * height * 0.0254);
    IPAkinematicViscocity = ipaKinV(Tbulkcoolant(i));
    ReIPA(i) = (Vcoolant(i) * HydraulicDiameter) / IPAkinematicViscocity;
    FrictionFactor(i) = darcy(roughness, HydraulicDiameter,ReIPA(i));

    % coolant leaving this cell (upstream) becomes bulkIn for next i‑1

    if(ipa(Pcoolant(i), Tbulkcoolant(i)) == "VAPOR")
        %disp("WARNING VAPOR IN COOOLANT LINES")
        Vapor(i) = 1;
    end
end

%% ---------------- P L O T S ----------------

%% ---------- temperatures + dual‑condition shading + legend labels ----
Tmax = 1100;                      % wall‑temperature limit (K)

axT  = subplot(2,2,1);  hold(axT,'on');

hHot   = plot(axT,x,Thw, 'DisplayName','Hot Wall');
hCold  = plot(axT,x,Tcw, 'DisplayName','Cold Wall');

if numel(Tbulkcoolant)==numel(x)+1
    xCool = [x, x(end)+simDx];
    hCool = plot(axT,xCool,Tbulkcoolant,'DisplayName','Coolant');
else
    hCool = plot(axT,x,Tbulkcoolant,   'DisplayName','Coolant');
end

hComb = plot(axT,x,T,   'DisplayName','Combustion');
hAw   = plot(axT,x,Taw, 'DisplayName','Adiabatic');
grid(axT,'on');

yl = ylim(axT);

% --- vapour shading (red) ---------------------------------------------
Vapor = logical(Vapor(:).');
d = diff([0 Vapor 0]);  s = find(d== 1);  e = find(d==-1)-1;
for k = 1:numel(s)
    patch('XData',[x(s(k)) x(e(k)) x(e(k)) x(s(k))], ...
          'YData',[yl(1)  yl(1)  yl(2)  yl(2)], ...
          'FaceColor',[1 0.6 0.6],'FaceAlpha',0.30,'EdgeColor','none', ...
          'Parent',axT);
end

% --- over‑temperature shading (orange) --------------------------------
over = Thw > Tmax;
d2 = diff([0 over 0]);  s2 = find(d2== 1);  e2 = find(d2==-1)-1;
for k = 1:numel(s2)
    patch('XData',[x(s2(k)) x(s2(k)) x(e2(k)) x(e2(k))], ...
          'YData',[yl(1)   yl(2)   yl(2)   yl(1)], ...
          'FaceColor',[1 0.8 0.4],'FaceAlpha',0.25,'EdgeColor','none', ...
          'Parent',axT);
end

% ---- dummy patches for legend keys -----------------------------------
hLegVap  = patch(NaN,NaN,[1 0.6 0.6],'FaceAlpha',0.30,'EdgeColor','none',...
                 'DisplayName','Vapour region');
hLegOver = patch(NaN,NaN,[1 0.8 0.4],'FaceAlpha',0.25,'EdgeColor','none',...
                 'DisplayName',sprintf('T_{hw} > %g K',Tmax));

legend(axT,[hHot hCold hCool hComb hAw hLegVap hLegOver],'Location','best');
hold(axT,'off');


% -------- wall‑adiabatic ΔT ------------------------------------------
subplot(2,2,2);
plot(x,Thw'-Taw); grid on;

% -------- Bartz h_g ---------------------------------------------------
subplot(2,2,3);
plot(x,he*1e-3); grid on;

% -------- radius profile ---------------------------------------------
subplot(2,2,4);
plot(x,stationrad,'k',x,-stationrad,'k'); axis equal;

% -------- coolant pressure (separate fig) ----------------------------
figure;
if numel(Pcoolant)==numel(x)+1
    xP = [x, x(end)+simDx];
    plot(xP,Pcoolant*0.000145038);
else
    plot(x,Pcoolant*0.000145038);
end
grid on;
