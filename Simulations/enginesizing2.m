% ---------- user inputs --------------------------------------------------
CombustionPressure = 300 * 6894.76;     % Pa
AmbientPressure    = 14.6959 * 6894.76; % Pa
CombustionGamma    = 1.2;

Thrust   = 250 * 4.44822;     % N
OFRatio  = [5];           % *** vector of mixture ratios ***

% ---------- look-up tables ----------------------------------------------
cea = readCEA('ceaoutput.txt');

% ---------- pre-plot setup ----------------------------------------------
figure
subplot(2,1,1); hold on; axis equal; grid on
xlabel('Axial position  x'), ylabel('Radius  r')
title('Nozzle / chamber profiles')

subplot(2,1,2); hold on; grid on
yyaxis left
ylabel('T  [K]')
xlabel('Axial position  x')
title('Static temperature along centre-line')

colors = lines(numel(OFRatio));      % colour map for each OF

% ---------- loop over each mixture ratio --------------------------------
for k = 1:numel(OFRatio)
    OF = OFRatio(k);

    % --- pull CEA data -------------
    Mexit   = ceaInterp(cea,'Mach_e','O_F',OF);
    aexit   = ceaInterp(cea,'a_e',  'O_F',OF);
    Cstar   = ceaInterp(cea,'Cstar','O_F',OF);
    Ae_At   = ceaInterp(cea,'Ae_At','O_F',OF);
    Tc      = ceaInterp(cea,'T_c',  'O_F',OF);
    gamma_e = ceaInterp(cea,'gamma_e','O_F',OF);

    % --- mass-flow, throat & exit ---
    Vexit   = Mexit * aexit;
    mdot    = Thrust / Vexit;
    At      = mdot * Cstar / CombustionPressure;        % m²
    rt      = sqrt(At/pi);                              % m
    re      = sqrt(Ae_At*At/pi);                        % m

    % --- rough 1.8 Dt nozzle length (in) --------------
    rt_in   = rt * 39.3701;  re_in = re * 39.3701;
    Lnoz_in = re_in * 1.8 * 2;

    % --- combustion-chamber profile (in) --------------
    CL_in = 3; CR_in = .65;                            % keeping your values
    x = 0:0.01:(CL_in + Lnoz_in);
    [r_in, ~] = combustionChamberProfile( ...
               x, CR_in, 0.25, 45, 0.25, rt_in, CL_in, ...
               15, re_in, 2, Lnoz_in);

    % --- Mach & temperature curve ---
    Mcurve = machFromArea(r_in, gamma_e, rt_in);
    Tcurve = Tc ./ (1 + (gamma_e-1)/2 .* Mcurve.^2);

    % ---- plot nozzle outline -------
    subplot(2,1,1)
    plot(x,  r_in,  'Color',colors(k,:),'LineWidth',1.2)
    plot(x, -r_in,  'Color',colors(k,:),'LineWidth',1.2)
    axis equal

    % ---- plot temperature ----------
    subplot(2,1,2)
    yyaxis left
    plot(x, Tcurve, 'Color',colors(k,:), 'LineWidth',1.2, ...
         'DisplayName',sprintf('O/F = %.1f',OF))
end

% ---------- finish twin-axis °F scale -----------------------------------
subplot(2,1,2)
yyaxis left
ylK = ylim;                                       % limits covering all curves
yyaxis right
ylim( (ylK - 273.15)*9/5 + 32 )
ylabel('T  [°F]')

legend('Location','best')
