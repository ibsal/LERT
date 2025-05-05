function tbl = parseCEAtransport(txtFile, toSI)
%PARSECEATRANSPORT  Convert a NASA-CEA TP output file to a tidy table.
%
%   tbl = parseCEAtransport('ceatransport.txt')            % keep CEA units
%   tbl = parseCEAtransport('ceatransport.txt', true)      % convert to SI
%
% Output columns            CEA-native units      →   SI when toSI == true
% ------------------------------------------------------------------------
%   OF        – mass ratio             [-]        unchanged
%   P_bar     – pressure               bar        Pa      (×1.0e5)
%   T_K       – temperature            K          K       (unchanged)
%   Mode      – 'eq' | 'fr'            —          —
%   rho       – density                kg/m³      kg/m³   (unchanged)
%   gamma     – specific-heat ratio    [-]        unchanged
%   a         – sonic velocity         m/s        m/s     (unchanged)
%   mu        – viscosity              milli-poise  → Pa·s (×1.0e-4)
%   k         – conductivity           mW/cm·K     → W/m·K (×0.1)
%   Cp        – specific heat          kJ/kg·K     → J/kg·K (×1.0e3)
%   Pr        – Prandtl number         [-]        unchanged
%
% The column names stay the same; only the numeric values change when
% toSI==true so downstream code does not break.

if nargin < 2, toSI = false; end

% ---------- read & parse (same logic as before) --------------------------
lines = strtrim(string(splitlines(fileread(txtFile))));
nL    = numel(lines);

rows  = {};
i     = 1;

isOF  = @(s) startsWith(s,'O/F=','IgnoreCase',true);
isP   = @(s) startsWith(s,'P,','IgnoreCase',true);
numV  = @(str) parseNums(str);
pad   = @(v,n) [v repmat(v(end),1,max(0,n-numel(v)))];
rx    = @(txt,label) regexp(txt,[label,'.*?([-+]?\d[\d.\-+Ee ]*)'], ...
                            'tokens','once','ignorecase');

while i<=nL
    if ~isOF(lines(i)),  i=i+1;  continue, end
    OF = str2double(regexp(lines(i),'O/F\s*=\s*([0-9.]+)','tokens','once'));
    i  = i+1;

    while i<=nL && ~isOF(lines(i))
        if ~isP(lines(i)), i=i+1; continue, end

        Pbar = numV(erase(lines(i),'P, BAR'));  Pbar = Pbar(1);
        i    = i+1;

        if i>nL || ~startsWith(lines(i),'T','IgnoreCase',true), break, end
        Tvec = numV(erase(lines(i),'T, K'));  nT = numel(Tvec);
        i    = i+1;

        chunk = "";
        j=i;
        while j<=nL && ~isP(lines(j)) && ~isOF(lines(j))
            chunk = chunk + newline + lines(j);
            j=j+1;
        end
        i = j;

        rho   = pad(numV(rx(chunk,'RHO')),     nT);
        gamma = pad(numV(rx(chunk,'GAMMAs')),  nT);
        avel  = pad(numV(rx(chunk,'SON VEL')), nT);
        mu    = pad(numV(rx(chunk,'VISC')),    nT);

        eqTok = regexp(chunk,'EQUILIBRIUM\s+REACTIONS(.*?)WITH\s+FROZEN', ...
                       'tokens','once','ignorecase');
        frTok = regexp(chunk,'FROZEN\s+REACTIONS(.*)$', ...
                       'tokens','once','ignorecase');
        if isempty(eqTok)||isempty(frTok), warning('Missing EQ/FR'); continue, end
        eqSec = eqTok{1};  frSec = frTok{1};

        Cp_eq = pad(numV(rx(eqSec,'Cp')),          nT);
        k_eq  = pad(numV(rx(eqSec,'CONDUCTIVITY')),nT);
        Pr_eq = pad(numV(rx(eqSec,'PRANDTL')),     nT);

        Cp_fr = pad(numV(rx(frSec,'Cp')),          nT);
        k_fr  = pad(numV(rx(frSec,'CONDUCTIVITY')),nT);
        Pr_fr = pad(numV(rx(frSec,'PRANDTL')),     nT);

        for t = 1:nT
            rows(end+1,:) = {OF,Pbar,Tvec(t),'eq', ...
                             rho(t),gamma(t),avel(t),mu(t), ...
                             k_eq(t),Cp_eq(t),Pr_eq(t)}; %#ok<AGROW>
            rows(end+1,:) = {OF,Pbar,Tvec(t),'fr', ...
                             rho(t),gamma(t),avel(t),mu(t), ...
                             k_fr(t),Cp_fr(t),Pr_fr(t)}; %#ok<AGROW>
        end
    end
end

tbl = cell2table(rows, 'VariableNames', ...
      {'OF','P_bar','T_K','Mode','rho','gamma','a','mu','k','Cp','Pr'});
tbl.Mode = categorical(tbl.Mode);

% ---------- optional unit conversion to SI ------------------------------
if toSI
    tbl.P_bar = tbl.P_bar * 1e5;     % bar  →  Pa
    tbl.mu    = tbl.mu    * 1e-4;    % mPoise → Pa·s
    tbl.k     = tbl.k     * 0.1;     % mW/cm·K → W/m·K
    tbl.Cp    = tbl.Cp    * 1e3;     % kJ/kg·K → J/kg·K
end
end
% ========================================================================
function v = parseNums(str)
    if isempty(str) || all(str==""), v=[]; return, end
    if iscell(str), str=str{1}; end
    str = regexprep(str,'(?<=\d)([+-])(?=\d)','E$1');
    v   = sscanf(str,'%f').';
end
