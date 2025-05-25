function [M, T, P] = nozzleMach(Pc, Pb, gamma, A, Tc)
%NOZZLEMACH  Compute 1‑D Mach‑number distribution along a nozzle.
%
%   M = nozzleMach(Pc, Pb, gamma, A)
%
%   Inputs
%   ───────────────────────────────────────────────
%   Pc     – Chamber / stagnation pressure  [Pa]
%   Pb     – Back (ambient) pressure        [Pa]
%   gamma  – Ratio of specific heats (Cp/Cv)
%   A      – Vector of cross‑sectional areas at each axial station [m²]
%            (include both converging and diverging sections)
%
%   Output
%   ───────────────────────────────────────────────
%   M      – Vector of Mach numbers at the same stations as A
%
%   The routine decides automatically whether the nozzle operates:
%     1) Entirely subsonic (unchoked)
%     2) Just choked (M=1 at throat, subsonic everywhere else)
%     3) Choked with isentropic supersonic flow to the exit
%     4) Choked with a normal shock in the divergent section
%
%   Assumptions
%   ───────────────────────────────────────────────
%   • 1‑D, steady, perfect‑gas, adiabatic, no friction.
%   • Normal shock, if present, obeys ideal shock relations.
%   • The throat is the minimum area in A.
%
%   Author: <your name>   Date: <date>

% ─── Pre‑processing ────────────────────────────────────────────────────────
A  = A(:).';                    % row vector
Nt = find(A == min(A), 1);      % throat index
At = A(Nt);
Ar = A ./ At;                    % area ratio A/A*

% helper ↦ invert area‑Mach
area2mach = @(AR, branch) fzero( ...
    @(M) ((1./M).* ((2/(gamma+1))*(1+(gamma-1)/2*M.^2)) ...
           .^((gamma+1)/(2*(gamma-1)))) - AR, ...
    branch);                    % branch=[1e-6, 0.9999999999] → sub, branch=3 → sup

% ─── Critical pressure ratio for choking ───────────────────────────────────
Pstar_by_P0 = (2/(gamma+1))^(gamma/(gamma-1));   % P* / P0
isChoked = (Pb/Pc) < Pstar_by_P0;

% ───  CASE 0: Un‑choked (fully subsonic)  ──────────────────────────────────
if ~isChoked
    M = arrayfun(@(ARi) area2mach(ARi, [1e-6, 0.9999999999]), Ar);  % subsonic root
    T          = Tc ./ (1 + (gamma-1)/2 .* M.^2);
    P          = Pc .* (1 + (gamma-1).*0.5 .* M.^2) ...
                                          .^(-gamma./(gamma-1));
    return
end

% ─── Downstream (exit) solutions for a choked nozzle ───────────────────────
Ae_At  = Ar(end);
Me_sub = area2mach(Ae_At, [1e-6, 0.9999999999]);     % sub‑branch
Me_sup = area2mach(Ae_At, [1.0000001, 20]);     % sup‑branch

Pe_sub = Pc * (1+(gamma-1)/2*Me_sub^2)^(-gamma/(gamma-1));
Pe_sup = Pc * (1+(gamma-1)/2*Me_sup^2)^(-gamma/(gamma-1));

% Decide among the three choked regimes
if Pb >= Pe_sub                       % (1) JUST CHOKED, fully subsonic
    M = arrayfun(@(ARi) area2mach(ARi, [1e-6, 0.9999999999]), Ar);
    M(Nt) = 1;
    disp("JUST CHOKED)")
    T          = Tc ./ (1 + (gamma-1)/2 .* M.^2);
    P          = Pc .* (1 + (gamma-1).*0.5 .* M.^2) ...
                                          .^(-gamma./(gamma-1));
    % throat exactly sonic
    return
elseif Pb <= Pe_sup                   % (2) ISENTROPIC SUPERSONIC
    M = zeros(size(A));
    for k = 1:numel(A)
        if k < Nt,  M(k) = area2mach(Ar(k), [1e-6, 0.9999999999]);      % converging (sub)
        elseif k==Nt, M(k)=1;
        else        , M(k) = area2mach(Ar(k), [1.0000001, 20]);    % diverging (sup)
        end
    end
    T          = Tc ./ (1 + (gamma-1)/2 .* M.^2);
    P          = Pc .* (1 + (gamma-1).*0.5 .* M.^2) ...
                                          .^(-gamma./(gamma-1));
    return
end

% ─── (3) NORMAL‑SHOCK case  ───────────────────────────────────────────────
% Find the station where a normal shock makes Pb match exit pressure.
% Sweep through divergent section until best match.
bestErr = inf;  idxShock = NaN;
disp("SHOCK")
Pstag2 = 0;
for j = Nt+1:numel(A)-1
    % upstream supersonic Mach at station j
    M1 = area2mach(Ar(j), [1.0000001, 20]);
    % Normal‑shock relations
    M2s = sqrt((1 + ((gamma-1)/2) * M1^2)/(gamma * M1^2 - (gamma -1)/2));
    P2_by_P1 = 1 + ((2*gamma)/(gamma+1)) * (M1^2 - 1);
    P1       = Pc * (1+(gamma-1)/2*M1^2)^(-gamma/(gamma-1));
    P2       = P1 * P2_by_P1;           % static P immediately after shock
    Pstarg2 = P2/((1 + 0.5*(gamma-1) * M2s^2)^((-1*gamma)/(gamma-1)));
    % --- NEW critical area after the shock -----------------------------
    F_sub  = @(M) (1./M) .* ((2/(gamma+1))*(1+(gamma-1)/2*M.^2)) ...
                              .^((gamma+1)/(2*(gamma-1)));
    Astar2 = A(j) / F_sub(M2s);           % <-- this line replaces AR_down = ...
    % downstream area ratios referenced to the NEW A*
    AR_down  = A(j+1:end)/Astar2;
    Me_vec   = arrayfun(@(ARi) area2mach(ARi,[1e-6,0.9999999]), AR_down);
    Pe       = Pstarg2 * (1 + 0.5 * (gamma-1) * Me_vec(end)^2)^((-1*gamma)/(gamma-1));
    err = abs(Pe - Pb);
    if err < bestErr
        M2 = M2s;
        bestErr = err;
        idxShock = j;
        M_downstream = Me_vec;
        Pstag2 = Pstarg2;
    end
end

if isnan(idxShock)
    error('Could not locate a normal‑shock position with the given Pb.');
end

% Build final Mach vector
M = zeros(size(A));
for k = 1:idxShock-1                % upstream: sonic + supersonic
    if k < Nt,  M(k) = area2mach(Ar(k), [1e-6, 0.9999999999]);
    elseif k==Nt, M(k)=1;
    else         M(k) = area2mach(Ar(k), [1.0000001, 20]);
    end
end

M(idxShock) = area2mach(Ar(idxShock), [1.0000001, 20]);  % M1 just before shock
M(idxShock+1) = M2;                          % immediately after shock
M(idxShock+2:end) = M_downstream(2:end);     % further subsonic expansion
T = zeros(size(M));
P = zeros(size(M));
T          = Tc ./ (1 + (gamma-1)/2 .* M.^2);
P(1:idxShock)          = Pc .* (1 + (gamma-1).*0.5 .* M(1:idxShock).^2) ...
                                          .^(-gamma./(gamma-1));
P((idxShock+1):end) = Pstag2 .* (1 + 0.5.*(gamma-1) .* M((idxShock+1):end).^2).^((-1*gamma)/(gamma-1));

end
