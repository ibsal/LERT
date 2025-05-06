function M = machFromArea(Rstation, gamma, Rthroat)
%MACHFROMAREA  Convert radius profile to Mach number profile (1‑D nozzle).
%
%   M = machFromArea(Rstation, gamma)         % throat inferred (min R)
%   M = machFromArea(Rstation, gamma, Rt)     % throat radius supplied
%
%   Rstation : vector of radii along axis, from chamber to exit (m)
%   gamma    : specific‑heat ratio (scalar)
%   Rthroat  : optional throat radius. If not given, min(Rstation).
%
%   The routine finds the throat index (first occurrence of R = Rthroat).
%   For all points *upstream* of that index it solves the **subsonic**
%   branch of the isentropic A/A* relation; for the throat Mach = 1; and
%   for points *downstream* it solves the **supersonic** branch.
%
%   Uses fzero with robust bracketing.
% -----------------------------------------------------------------------

    narginchk(2,3);
    Rstation = Rstation(:);                % column vector
    if nargin < 3 || isempty(Rthroat)
        Rthroat = min(Rstation);           % infer throat
    end

    % first index where radius equals the throat (within tol)
    tol      = 1e-9;
    throatIdx = find(abs(Rstation - Rthroat) < tol, 1, 'first');
    if isempty(throatIdx)
        error('machFromArea:NoThroat','Rthroat not found in Rstation vector');
    end

    % preallocate
    M = zeros(size(Rstation));

    % anonymous for A/A* in terms of Mach
    Afun = @(M) (1./M).*((2/(gamma+1)).*(1 + (gamma-1)/2.*M.^2)).^((gamma+1)/(2*(gamma-1)));

    % iterate through stations
    for i = 1:numel(Rstation)
        Ar = (Rstation(i)/Rthroat)^2;      % A/A* via radius^2
        if i < throatIdx                   % upstream -> subsonic
            if Ar == 1                     % rare flat plateau before min
                M(i) = 1;
            else
                fun = @(M) Afun(M) - Ar;
                M(i) = fzero(fun, [1e-3 0.99999999999999]);
            end
        elseif i == throatIdx              % throat
            M(i) = 1;
        else                               % downstream -> supersonic
            fun = @(M) Afun(M) - Ar;
            M(i) = fzero(fun, [1.000000000001 10]);
        end
    end
end