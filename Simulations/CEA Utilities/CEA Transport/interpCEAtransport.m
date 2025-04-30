function val = interpCEAtransport(tbl, Pbar, TK, OF, prop , mode , method)
%INTERPCEATRANSPORT  Interpolate any CEA transport/thermo property.
%
%   v = interpCEAtransport(tbl, Pbar, TK, OF, prop, mode, method)
%
% Required inputs
%   tbl    – table returned by parseCEAtransport
%   Pbar   – pressure  (scalar or vector)   [bar]
%   TK     – temperature                    [K]
%   OF     – mixture ratio                  [-]
%   prop   – string:  'rho','gamma','a','mu','k','Cp','Pr'
%
% Optional
%   mode   – 'eq' (default) or 'fr'.  Ignored for rho, gamma, a, mu.
%   method – 'linear' | 'nearest' | 'natural'  (default 'linear')
%
% Example
%   mu40  = interpCEAtransport(tbl,40,2600,5.5,'mu');
%   CpFr  = interpCEAtransport(tbl,50,3000,6,'Cp','fr');

arguments
    tbl       table
    Pbar      double
    TK        double
    OF        double
    prop      char   {mustBeMember(prop,{'rho','gamma','a','mu','k','Cp','Pr'})}
    mode      char   {mustBeMember(mode,{'eq','fr'})} = 'eq'
    method    char   {mustBeMember(method,{'linear','nearest','natural'})} = 'linear'
end

if ismember(prop,{'k','Cp','Pr'})
    sub = tbl(tbl.Mode==mode,:);
else
    sub = tbl;      % property does not depend on reaction model
end

F = scatteredInterpolant(sub.P_bar, sub.T_K, sub.OF, sub.(prop) , ...
                         method, method);
val = F(Pbar, TK, OF);
end
