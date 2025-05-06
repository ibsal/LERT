function val = interpCEAtransport(tbl, Pbar, TK, OF, prop, mode, method)
    arguments
      tbl       table
      Pbar      double
      TK        double
      OF        double
      prop      char   {mustBeMember(prop,{'rho','gamma','a','mu','k','Cp','Pr'})}
      mode      char   {mustBeMember(mode,{'eq','fr'})} = 'eq'
      method    char   {mustBeMember(method,{'linear','nearest','natural'})} = 'linear'
    end

    %----------------------------------------------------------------------%
    % Build a key based on prop, mode, method so we only do each once:
    persistent F_cache
    if isempty(F_cache)
      F_cache = containers.Map;
    end
    key = sprintf('%s_%s_%s', prop, mode, method);

    % If we haven’t built this one yet, do it now and stash it:
    if ~F_cache.isKey(key)
      if ismember(prop,{'k','Cp','Pr'})
        sub = tbl(tbl.Mode==mode,:);
      else
        sub = tbl;
      end

      F = scatteredInterpolant(...
            sub.P_bar, sub.T_K, sub.OF, sub.(prop), ...
            method, method);
      F_cache(key) = F;
    else
      F = F_cache(key);
    end
    %----------------------------------------------------------------------%

    % Now just call the (cached) interpolant:
    val = F(Pbar, TK, OF);
end
