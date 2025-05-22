function val = ceaGridInterp(csvFile, field, Pc_bar, OF)
% Fast scatteredInterpolant cache for any property in the CSV grid
persistent Fns CacheFile
if isempty(Fns) || ~strcmp(CacheFile,csvFile)
    T = readtable(csvFile);
    CacheFile = csvFile;
    props = T.Properties.VariableNames;
    for p = props
        if ismember(p{1},{'Pc','O_F'}), continue, end
        Fns.(p{1}) = scatteredInterpolant(T.Pc, T.O_F, T.(p{1}), 'natural');
    end
end
if ~isfield(Fns,field)
    error('ceaGridInterp:badfield','Field "%s" not found in %s',field,csvFile);
end
val = Fns.(field)(Pc_bar, OF);
end
