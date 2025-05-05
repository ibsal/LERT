function val = ceaInterp(cea,y,x,xq)
    if ~isfield(cea,y)||~isfield(cea,x), error('ceaInterp:BadField','Field not present'); end
    [xs,idx]=sort(cea.(x)); ys=cea.(y)(idx);
    val=interp1(xs,ys,xq,'linear','extrap'); if isscalar(xq), val=val(1); end
end
