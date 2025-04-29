% ---------------- interpolation utility ----------------------------------
function val = ceaInterp(cea,yField,xField,xq)
    arguments, cea struct, yField (1,:) char, xField (1,:) char, xq double, end
    if ~isfield(cea,yField)||~isfield(cea,xField)
        error("ceaInterp:BadField","Field not present in CEA struct.");
    end
    [xs,idx] = sort(cea.(xField));
    ys       = cea.(yField)(idx);
    val      = interp1(xs,ys,xq,'linear','extrap');
    if isscalar(xq), val = val(1); end
end
