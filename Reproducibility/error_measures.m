function m = error_measures(p, t, v, stats)
%ERROR_MEASURES Maximum nodal Euclidean error and endpoint error.
% Subtract in the working arithmetic BEFORE converting high-precision values.
% No full-interval error is assigned to an incomplete run.
if nargin == 4 && ~stats.completed
    m.endpoint = NaN; m.maximum = NaN; m.signed_endpoint = NaN;
    return
end
e = double(p.exact(t)-v);
nodal = sqrt(sum(e.^2,1));
m.endpoint = nodal(end);
m.maximum = max(nodal);
m.signed_endpoint = e(:,end);
end
