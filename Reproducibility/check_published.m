function check_published(name, errors, rates, counts, high_precision)
%CHECK_PUBLISHED Compare computed data with published rounding intervals.
% Double precision gets a small, explicit roundoff allowance.
if nargin < 4, counts = []; end
if nargin < 5, high_precision = false; end
r = published_values(name);
if ~isequal(size(errors),size(r.errors)) || ~isequal(size(rates),size(r.rates))
    error('MARBF:ReferenceDimensions','Unexpected array size for %s.',name);
end
half_unit = 0.500001*10.^(floor(log10(abs(r.errors)))-r.error_digits);
if high_precision, roundoff = 0; else, roundoff = 5e-14; end
if any(~isfinite(errors(:))) || any(abs(errors(:)-r.errors(:)) > half_unit(:)+roundoff)
    error('MARBF:ErrorMismatch','Errors differ from the manuscript for %s.',name);
end
valid = isfinite(r.rates);
% 0.5 unit is the rounding interval; doubles allow one final displayed unit.
if high_precision, tol = 5.0001e-5; else, tol = 1e-4; end
if any(~isfinite(rates(valid))) || any(abs(rates(valid)-r.rates(valid)) > tol)
    error('MARBF:RateMismatch','Rates differ from the manuscript for %s.',name);
end
if ~isempty(r.counts) && ~isequal(counts,r.counts)
    error('MARBF:CountMismatch','Activation counts differ for %s.',name);
end
fprintf('Manuscript regression check passed: %s\n',name);
end
