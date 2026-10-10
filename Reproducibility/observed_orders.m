function rates = observed_orders(N, errors)
%OBSERVED_ORDERS Use unrounded errors; allow non-dyadic companion meshes.
% errors has one row per mesh and one column per method/threshold.
N = N(:);
if size(errors,1) ~= numel(N)
    error('MARBF:Dimensions','errors must have one row per mesh.');
end
rates = NaN(size(errors));
for j = 2:numel(N)
    rates(j,:) = log(errors(j-1,:)./errors(j,:))/log(N(j)/N(j-1));
end
end
