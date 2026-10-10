function S = compensated_cumsum(X)
%COMPENSATED_CUMSUM Kahan cumulative summation, independently for each row.
% Used on scaled residuals, not on nearly equal O(1) numerical states.
S = zeros(size(X)); total = zeros(size(X,1),1); carry = total;
for k = 1:size(X,2)
    y = X(:,k)-carry;
    next = total+y;
    carry = (next-total)-y;
    total = next;
    S(:,k) = total;
end
end
