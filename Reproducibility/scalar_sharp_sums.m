function out = scalar_sharp_sums(N, one)
%SCALAR_SHARP_SUMS Exact residual sums for F=t+t^2 and delta=h, in VPA.
% Returns N^3 times signed errors at t=0 and t=1/2 (AB, AM/PC).
% Partial fractions of sharpABres/sharpAMres reduce the sums to harmonics.
% H_m=psi(m+1)-psi(1) evaluates the SAME finite sums efficiently at 80 digits.
% This avoids both long symbolic loops and subtraction of close O(1) states.
if mod(N,2) ~= 0 || N < 8, error('MARBF:Mesh','Use even N >= 8.'); end
M = N/2; h = one/N;
H = @(m) psi((m+1)*one)-psi(one);
c = [one,5*one/3]; factors = [-41*one/12,one/4];
out.zero = repmat(0*one,1,2); out.endpoint = out.zero;
for k = 1:2
    % At the node 0 the last correction center is t=-2h. No active term
    % contributes yet when delta=h; both placements give the same zero error.
    at_zero = (1+c(k)*h)*(1-H(M-1))+(1-c(k)*h)*(H(N-2)-H(M));
    % All centers -1/2+h,...,1/2-2h are summed. Only t=-h and t=0
    % are active: at -h, R_delta(F)=1-N; the t=0 residual is zero.
    at_end = (1+c(k)*h)*(1-one/(M-1)) ...
        +(1-c(k)*h)*(H(N-2)-H(M)+H(3*M-2)-H(N)) ...
        +(h-1)*(1+(c(k)-2)*h);
    out.zero(k) = factors(k)*at_zero;
    out.endpoint(k) = factors(k)*at_end;
end
out.original_zero = out.zero;
chi_sum = 1+2*h-h^2;
out.original_endpoint = out.endpoint+chi_sum*[5*one/6,-one/6];
end
