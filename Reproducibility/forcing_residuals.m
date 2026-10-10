function out = forcing_residuals(N, problem, method, placement, delta)
%FORCING_RESIDUALS Exact quadratic-forcing residuals, scaled by h^3.
% This is an algebraically equivalent error computation, NOT a fitted model.
% Quadratic forcing is integrated exactly by each order-three polynomial base.
% For state-independent f, errors equal cumulative local residuals and AM=PC.
% The direct 80-digit recurrence identity is checked by script_verify_80digits.
% problem: 'forcing', 'vector_log', 'vector_cancel'.
% ONLY scalar forcing admits the manuscript's original safeguard extension.
validateattributes(N,{'numeric'},{'scalar','integer','>=',4});
if mod(N,2) ~= 0, error('MARBF:EvenMesh','This diagnostic needs an even N.'); end
if ~any(strcmp(method,{'AB','AM','PC'})), error('MARBF:Method','Unknown method.'); end
if ~any(strcmp(placement,{'original','order-preserving'}))
    error('MARBF:Safeguard','Residual diagnostic requires a safeguard.');
end
switch problem
    case 'forcing',       g = 1; q = 1;
    case 'vector_log',    g = [1;0]; q = [1;1];
    case 'vector_cancel', g = [1;0]; q = [0;1];
    otherwise, error('MARBF:Forcing','Unknown quadratic forcing.');
end
if numel(g) > 1 && strcmp(placement,'original')
    error('MARBF:VectorMethod','Original safeguards are scalar in the manuscript.');
end
h = 1/N;
t = -1/2+(1:N-2)/N;             % Denominators n=1,...,N-2, NOT update endpoints.
F = g*t+q*(t.^2);
norm2 = sum(F.^2,1);
den = max(delta^2,norm2);
Rq = (q.'*F)./den;              % <R_delta(F),q>; one radial denominator.
if strcmp(method,'AB')
    residual = -(41/12)*h*bsxfun(@times,Rq,bsxfun(@plus,g,q*(2*t+h)));
else
    residual = (1/4)*h*bsxfun(@times,Rq,bsxfun(@plus,g,q*(2*t+5*h/3)));
end
if strcmp(placement,'original')
    chi = 1-norm2./den;
    if strcmp(method,'AB'), c = 5/6; else, c = -1/6; end
    residual = residual+c*q*chi;
end
scaled = [zeros(numel(g),3),compensated_cumsum(residual)];
out.zero = scaled(:,N/2+1);     % Three zero starting errors, at j=0,1,2.
out.endpoint = scaled(:,end);
out.maximum = max(sqrt(sum(scaled.^2,1)));
out.scaled_errors = scaled;
out.h = h;
end
