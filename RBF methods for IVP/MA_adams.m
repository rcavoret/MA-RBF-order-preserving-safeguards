function [t, v, stats] = MA_adams(p, N, method, placement, delta, startup)
%MA_ADAMS Computed MA-RBF AB, AM and single-correction PC recurrences.
% [t,v,stats] = MA_adams(p,N,method,placement,delta,startup)
%   method:    'AB', 'AM', or 'PC'.
%   placement: 'none', 'original', or 'order-preserving'.
%   delta:     fixed numeric threshold for this run (ignored for 'none').
%   startup:   'exact' or two classical 'rk4' steps.
% p comes from manuscript_problem. Columns of v correspond to entries of t.
% Original maps: manuscript labels MA-ab2-safe / MA-am1-safe / ABM_MA_scheme.
% Order-preserving maps: OPAB / OPAM / OPPC; vector extension: vectorAB/AM.
% n=1 is the first correction: k=3 below stores A=f_2,B=f_1,C=f_0.
if nargin < 6, startup = 'exact'; end
validateattributes(N, {'numeric'}, {'scalar','integer','>=',3});
if ~any(strcmp(method,{'AB','AM','PC'}))
    error('MARBF:Method','method must be AB, AM, or PC.');
end
if ~any(strcmp(placement,{'none','original','order-preserving'}))
    error('MARBF:Safeguard','Unknown safeguard placement.');
end
if ~any(strcmp(startup,{'exact','rk4'}))
    error('MARBF:Startup','startup must be exact or rk4.');
end
if ~strcmp(placement,'none') && ~logical(delta > 0)
    error('MARBF:Threshold','A safeguarded run needs delta > 0.');
end
d = numel(p.u0);
if d > 1 && ~strcmp(placement,'order-preserving')
    error('MARBF:VectorMethod','Systems use the order-preserving radial extension.');
end
h = (p.b-p.a)/N;
t = p.a+(0:N)*h;
v = repmat(0*p.u0(:),1,N+1);
v(:,1) = p.u0(:);
for j = 1:2
    if strcmp(startup,'exact')
        v(:,j+1) = p.exact(t(j+1));
    else
        k1 = p.f(t(j),v(:,j));
        k2 = p.f(t(j)+h/2,v(:,j)+h*k1/2);
        k3 = p.f(t(j)+h/2,v(:,j)+h*k2/2);
        k4 = p.f(t(j)+h,v(:,j)+h*k3);
        v(:,j+1) = v(:,j)+h*(k1+2*k2+2*k3+k4)/6;
    end
end
C = p.f(t(1),v(:,1)); B = p.f(t(2),v(:,2));
A = p.f(t(3),v(:,3));
stats.completed = true;
stats.stop_n = NaN;                 % Manuscript's ZERO-based denominator index.
stats.completed_steps = N;
stats.activations = 0;
stats.max_newton_iterations = 0;
stats.max_equation_residual = 0;
stats.reason = '';
for k = 3:N
    D = A-2*B+C;
    if strcmp(placement,'none')
        if logical(B == 0)
            stats.completed = false;
            stats.stop_n = k-2;
            stats.reason = 'undefined reciprocal';
            break
        end
        s = D/B;
    else
        normB2 = sum(B.*B);          % Euclidean radial reciprocal, not componentwise.
        active = logical(normB2 < delta^2);
        stats.activations = stats.activations+double(active);
        if active, denominator = delta^2; else, denominator = normB2; end
        s = sum(B.*D)/denominator;
    end
    x = v(:,k);
    op = strcmp(placement,'order-preserving');
    if ~strcmp(method,'AM')
        if op
            pred = x+h*(23*A-16*B+5*C)/12+41*h*s*(A-B)/24;
        else
            pred = x+h*(3*A-B)/2+h*s*(41*A-31*B)/24;
        end
    end
    if strcmp(method,'AB')
        z = pred;
    else
        a = h*(1/2-s/24);
        if op
            q = x+h*(5*A+2*B-C)/12-h*s*(A-2*B)/24;
        else
            q = x+a*A;
        end
        if strcmp(method,'PC')
            z = q+a*p.f(t(k+1),pred); % Exactly ONE corrector evaluation.
        else
            guess = x+h*A;          % Forward Euler guess, as in Table 1(b).
            [z, nit, res] = solve_implicit_step(p,t(k+1),q,a,guess);
            stats.max_newton_iterations = max(stats.max_newton_iterations,nit);
            stats.max_equation_residual = max(stats.max_equation_residual,res);
        end
    end
    if any(~isfinite(double(z)))
        stats.completed = false;
        stats.stop_n = k-2;
        stats.reason = 'nonfinite update';
        break
    end
    v(:,k+1) = z;
    C = B; B = A;
    A = p.f(t(k+1),z);              % Store the CORRECTED PC value, never the predictor.
end
if ~stats.completed
    stats.completed_steps = k-1;
    t = t(1:k);
    v = v(:,1:k);
end
end
