function [z, iterations, residual] = solve_implicit_step(p, t, q, a, guess)
%SOLVE_IMPLICIT_STEP Solve z=q+a*f(t,z) using the manuscript's local root.
% Quadratic and affine cases are solved directly; Example 1 uses Newton.
% No clipping, restart, denominator floor, or step rejection is used.
iterations = 0;
switch p.solve_kind
    case 'quadratic'
        c = p.quadratic_coefficient(t);
        disc = 1-4*a*c*q;
        if ~logical(disc >= 0)
            error('MARBF:QuadraticRoot','The local quadratic root is not real.');
        end
        z = 2*q/(1+sqrt(disc));
    case 'affine'
        d = numel(q);
        z = (eye(d)-a*p.J) \ (q+a*p.affine_term(t));
    case 'forcing'
        z = q+a*p.f(t,guess);
    case 'newton'
        z = guess;
        converged = false;
        for iterations = 1:p.max_newton_iterations
            G = z-q-a*p.f(t,z);
            M = eye(numel(z))-a*p.jacobian(t,z);
            dz = -(M\G);
            z = z+dz;
            if max(abs(double(dz))) <= double(p.newton_tolerance)
                converged = true;
                break
            end
        end
        if ~converged
            error('MARBF:NewtonFailure', ...
                'Newton did not meet the increment tolerance at t=%.17g.',double(t));
        end
    otherwise
        error('MARBF:ImplicitSolver','Unknown implicit solver: %s',p.solve_kind);
end
if any(~isfinite(double(z)))
    error('MARBF:NonfiniteRoot','The implicit step produced a nonfinite value.');
end
residual = max(abs(double(z-q-a*p.f(t,z))));
end
