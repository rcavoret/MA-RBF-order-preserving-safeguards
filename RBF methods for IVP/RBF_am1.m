function [err, t, v, stats] = RBF_am1(f, df, exact, a, b, u0, N, kernel)
%RBF_AM1 Standalone AM1 comparison in Table 1(b), not a predictor-corrector.
% kernel is 'classical', 'MQ', or 'MA'. All three use exact v_1 and v_2.
% MQ and MA use the SAME Taylor-expanded weight (AM1_common_weight),
% so identical histories and the same nonlinear root give identical values.
if ~any(strcmp(kernel,{'classical','MQ','MA'}))
    error('MARBF:Kernel','kernel must be classical, MQ, or MA.');
end
h = (b-a)/N;
t = a+(0:N)'*h;
v = zeros(N+1,1);
v(1) = u0; v(2) = exact(t(2)); v(3) = exact(t(3));
p.f = f; p.jacobian = df; p.solve_kind = 'newton';
p.newton_tolerance = 1e-14; p.max_newton_iterations = 30;
C = f(t(1),v(1)); B = f(t(2),v(2)); A = f(t(3),v(3));
stats.max_newton_iterations = 0; stats.max_equation_residual = 0;
for k = 3:N
    gamma = h/2;
    if ~strcmp(kernel,'classical')
        if B == 0, error('MARBF:UndefinedReciprocal','Zero f_n in AM1 comparison.'); end
        D = A-2*B+C;
        gamma = h*(1/2-D/(24*B));
    end
    q = v(k)+gamma*A;
    [v(k+1),nit,res] = solve_implicit_step(p,t(k+1),q,gamma,v(k)+h*A);
    stats.max_newton_iterations = max(stats.max_newton_iterations,nit);
    stats.max_equation_residual = max(stats.max_equation_residual,res);
    C = B; B = A; A = f(t(k+1),v(k+1));
end
err = abs(v(end)-exact(b));
end
