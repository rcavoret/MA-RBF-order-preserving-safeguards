function check_invariants()
%CHECK_INVARIANTS Deterministic algebra, indexing, root and regression unit tests.
% No random inputs, optional toolboxes, or stored generated data are needed.
for h = [0.001,0.03]
    for delta = [0.1,0.2]
        for B = [-0.7,-0.1,0,0.06,0.2,0.8]
            A = 0.7*B+0.13; C = -0.3*B+0.04; E = 0.15; x = 1.2;
            D = A-2*B+C; R = B/max(delta^2,B^2); s = R*D; chi = 1-B*R;
            Porig = x+h*(3*A-B)/2+h*s*(41*A-31*B)/24;
            Pop = x+h*(23*A-16*B+5*C)/12+41*h*s*(A-B)/24;
            Qorig = x+h*(1/2-s/24)*(E+A);
            Qop = x+h*(6*E+5*A+2*B-C)/12-h*s*(E+A-2*B)/24;
            assert(abs(Pop-Porig-5*h*chi*D/12) < 2e-14);
            assert(abs(Qop-Qorig+h*chi*D/12) < 2e-14);
        end
    end
end
% Compare the preserved upstream MA routines with the factorized implementation.
p = manuscript_problem('example1');
for N = [16,32,64]
    methods = {'AB','PC'}; upstream_codes = [3,4];
    for k = 1:2
        old = RBF_adams(p.f,p.exact,p.a,p.b,p.u0,N,upstream_codes(k));
        [t,v,s] = MA_adams(p,N,methods{k},'none',0,'exact');
        m = error_measures(p,t,v,s);
        assert(abs(old-m.endpoint) < 5e-13);
    end
end
% Direct recurrence versus residual accumulation: catches the n -> n+2 offset.
for name = {'forcing','vector_log','vector_cancel'}
    p = manuscript_problem(name{1}); N = 32;
    for method = {'AB','AM','PC'}
        [t,v,s] = MA_adams(p,N,method{1},'order-preserving',1/N,'exact');
        residual = forcing_residuals(N,name{1},method{1},'order-preserving',1/N);
        direct = (p.exact(t)-v)*N^3;
        gap = max(abs(direct(:)-residual.scaled_errors(:)));
        assert(gap < 2e-9,'Direct/residual accumulation index or sign mismatch.');
        assert(s.completed);
    end
end
% Verify the linearization of the ACTUAL safeguarded maps, not formal targets.
h = 0.03; lambda = -1; z = h*lambda; eps0 = 1e-6; history = [0.4,-0.2,0.3];
p.a = 0; p.b = 3*h; p.u0 = eps0*history(1);
p.exact = @(t) eps0*(history(1)+(history(2)-history(1))*(t/h) ...
    +(history(3)-2*history(2)+history(1))*(t/h).*(t/h-1)/2);
p.f = @(t,x) lambda*x; p.J = lambda; p.affine_term = @(t) 0;
p.solve_kind = 'affine';
for placement = {'original','order-preserving'}
    op = strcmp(placement{1},'order-preserving');
    if op
        predicted = history(3)+z*(23*history(3)-16*history(2)+5*history(1))/12;
        base_am = (history(3)+z*(5*history(3)+2*history(2)-history(1))/12)/(1-z/2);
        base_pc = history(3)+z*(6*predicted+5*history(3)+2*history(2)-history(1))/12;
    else
        predicted = history(3)+z*(3*history(3)-history(2))/2;
        base_am = history(3)*(1+z/2)/(1-z/2);
        base_pc = history(3)+z*(predicted+history(3))/2;
    end
    expected = [predicted,base_am,base_pc]; methods = {'AB','AM','PC'};
    for k = 1:3
        [~,v,~] = MA_adams(p,3,methods{k},placement{1},0.2,'exact');
        assert(abs(v(end)/eps0-expected(k)) < 1e-9);
    end
end
% Root formulas used in Figure 1 and analytical negative-real checks.
for z = [-2,-0.5,0.2+0.7i,0.5-0.8i,2+1i]
    b = -(1+3*z/2+41*z^3/24); c = z/2+31*z^3/24;
    fast = max(abs([(-b+sqrt(b^2-4*c))/2,(-b-sqrt(b^2-4*c))/2]));
    assert(abs(fast-max(abs(roots([1,b,c])))) < 2e-12);
    b = -(1+3*z/2-7*z^3/24); c = z/2-17*z^3/24;
    fast = max(abs([(-b+sqrt(b^2-4*c))/2,(-b-sqrt(b^2-4*c))/2]));
    assert(abs(fast-max(abs(roots([1,b,c])))) < 2e-12);
end
ab3_poly = @(z) [1,-(1+23*z/12),4*z/3,-5*z/12];
assert(max(abs(roots(ab3_poly(-6/11+1e-4)))) < 1);
assert(max(abs(roots(ab3_poly(-6/11-1e-4)))) > 1);
euler_poly = @(z) [1,-(1+3*z/2-z^2/2),z/2-z^2/2];
assert(max(abs(roots(euler_poly(1-sqrt(3)+1e-4)))) < 1);
assert(max(abs(roots(euler_poly(1-sqrt(3)-1e-4)))) > 1);
fprintf('Algebra, upstream equivalence, residual indexing and stability checks passed.\n');
end
