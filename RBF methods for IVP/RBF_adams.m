function err = RBF_adams(f, u_ex, a, b, u0, N, Meth)
%
% Usage:    err = MQ_ab2(f, u_ex, a, b, u0, N)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t) with
%           an initial condition u(a) = u0 and computes the global error
%           using a RBF Adams method of the third order
% Input:    f = given function f(t,u) of the problem
%        u_ex = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
%        Meth = RBF Adams method: '1' MQ-RBF two-step Adams Bashforth, 
%               '2'MQ-RBF two step Adams-Bashforth one-step Adams Moulton
%               predictor-corrector, '3' MA-RBF two-step Adams Bashforth, 
%               '4' MA-RBF two step Adams-Bashforth one-step Adams Moulton
%               predictor-corrector
% Output: err = global error at the final point t = b
%
h = (b-a)/N;
t = linspace(a, b, N+1)';
u = zeros(N+1, 1); 
u(1) = u0;
u(2) = u_ex(t(2));
u(3) = u_ex(t(3));
t3 = f(t(1), u(1)); 
t2 = f(t(2), u(2)); 
for n = 3:N
    t1 = f(t(n), u(n));
    if Meth == 1 % MQ ab2
        eps2 = (t1 - 2*t2 + t3) / (h^2 * t2);
        u(n+1) = u(n) + (3/2*h - 7/24*eps2*h^3)*t1 + (-1/2*h + 17/24*eps2*h^3)*t2;
    elseif Meth == 2 % MQ ab2am1
        eps2 = (t1 - 2*t2 + t3) / (h^2 * t2);
        u_p= u(n) + (3/2*h - 7/24*eps2*h^3)*t1 + (-1/2*h + 17/24*eps2*h^3)*t2;
        t0 = f(t(n+1),u_p);
        u(n+1) = u(n) + (h/2 - eps2/24*h^3) * (t0 + t1); 
    elseif Meth == 3 % MA ab2
        eps2 = -(3*t1 - 6*t2 + 3*t3) / (h^2 * t2);
        u(n+1) = u(n) + (3/2*h - 41/72*eps2*h^3)*t1 + (-1/2*h + 31/72*eps2*h^3)*t2;
    elseif Meth == 4 % MA ab2am1
        eps2 = -(3*t1 - 6*t2 + 3*t3) / (h^2 * t2);
        u_p= u(n) + (3/2*h - 41/72*eps2*h^3)*t1 + (-1/2*h + 31/72*eps2*h^3)*t2;
        t0 = f(t(n+1),u_p);
        u(n+1) = u(n) + (h/2 + eps2/72*h^3) * (t0 + t1); 
    else 
        warning('Method not available')
    end
    if n < N
       t3 = t2; 
       t2 = t1; 
    else
       break
    end
end
err = abs(u(end) - u_ex(b));

