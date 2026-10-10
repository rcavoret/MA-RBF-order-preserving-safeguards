function [err, t, u] = ab2(f, u_ex, a, b, u0, N)
%
% Usage:    err = ab2(f, u_ex, a, b, u0, N)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t)
%           with an initial condition u(a) = u0 and computes the global 
%           error using the two-step Adams-Bashforth method
% Input:    f = given function f(t,u) of the problem
%        u_ex = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
% Output: err = endpoint error; optional t and u contain the full trajectory
%
h = (b - a) / N;
t = a + (0:N)' * h;
u = zeros(N+1, 1);
u(1) = u0;
u(2) = u_ex(t(2));
t2 = f(t(1), u(1));
for n = 2:N
    t1 = f(t(n), u(n));
    u(n+1) = u(n) + (h/2) * (3*t1 - t2);
    if n < N
       t2 = t1; 
    else
       break    
    end
end
err = abs(u(end) - u_ex(b));
