function [err, t, u] = euler(f, u_ex, a, b, u0, N)
%
% Usage:    err = euler(f, u_ex, a, b, u0, N)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t) with
%           an initial condition u(a) = u0 and computes the global error
%           using the Euler method
% Input:    f = given function f(t,u) of the problem
%        u_ex = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
% Output: err = endpoint error; optional t and u contain the full trajectory
%
h = (b-a) / N;
t = a + (0:N)' * h;
u = zeros(N+1, 1); 
u(1) = u0;
for n = 1:N
    u(n+1) = u(n) + h*f(t(n),u(n));
end
err = abs(u(end) - u_ex(b));


