function err = RBF_euler(f, u_ex, a, b, u0, N, Meth)
%
% Usage:    err = RBF_euler(f, u_ex, a, b, u0, N)
% Purpose:  it solves the differential equation u'(t) = f(u(t),t) with
%           an initial condition u(a) = u0 and computes the global error 
%           using a RBF Euler method of the second order
% Input:    f = given function f(t,u) of the problem
%        u_ex = exact solution
%           a = initial point
%           b = end point 
%          u0 = initial value
%           N = controls the step size h
%        Meth = RBF Euler method: '1' MQ-RBF, '2' MA-RBF
% Output: err = global error at the final point t = b
%
h = (b-a)/N;
t = linspace(a,b,N+1)';
u = zeros(N+1, 1); 
u(1) = u0;
u(2) = u_ex(t(2));
for n = 2:N
    fn = f(t(n),u(n)); 
    fn_1 = f(t(n-1),u(n-1));
    if Meth == 1
        eps2 = (fn - fn_1) / (h * u(n));
        u(n+1) = (1 + (eps2*h^2)/2) * (u(n) + h*fn);
    elseif Meth == 2
        eps2 = -(3*fn - 3*fn_1) / (h * u(n));
        u(n+1) = (1 - (eps2*h^2)/6) * u(n) + (h + (eps2*h^3)/6) * fn;
    else
        warning('Method not available')
    end
end
err = abs(u(end) - u_ex(b));
