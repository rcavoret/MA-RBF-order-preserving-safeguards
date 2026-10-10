function p = manuscript_problem(name, one)
%MANUSCRIPT_PROBLEM Exact IVPs in Section 5; ONE fixes the arithmetic type.
% Use one=1 for double; one=vpa('1') after digits(80) for 80-digit arithmetic.
% Rational data are formed AFTER selecting the type, never from rounded doubles.
if nargin < 2, one = 1; end
p.name = name;
p.newton_tolerance = one/10^14;
p.max_newton_iterations = 30;
switch name
    case 'example1'
        p.a = 0*one; p.b = one; p.u0 = one;
        p.f = @(t,x) -x.^2;
        p.exact = @(t) one./(one+t);
        p.jacobian = @(t,x) -2*x;
        p.solve_kind = 'newton';
    case 'example2'
        p.a = one; p.b = 3*one; p.u0 = 2*one;
        p.f = @(t,x) (2*t.^2-x)./(t.^2.*x-t);
        p.exact = @(t) one./t+sqrt(one./t.^2+4*t-4);
        p.jacobian = @(t,x) (t-2*t.^4)./(t.^2.*x-t).^2;
        p.solve_kind = 'newton';
    case {'example3','control'}
        p.a = -one/2; p.b = one/2;
        if strcmp(name,'control'), p.a = one/4; end
        p.f = @(t,x) t.*(one+t).*x.^2;
        p.exact = @(t) one./(one-t.^2/2-t.^3/3);
        p.u0 = p.exact(p.a);
        p.jacobian = @(t,x) 2*t.*(one+t).*x;
        p.quadratic_coefficient = @(t) t.*(one+t);
        p.solve_kind = 'quadratic';
    case 'forcing'
        p.a = -one/2; p.b = one/2;
        p.f = @(t,x) t+t.^2;
        p.exact = @(t) one+t.^2/2+t.^3/3;
        p.u0 = p.exact(p.a);
        p.solve_kind = 'forcing';
    case {'vector_log','vector_cancel','coupled'}
        p.a = -one/2; p.b = one/2;
        p.exact = @(t) [one+t.^2/2+t.^3/3; 2*one+t.^3/3];
        p.forcing = @(t) [t+t.^2; t.^2];
        if strcmp(name,'vector_cancel')
            p.exact = @(t) [one+t.^2/2; 2*one+t.^3/3];
            p.forcing = @(t) [t; t.^2];
        end
        p.f = @(t,x) p.forcing(t);
        p.solve_kind = 'forcing';
        if strcmp(name,'coupled')
            p.J = [-one, 3*one/10; -one/5, -7*one/10];
            exact = p.exact; F = p.forcing; J = p.J;
            p.affine_term = @(t) F(t)-J*exact(t);
            p.f = @(t,x) F(t)+J*(x-exact(t));
            p.jacobian = @(t,x) J;
            p.solve_kind = 'affine';
        end
        p.u0 = p.exact(p.a);
    otherwise
        error('MARBF:Problem','Unknown manuscript problem: %s',name);
end
end
