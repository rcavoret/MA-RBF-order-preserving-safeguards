function data = formal_stability_data()
%FORMAL_STABILITY_DATA Figure 1: FORMAL targets, not safeguarded-map linearizations.
% Meshes are inherited from the upstream plotting driver. All disconnected
% components are retained; no assumption that AB stability sets are nested.
data.euler_x = linspace(-3,3,400); data.euler_y = data.euler_x;
[X,Y] = meshgrid(data.euler_x,data.euler_y); Z = X+1i*Y;
data.euler_classical = 1+Z;
data.euler_MQ = 1+Z+Z.^2/2+Z.^3/2;
data.euler_MA = 1+Z+Z.^2/2-Z.^3/2;
data.ab_x = linspace(-3,3,800); data.ab_y = data.ab_x;
[X,Y] = meshgrid(data.ab_x,data.ab_y); Z = X+1i*Y;
% r^2+b(z)r+c(z). Explicit quadratic roots give the complete root condition.
b = -(1+3*Z/2-7*Z.^3/24); c = Z/2-17*Z.^3/24;
data.ab_MQ_radius = max(abs((-b+sqrt(b.^2-4*c))/2),abs((-b-sqrt(b.^2-4*c))/2));
b = -(1+3*Z/2+41*Z.^3/24); c = Z/2+31*Z.^3/24;
data.ab_MA_radius = max(abs((-b+sqrt(b.^2-4*c))/2),abs((-b-sqrt(b.^2-4*c))/2));
data.am_x = linspace(-6,6,600); data.am_y = data.am_x;
[X,Y] = meshgrid(data.am_x,data.am_y); Z = X+1i*Y;
data.am_classical_stable = abs((1+Z/2)./(1-Z/2)) <= 1;
data.am_RBF_stable = abs((1+Z/2-Z.^3/24)./(1-Z/2+Z.^3/24)) <= 1;
end
