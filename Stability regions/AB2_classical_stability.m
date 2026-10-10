function h = AB2_classical_stability()
%
% Usage:   AB2_classical_stability
% Purpose: the function plots the stability region of the 
%          classical two-step Adams-Bashforth method
% Input:   none 
% Output:  optional graphics handle h; the function 
%          produces a contour plot of the stability region 
%          in the complex plane, showing where the method is 
%          stable
%
theta = linspace(0, 2*pi, 1000);
r = exp(1i*theta);
z = (r.^2 - r) ./ (1.5*r - 0.5);
h = plot(real(z), imag(z),'Color', 'b', 'LineWidth', 2);
grid on;
axis equal;
xlabel('Re(h\lambda)');
ylabel('Im(h\lambda)');

