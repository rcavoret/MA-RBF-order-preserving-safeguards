function h = multistep_boundary(x, y, coeffs, color, stable)
%MULTISTEP_BOUNDARY Original root-condition boundary with optional computed mask.
% The optional mask avoids repeated scalar ROOTS calls for Figure 1's quadratics.
% Without it, the original coefficient-function interface remains available.
[X,Y] = meshgrid(x,y);
if nargin < 5
    Z = X+1i*Y; stable = false(size(Z));
    for j = 1:numel(Z)
        r = roots(coeffs(Z(j)));
        stable(j) = all(abs(r) <= 1);
    end
end
[~,h] = contour(X,Y,double(stable),[0.5 0.5],color,'LineWidth',2);
axis equal; grid on;
xlabel('\Re(h\lambda)'); ylabel('\Im(h\lambda)');
end
