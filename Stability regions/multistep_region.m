function p = multistep_region(x, y, coeffs, face_color, edge_color, alpha_val, stable)
%MULTISTEP_REGION Filled root-condition region, preserving the original interface.
% An optional logical mask skips scalar ROOTS calls. Closed contour polygons
% replace the unsupported assignment of FaceAlpha to a contour-group object.
% POLYSHAPE (base MATLAB, R2017b+) handles multiple components and holes.
if nargin < 7
    [X,Y] = meshgrid(x,y); Z = X+1i*Y; stable = false(size(Z));
    for j = 1:numel(Z)
        r = roots(coeffs(Z(j)));
        stable(j) = all(abs(r) <= 1);
    end
end
% A one-cell exterior closes components that meet the plotting window.
dx = x(2)-x(1); dy = y(2)-y(1);
xp = [x(1)-dx,x,x(end)+dx]; yp = [y(1)-dy,y,y(end)+dy];
S = zeros(numel(yp),numel(xp)); S(2:end-1,2:end-1) = double(stable);
C = contourc(xp,yp,S,[0.5 0.5]); xx = []; yy = []; j = 1;
while j < size(C,2)
    count = C(2,j); nodes = C(:,j+(1:count));
    xx = [xx,nodes(1,:),NaN]; yy = [yy,nodes(2,:),NaN]; %#ok<AGROW>
    j = j+count+1;
end
if isempty(xx)
    p = patch(NaN,NaN,face_color,'EdgeColor',edge_color,'FaceAlpha',alpha_val);
else
    region = polyshape(xx,yy,'Simplify',true);
    p = plot(region,'FaceColor',face_color,'EdgeColor',edge_color, ...
        'FaceAlpha',alpha_val,'LineWidth',2);
end
axis square; grid on;
xlabel('\Re(h\lambda)'); ylabel('\Im(h\lambda)');
end
