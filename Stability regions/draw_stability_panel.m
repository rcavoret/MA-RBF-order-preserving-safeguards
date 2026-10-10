function draw_stability_panel(data, panel)
%DRAW_STABILITY_PANEL Plot one manuscript panel into the current axes.
switch panel
    case 1
        h1 = onestep_boundary(data.euler_x,data.euler_y,data.euler_classical,'k'); hold on
        h2 = onestep_boundary(data.euler_x,data.euler_y,data.euler_MQ,'b');
        h3 = onestep_boundary(data.euler_x,data.euler_y,data.euler_MA,'g');
        legend([h1,h2,h3],{'Euler','MQ-RBF Euler','MA-RBF Euler'},'Location','northeast');
        axis([-3,3,-3,3]);
    case 2
        h1 = AB2_classical_stability(); hold on
        h2 = multistep_boundary(data.ab_x,data.ab_y,[],'r',data.ab_MQ_radius <= 1);
        h3 = multistep_boundary(data.ab_x,data.ab_y,[],'g',data.ab_MA_radius <= 1);
        legend([h1,h2,h3],{'AB2','MQ-RBF AB2','MA-RBF AB2'},'Location','northeast');
        axis([-2,1,-1.5,1.5]);
    case 3
        classical = data.am_classical_stable; rbf = data.am_RBF_stable;
        % Precompose the upstream alpha=0.7 colours. Opaque polygons keep PDF
        % and EPS vector output portable without renderer-dependent transparency.
        pink = 0.7*[1,0.75,0.75]+0.3*[1,1,1];
        blue = 0.7*[0.6,0.8,1]+0.3*[1,1,1];
        overlap = 0.7*[0.6,0.8,1]+0.3*pink;
        hold on
        multistep_region(data.am_x,data.am_y,[],pink,'none',1,classical & ~rbf);
        multistep_region(data.am_x,data.am_y,[],blue,'none',1,rbf & ~classical);
        multistep_region(data.am_x,data.am_y,[],overlap,'none',1,rbf & classical);
        multistep_boundary(data.am_x,data.am_y,[],'r',classical);
        multistep_boundary(data.am_x,data.am_y,[],'b',rbf);
        h1 = patch(NaN,NaN,pink,'EdgeColor','r','LineWidth',2);
        h2 = patch(NaN,NaN,blue,'EdgeColor','b','LineWidth',2);
        legend([h1,h2],{'AM1','MQ/MA-RBF AM1'},'Location','northeast');
        axis([-6,6,-6,6]);
    otherwise
        error('MARBF:Panel','Panel must be 1, 2, or 3.');
end
axis square; grid on; box on;
set(gca,'FontName','Times','FontSize',10);
xlabel('Re(h\lambda)'); ylabel('Im(h\lambda)');
hold off
end
