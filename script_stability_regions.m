% Figure 1(a-c): formal-target Euler, AB2 and AM1 stability diagnostics.
% One run writes the COMPLETE Figure 1 and all three standalone panels.
addpath(fileparts(mfilename('fullpath')));
main();

function main()
paths = setup_package();
data = formal_stability_data();
titles = {'(a) Euler methods','(b) AB2 methods','(c) AM1 methods'};
for k = 1:3
    fig = figure('Color','w','Visible','off','Units','centimeters', ...
        'Position',[1,1,16,12.8]);
    closer = onCleanup(@() close(fig));
    axes('Parent',fig,'Position',[0.14,0.13,0.78,0.79]);
    draw_stability_panel(data,k);
    save_publication_figure(fig,paths.figures,sprintf('figure1%c','a'+k-1));
    clear closer
end
fig = figure('Color','w','Visible','off','Units','centimeters', ...
    'Position',[1,1,24,8.8]);
closer = onCleanup(@() close(fig)); %#ok<NASGU>
for k = 1:3
    axes('Parent',fig,'Position',[0.055+(k-1)*0.325,0.16,0.275,0.75]);
    draw_stability_panel(data,k);
    set(gca,'FontSize',8);
    legend_handle = legend;
    set(legend_handle,'FontSize',6.5);
    title(titles{k},'FontWeight','normal','FontSize',9);
end
save_publication_figure(fig,paths.figures,'figure1');
metadata = struct('matlab_version',version,'created',datestr(now,30), ...
    'interpretation','formal-target diagnostics, not stability regions of the active safeguarded maps');
save(fullfile(paths.data,'figure1_stability_data.mat'),'data','metadata');
end
