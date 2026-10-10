function save_publication_figure(fig, directory, stem)
%SAVE_PUBLICATION_FIGURE Vector PDF/EPS and 600-dpi PNG, with explicit page size.
set(fig,'PaperPositionMode','auto','InvertHardcopy','off','Color','w');
position = get(fig,'Position'); units = get(fig,'Units');
set(fig,'Units','centimeters'); cm = get(fig,'Position');
set(fig,'PaperUnits','centimeters','PaperSize',cm(3:4), ...
    'PaperPosition',[0,0,cm(3:4)]);
set(fig,'Units',units,'Position',position);
print(fig,fullfile(directory,[stem '.pdf']),'-dpdf','-painters');
print(fig,fullfile(directory,[stem '.eps']),'-depsc2','-painters');
print(fig,fullfile(directory,[stem '.png']),'-dpng','-r600');
fprintf('Saved figure: %s\n',fullfile(directory,[stem '.pdf']));
end
