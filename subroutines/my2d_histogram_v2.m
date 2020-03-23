function values =  my2d_histogram_v2(X,Y,edges_x,edges_y)

edges{1} = edges_x;
edges{2} = edges_y;

% bins the data
values = hist3([X(:) Y(:)],'Edges',edges);
values = values(1:end-1,1:end-1);

[X_plot, Y_plot] = meshgrid(edges{1}(1:end-1),edges{2}(1:end-1));

% temp = log10(values);
temp = values;

temp(abs(temp)==Inf)=0;
hold on
surf(X_plot',Y_plot',temp), shading flat
view(2)
% set(gca,'XScale','log')
% set(gca,'YScale','log')
set(gca,'XLim', [min(edges{1}(:)) max(edges{1}(1:end-1))])
set(gca,'YLim', [min(edges{2}(:)) max(edges{2}(1:end-1))])
colormap('gray')

caxis([0 max(max(values))])
