function h_fig = plotting_Urich_random_hairpins_20191024(GCF_name,short_name,...
    GC_content,...
    U_rich_hairpins,random_hairpins,cut,...
    bool_cuts,stop_to_stem_distance,...
    MFE_thresholds,fraction_in_rand,fraction_in_U,dG_cuts)
    

%% unpacking variables

field_names = fieldnames(U_rich_hairpins);
for l = 1:length(field_names)
    str_oi = sprintf('%s_Us = U_rich_hairpins.%s;',...
        field_names{l},field_names{l});
    eval(str_oi);
end
bool_U_n_hairpins = U_rich_hairpins.n_hairpins==1;

field_names = fieldnames(random_hairpins);
for l = 1:length(field_names)
    str_oi = sprintf('%s_rand = random_hairpins.%s;',...
        field_names{l},field_names{l});
    eval(str_oi);
end
bool_rand_n_hairpins = random_hairpins.n_hairpins==1;

field_names = fieldnames(cut);
for i = 1:length(field_names)
    eval(sprintf('%s = cut.%s;',field_names{i},field_names{i}));
end


h_fig = figure;

n_rows = 7;
n_cols = 18;
shift1 = 2*n_cols;
shift = 4*n_cols;


% distribution of # of poly-U regions. 
subplot(n_rows,n_cols,[2*n_cols+[1 2]]+shift1)

hold on
[Y,X] = mycdfcalc(consecutive_Us_Us);
stairs(X,Y*length(consecutive_Us_Us),'-k')
default_plot(0);
xlabel('# consecutive Us');
ylabel('cumul. # positions');
y_max = max([3500 length(consecutive_Us_Us)*1.1]);
set(gca,'XLim',[3.5 12],'YLim',[0 y_max]); %length(consecutive_Us)*1.1]);


% distribution of distance between U tract and stem. 
subplot(n_rows,n_cols,[4*n_cols+[1 2]]+shift1)

hold on
[Y,X] = mycdfcalc(distance_stem_U_Us);
stairs(X,Y,'-k')
plot(distance_stem_end_thresh(1)*[1 1],[0 1],'--b','LineWidth',2)
default_plot(0);
xlabel('Dist. 3'' stem to U tract start')
set(gca,'YLim',[0 1],'XLim',[0 15])
ylabel('CDF')



%% loop and bp: U rich
% 2D distribution of n_bp and loop size
subplot(n_rows,n_cols,[3*n_cols+4 3*n_cols+5 4*n_cols+4 4*n_cols+5]+shift1)   

hold on
n_bp_edges = 0:20;
loop_size_edges = 0:25;
binned_values = my2d_histogram_v2(loop_size_Us(bool_U_n_hairpins),n_bp_Us(bool_U_n_hairpins),loop_size_edges,n_bp_edges);
max_z = max(max(binned_values));
plot3([loop_low(1) loop_low(1) loop_high(1) loop_high(1) loop_low(1)],...
    [bp_low(1) bp_high(1) bp_high(1) bp_low(1) bp_low(1)],max_z*1.1*[1 1 1 1 1],'--b','LineWidth',2);
if max_z>0
    set(gca,'ZLim',[0 max_z*1.1]);
end
xlabel('Loop size')
ylabel('N_{bp}')
set(gca,'FontSize',16)
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');

% 1D distribution n_bp
[n,edges,~] = histcounts(n_bp_Us(bool_U_n_hairpins),n_bp_edges);
subplot(n_rows,n_cols,[3*n_cols+6 4*n_cols+6]+shift1)
hold on
stairs([0 n/sum(n)],edges,'-k')
plot([0 1 NaN 0 1],[bp_low(1) bp_low(1) NaN bp_high(1) bp_high(1)],'--b','LineWidth',2);
xlabel('p(N_{bp})')
set(gca,'XLim',[0 1.1*max(n/sum(n))],'YLim',ylims)
default_plot(0)
set(gca,'YTickLabels','');

% 1D distribution loop size
[n,edges,~] = histcounts(loop_size_Us(bool_U_n_hairpins),loop_size_edges);
subplot(n_rows,n_cols,[2*n_cols+4 2*n_cols+5]+shift1)
hold on
stairs(edges,[0 n/sum(n)],'-k')
plot([loop_low(1) loop_low(1) NaN loop_high(1) loop_high(1)],[0 1 NaN 0 1],'--b','LineWidth',2);
ylabel('p(Loop)')
set(gca,'YLim',[0 1.1*max(n/sum(n))],'XLim',xlims)
default_plot(0)
set(gca,'XTickLabels','');


%% loop and bp: random

% 2D distribution of n_bp and loop size
subplot(n_rows,n_cols,[3*n_cols+4 3*n_cols+5 4*n_cols+4 4*n_cols+5]+shift1-shift)   
hold on
n_bp_edges = 0:20;
loop_size_edges = 0:25;
binned_values = my2d_histogram_v2(loop_size_rand(bool_rand_n_hairpins),n_bp_rand(bool_rand_n_hairpins),loop_size_edges,n_bp_edges);
max_z = max(max(binned_values));
plot3([loop_low(1) loop_low(1) loop_high(1) loop_high(1) loop_low(1)],...
    [bp_low(1) bp_high(1) bp_high(1) bp_low(1) bp_low(1)],max_z*1.1*[1 1 1 1 1],'--b','LineWidth',2);
if max_z>0
    set(gca,'ZLim',[0 max_z*1.1]);
end
xlabel('Loop size')
ylabel('N_{bp}')
set(gca,'FontSize',16)
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');

% 1D distribution n_bp
[n,edges,~] = histcounts(n_bp_rand(bool_rand_n_hairpins),n_bp_edges);
subplot(n_rows,n_cols,[3*n_cols+6 4*n_cols+6]+shift1-shift)
hold on
stairs([0 n/sum(n)],edges,'-k')
plot([0 1 NaN 0 1],[bp_low(1) bp_low(1) NaN bp_high(1) bp_high(1)],'--b','LineWidth',2);
xlabel('p(N_{bp})')
set(gca,'XLim',[0 1.1*max(n/sum(n))],'YLim',ylims)
default_plot(0)
set(gca,'YTickLabels','');

% 1D distribution loop size
[n,edges,~] = histcounts(loop_size_rand(bool_rand_n_hairpins),loop_size_edges);
subplot(n_rows,n_cols,[2*n_cols+4 2*n_cols+5]+shift1-shift)
hold on
stairs(edges,[0 n/sum(n)],'-k')
plot([loop_low(1) loop_low(1) NaN loop_high(1) loop_high(1)],[0 1 NaN 0 1],'--b','LineWidth',2);
ylabel('p(Loop)')
set(gca,'YLim',[0 1.1*max(n/sum(n))],'XLim',xlims)
default_plot(0)
set(gca,'XTickLabels','');



%% dG and f_stem: U rich

load('categorical_colormap_v2.mat');
cat_col2 = cat_col2(5:end,:);

df = 0.05;
d_dG = 0.5;

fraction_edges = 0:df:(1+2*df);
MFE_edges = 0:d_dG:35;

% 2D distribution of \Delta G and f_stem
subplot(n_rows,n_cols,[3*n_cols+8 3*n_cols+9 4*n_cols+8 4*n_cols+9]+shift1)
hold on

binned_values =my2d_histogram_v2(fraction_in_stem_Us(bool_U_n_hairpins),-MFE_Us(bool_U_n_hairpins),fraction_edges,MFE_edges);
max_z_dG_U = max(max(binned_values));
for i = 1:length(cut.frac_low)
    plot3([cut.frac_low(i) cut.frac_low(i) 1+df],...
        [max(MFE_edges) cut.MFE(i) cut.MFE(i)],1.1*max_z_dG_U*ones(3,1),'--','LineWidth',2,'Color',cat_col2(i,:));
end
xlabel('Fraction paired in stem')
ylabel('-\DeltaG_{hairpin}')
if max_z_dG_U>0
    set(gca,'XLim',[0 max(fraction_edges)-df],'ZLim',[0 max_z_dG_U*1.1])
end
default_plot(0);
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');


% 1D distribution MFE
[n,edges,~] = histcounts(-MFE_Us(bool_U_n_hairpins),MFE_edges);
subplot(n_rows,n_cols,[3*n_cols+10 4*n_cols+10]+shift1)
hold on
stairs([0 n/sum(n)],edges,'-k')
for i = 1:length(cut.frac_low)
    plot([0 1],[cut.MFE(i) cut.MFE(i)],'--','LineWidth',2,'Color',cat_col2(i,:));
end
xlabel('p(-\DeltaG)')
if max(n)>0
    set(gca,'XLim',[0 1.1*max(n/sum(n))],'YLim',ylims)
end
default_plot(0)
set(gca,'YTickLabels','');


% 1D distribution f_stem
[n,edges,~] = histcounts(fraction_in_stem_Us(bool_U_n_hairpins),fraction_edges);
subplot(n_rows,n_cols,[2*n_cols+8 2*n_cols+9]+shift1)
hold on
stairs(edges-df,[0 n/sum(n)],'-k')
for i = 1:length(cut.frac_low)
    plot([cut.frac_low(i) cut.frac_low(i)],[0 1],'--','LineWidth',2,'Color',cat_col2(i,:));
end
ylabel('p(f_{stem})')
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[0 1.1*max(n/sum(n))],'XLim',xlims)
default_plot(0)    
set(gca,'XTickLabels','');




%% dG and f_stem random


% 2D distribution of \Delta G and f_stem
subplot(n_rows,n_cols,[3*n_cols+8 3*n_cols+9 4*n_cols+8 4*n_cols+9]+shift1-shift)
hold on
binned_values =my2d_histogram_v2(fraction_in_stem_rand(bool_rand_n_hairpins),-MFE_rand(bool_rand_n_hairpins),fraction_edges,MFE_edges);
max_z_dG_rand = max(max(binned_values));
for i = 1:length(cut.frac_low)
    plot3([cut.frac_low(i) cut.frac_low(i) 1+df],...
        [max(MFE_edges) cut.MFE(i) cut.MFE(i)],1.1*max_z_dG_rand*ones(3,1),'--','LineWidth',2,'Color',cat_col2(i,:));
end
xlabel('Fraction paired in stem')
ylabel('-\DeltaG_{hairpin}')
if max_z_dG_rand>0
    set(gca,'XLim',[0 max(fraction_edges)-df],'ZLim',[0 max_z_dG_rand*1.1])
end
default_plot(0);
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');


% 1D distribution MFE
[n,edges,~] = histcounts(-MFE_rand(bool_rand_n_hairpins),MFE_edges);
subplot(n_rows,n_cols,[3*n_cols+10 4*n_cols+10]+shift1-shift)
hold on
stairs([0 n/sum(n)],edges,'-k')
for i = 1:length(cut.frac_low)
    plot([0 1],[cut.MFE(i) cut.MFE(i)],'--','LineWidth',2,'Color',cat_col2(i,:));
end
xlabel('p(-\DeltaG)')
if max(n)>0
    set(gca,'XLim',[0 1.1*max(n/sum(n))],'YLim',ylims)
end
default_plot(0)
set(gca,'YTickLabels','');


% 1D distribution f_stem
[n,edges,~] = histcounts(fraction_in_stem_rand(bool_rand_n_hairpins),fraction_edges);
subplot(n_rows,n_cols,[2*n_cols+8 2*n_cols+9]+shift1-shift)
hold on
stairs(edges-df,[0 n/sum(n)],'-k')
for i = 1:length(cut.frac_low)
    plot([cut.frac_low(i) cut.frac_low(i)],[0 1],'--','LineWidth',2,'Color',cat_col2(i,:));
end
ylabel('p(f_{stem})')
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[0 1.1*max(n/sum(n))],'XLim',xlims)
default_plot(0)    
set(gca,'XTickLabels','');


cmap = flipud(gray);
colormap(cmap);


%% stop to stem distnace

subplot(n_rows,n_cols,[(n_cols-3):n_cols n_cols+((n_cols-3):n_cols) 2*n_cols+(n_cols-3:n_cols) ])
hold on
for i = 1:length(cut.frac_low)
    [Y,X] = mycdfcalc(stop_to_stem_distance{i}(~isnan(stop_to_stem_distance{i})));
    stairs(X,Y,'-','Color',cat_col2(i,:),'LineWidth',1);
end
full_stop_stem_distance = NaN(size(stop_to_stem_distance{1}));
for i = 1:length(stop_to_stem_distance{1})
    if ~isnan(stop_to_stem_distance{1}(i))
        full_stop_stem_distance(i) = stop_to_stem_distance{1}(i);
    elseif ~isnan(stop_to_stem_distance{2}(i))
        full_stop_stem_distance(i) = stop_to_stem_distance{2}(i);
    end
end
[Y,X] = mycdfcalc(full_stop_stem_distance(~isnan(full_stop_stem_distance)));
stairs(X,Y,'-k','LineWidth',1);
    

    
xlabel('Stop to stem distance (nt)')
ylabel('CDF')
default_plot(0)
set(gca,'XLim',[-20 150],'YLim',[0 1]);

% venn diagram of the overlap between the targets

n_single = [];
legend_text = [];
for i = 1:length(cut.frac_low)
    n_single(i) = sum(~isnan(stop_to_stem_distance{i}));
    legend_text{i} = sprintf('cut %d: n_{%d}=%d, q_{10}=%d nt, q_{50}=%d nt',...
        i,i,n_single(i),...
        ceil(quantile(stop_to_stem_distance{i}(~isnan(stop_to_stem_distance{i})),[0.1 0.5])));
end

legend(legend_text,'Location','SouthEast','FontSize',10);


n_pairs = [];
pairs = [1 2];
for i = 1:size(pairs,1)
    n_pairs(i) = sum(~isnan(stop_to_stem_distance{pairs(i,1)}) & ...
        ~isnan(stop_to_stem_distance{pairs(i,2)}));
end

% n_triple = sum(~isnan(stop_to_stem_distance{1}) & ...
%         ~isnan(stop_to_stem_distance{2}) & ...
%         ~isnan(stop_to_stem_distance{3}));

% subplot(n_rows,n_cols,[4*n_cols-1 4*n_cols 5*n_cols-1 5*n_cols])
% dy_text = 0.15;
% offset_text = 0.5;
title_str = [];
% for i = 1:length(cut.frac_low)
%     title_str = [title_str sprintf('n_{%d} = %d, ',i,n_single(i))];
% %     text(0,offset_text-dy_text*i,sprintf('n_{%d} = %d\n',i,n_single(i)),'HorizontalAlignment','left','VerticalAlignment','cap','FontSize',14);
% end

for i = 1:size(pairs,1)
    title_str = [title_str sprintf('n_{%d and %d} = %d, ',pairs(i,1),pairs(i,2),n_pairs(i))];
%     text(0,offset_text-3*dy_text-dy_text*i,sprintf('n_{%d%d} = %d,\n',pairs(i,1),pairs(i,2),n_pairs(i)),'HorizontalAlignment','left','VerticalAlignment','cap','FontSize',14);
end
title_str = [title_str sprintf('n_{1 or 2} = %d.',sum(~isnan(full_stop_stem_distance)))];
% text(0,offset_text-3*dy_text-dy_text*(i+1),sprintf('n_{123} = %d\n',n_triple),'HorizontalAlignment','left','VerticalAlignment','cap','FontSize',14);
% axis off
title(title_str,'FontWeight','normal','FontSize',12)
% try
% [H,S] = venn(n_single, [n_pairs n_triple],'FaceColor',{cat_col2(1,:), cat_col2(2,:),cat_col2(3,:) });
% catch
%     bla
% end
% axis equal
% axis off
% A is a three element vector [c1 c2 c3], 
%and I is a four element vector [i12 i13 i23 i123], specifiying the 
%two-circle intersection areas i12, i13, i23, and the three-circle
%intersection i123.


%% thresholding intermediate data

frac_lims = [1E-3 1];
subplot(n_rows,n_cols,[6*n_cols-6 6*n_cols-5 7*n_cols-6 7*n_cols-5])
hold on
for j = 1:length(cut.bp_high)
    stairs(MFE_thresholds,fraction_in_U(:,j),'-','Color',cat_col2(j,:),'LineWidth',1.5);
    plot([dG_cuts(j) dG_cuts(j)],frac_lims,'--','Color',cat_col2(j,:))
end
xlabel('-\DeltaG_{hairpin} threshold')
ylabel('Fraction in Us')
default_plot(0);
set(gca,'YScale','log','XLim',[min(MFE_thresholds) max(MFE_thresholds)],'YLim',frac_lims)


subplot(n_rows,n_cols,[2*n_cols-6 2*n_cols-5 3*n_cols-6 3*n_cols-5])
hold on
for j = 1:length(cut.bp_high)
    stairs(MFE_thresholds,fraction_in_rand(:,j),'-','Color',cat_col2(j,:),'LineWidth',1.5);
    plot([dG_cuts(j) dG_cuts(j)],frac_lims,'--','Color',cat_col2(j,:))
end


xlabel('-\DeltaG_{hairpin} threshold')
ylabel('Fraction in Rand')
default_plot(0);
set(gca,'YScale','log','XLim',[min(MFE_thresholds) max(MFE_thresholds)],'YLim',frac_lims)


subplot(n_rows,n_cols,[5*n_cols-3 5*n_cols-2 5*n_cols-1 5*n_cols 6*n_cols-3 6*n_cols-2 6*n_cols-1 6*n_cols 7*n_cols-3 7*n_cols-2 7*n_cols-1 7*n_cols])
hold on
for j = 1:length(cut.bp_high)
    stairs(fraction_in_rand(:,j),fraction_in_U(:,j),'-','Color',cat_col2(j,:),'LineWidth',1);
    %         plot([dG_cut dG_cut],frac_lims,':','Color',cat_col2(j,:))
end
xlabel('Fraction in Rand')
ylabel('Fraction in Us')
default_plot(1)
set(gca,'XLim',[1E-4 1],'YLim',[1E-4 1])


%% text results 

subplot(n_rows,n_cols,1)
hold on
set(gca,'XLim',[0 1],'YLim',[0 1]);
strain_info = sprintf('Genbank ID: %s\nSpecies: %s\n',GCF_name,strrep(short_name,'_',' '));
strain_info = [strain_info sprintf('GC cont.: %.2f\n',GC_content)];
text(-2,1,strain_info,'HorizontalAlignment','left','VerticalAlignment','cap','Interpreter','none','FontSize',14);

for j = 1:length(cut.frac_low)
    try
    text(-2,0.5-0.17*j,sprintf('Fraction pass cuts %d : %.2f  (F: %d/%d, R: %d/%d).',...
        j,sum(~isnan(stop_to_stem_distance{j}))/length(bool_cuts{j}),...
        sum(~isnan(stop_to_stem_distance{j}) & strand_Us==1),sum(strand_Us==1),...
        sum(~isnan(stop_to_stem_distance{j}) & strand_Us==0),sum(strand_Us==0)),'Color',cat_col2(j,:),...
        'FontSize',14);
    catch
       bla 
    end
end
axis off


set(gcf,'Position',[  -364        1256        2547        1174]);

