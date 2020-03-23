function [dG_cut,MFE,n_bp,loop_size,fraction_in_stem,n_hairpin] = plot_hairpin_properties_summary_random(...
    all_hairpin_param_f,all_hairpin_param_r,...
    terminator_properties_f,terminator_properties_r, cut,...
    species,short_species_name,GC_content,...
    now_str,ind_species,n_hairpin_plot_f,n_hairpin_plot_r,f_pass_random)


% unpack the data
field_names = fieldnames(all_hairpin_param_f);
for i = 1:length(field_names)
    eval(sprintf('%s = [all_hairpin_param_f(1).%s all_hairpin_param_r(1).%s];',...
        field_names{i},field_names{i},field_names{i}));
end

% loop_size = [all_hairpin_param_f(1).loop_size all_hairpin_param_r(1).loop_size];
% n_bp = all_hairpin_param(i).n_bp;
% fraction_in_stem = all_hairpin_param(i).fraction_in_stem;
% MFE_hairpin = all_hairpin_param(i).MFE;
% distance_stem_U = all_hairpin_param(i).distance_stem_U;


% unpack the hairpin parameter cut variable
field_names = fieldnames(cut);
for i = 1:length(field_names)
    eval(sprintf('%s = cut.%s;',field_names{i},field_names{i}));
end


% restricting to regions with a single hairpin
bool_f = n_hairpin_plot_f==1;
bool_r = n_hairpin_plot_r==1;
n_hairpin = [n_hairpin_plot_f n_hairpin_plot_r];
bool = [bool_f'; bool_r']';

try
dG_cut_all = 0:0.1:20;
fraction_pass = [];
for i = 1:length(dG_cut_all)
    global_bool = -MFE>dG_cut_all(i) &...
        n_bp<=bp_high & ...
        n_bp>=bp_low & ...
        loop_size <= loop_high & ...
        loop_size >= loop_low & ...
        fraction_in_stem >= frac_low & bool;
    fraction_pass(i) = sum(global_bool)/length(global_bool);
end
catch
   bla 
end

ind = find(abs(fraction_pass-f_pass_random)==min(abs(fraction_pass-f_pass_random)),1,'first');
dG_cut = dG_cut_all(ind);


hfig = figure;



subplot(6,10,[7:10 17:20])
hold on
stairs(dG_cut_all,fraction_pass,'-k')
plot(dG_cut*[1 1],[1E-3 1],'--r');
default_plot(0)
set(gca,'YScale','log','YLim',[1E-3 1])
xlabel('\DeltaG cut')
ylabel('Fraction pass');



% printing summary statistics
subplot(6,10,1)
hold on
set(gca,'XLim',[0 1],'YLim',[0 1]);
fraction_overall_pass = (length(terminator_properties_f) + length(terminator_properties_r))/(sum(bool_f) + sum(bool_r));
summary_stats = [sprintf('species: %s, %s\n',species,short_species_name)];
summary_stats = [summary_stats sprintf('GC cont.: %.2f\n',GC_content)];
summary_stats = [summary_stats sprintf('Fraction pass cuts: %.3f  (For: %d/%d, Rev: %d/%d)',...
    fraction_overall_pass,...
    length(terminator_properties_f),sum(bool_f),...
    length(terminator_properties_r),sum(bool_r))];
text(0,1,summary_stats,'HorizontalAlignment','left','VerticalAlignment','cap','Interpreter','none','FontSize',14);
axis off



% 2D distribution of n_bp and loop size
subplot(6,10,[34 35 44 45]+10)   
hold on
n_bp_edges = 0:20;
loop_size_edges = 0:25;
binned_values = my2d_histogram_v2(loop_size(bool),n_bp(bool),loop_size_edges,n_bp_edges);
max_z = max(max(binned_values));
plot3([loop_low loop_low loop_high loop_high loop_low],...
    [bp_low bp_high bp_high bp_low bp_low],max_z*1.1*[1 1 1 1 1],'--r','LineWidth',2);
if max_z>0
    set(gca,'ZLim',[0 max_z*1.1]);
end
xlabel('Loop')
ylabel('N_{bp}')
set(gca,'FontSize',16)
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');

% 1D distribution n_bp
[n,edges,~] = histcounts(n_bp(bool),n_bp_edges);
subplot(6,10,[36 46]+10)
hold on
stairs([0 n/sum(n)],edges,'-k')
plot([0 1 NaN 0 1],[bp_low bp_low NaN bp_high bp_high],'--r','LineWidth',2);
xlabel('p(N_{bp})')
set(gca,'XLim',[0 1.1*max(n/sum(n))],'YLim',ylims)
default_plot(0)
set(gca,'YTickLabels','');

% 1D distribution loop size
[n,edges,~] = histcounts(loop_size(bool),loop_size_edges);
subplot(6,10,[24 25]+10)
hold on
stairs(edges,[0 n/sum(n)],'-k')
plot([loop_low loop_low NaN loop_high loop_high],[0 1 NaN 0 1],'--r','LineWidth',2);
ylabel('p(Loop)')
set(gca,'YLim',[0 1.1*max(n/sum(n))],'XLim',xlims)
default_plot(0)
set(gca,'XTickLabels','');



% 2D distribution of \Delta G and f_stem
subplot(6,10,[38 39 48 49]+10)
hold on
fraction_edges = 0:0.1:1.2;
MFE_edges = 0:21;
binned_values =my2d_histogram_v2(fraction_in_stem(bool),-MFE(bool),fraction_edges,MFE_edges);
max_z = max(max(binned_values));
plot3([frac_low frac_low 1.1],...
    [21 dG_cut dG_cut],1.1*max_z*[1 1 1],'--r','LineWidth',2);
xlabel('Fraction paired in stem')
ylabel('-\DeltaG_{hairpin}')
if max_z>0
    set(gca,'XLim',[0 1.1],'ZLim',[0 max_z*1.1])
end
default_plot(0);
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');


% 1D distribution MFE
[n,edges,~] = histcounts(-MFE(bool),MFE_edges);
subplot(6,10,[40 50]+10)
hold on
stairs([0 n/sum(n)],edges,'-k')
plot([0 1],[dG_cut dG_cut],'--r','LineWidth',2);
xlabel('p(-\DeltaG)')
if max(n)>0
    set(gca,'XLim',[0 1.1*max(n/sum(n))],'YLim',ylims)
end
default_plot(0)
set(gca,'YTickLabels','');


% 1D distribution f_stem
[n,edges,~] = histcounts(fraction_in_stem(bool),fraction_edges);
subplot(6,10,[28 29]+10)
hold on
stairs(edges,[0 n/sum(n)],'-k')
plot([frac_low frac_low],[0 1],'--r','LineWidth',2);
ylabel('p(f_{stem})')
set(gca,'XLim',[0 1.1]);
set(gca,'YLim',[0 1.1*max(n/sum(n))],'XLim',xlims)
default_plot(0)    
set(gca,'XTickLabels','');


set(gcf,'Position',[ 200         127        1443         971]);


saveas(hfig,sprintf('random_%d_%s_%s_%s.png',ind_species,species,short_species_name,now_str));
close(hfig);

    
end





