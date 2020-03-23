function [bool_cuts,dG_cuts,cut,fraction_in_U,fraction_in_rand] = get_final_cut_20191019(MFE_thresholds,...
   U_rich_hairpins,random_hairpins,cut,plot_bool)


   
%% unpacking variables

field_names = fieldnames(U_rich_hairpins);
for l = 1:length(field_names)
    str_oi = sprintf('%s_Us = U_rich_hairpins.%s;',...
        field_names{l},field_names{l});
    eval(str_oi);
end

field_names = fieldnames(random_hairpins);
for l = 1:length(field_names)
    str_oi = sprintf('%s_rand = random_hairpins.%s;',...
        field_names{l},field_names{l});
    eval(str_oi);
end


   
%% thresholding

fraction_in_U = NaN(length(MFE_thresholds),length(cut.bp_high));
fraction_in_rand = NaN(length(MFE_thresholds),length(cut.bp_high));


for j = 1:length(MFE_thresholds)
    
    fraction_in_U(j,:) = number_threshold_v2(loop_size_Us,n_bp_Us,...
        MFE_Us,fraction_in_stem_Us,n_hairpins_Us,cut,distance_stem_U_Us,1,...
        MFE_thresholds(j))/length(loop_size_Us);
    
    fraction_in_rand(j,:) = number_threshold_v2(loop_size_rand,n_bp_rand,...
        MFE_rand,fraction_in_stem_rand,n_hairpins_rand,cut,[],0,...
        MFE_thresholds(j))/length(loop_size_rand);
    
end

dG_cuts = [];
for i = 1:length(cut.f_pass)
    if cut.f_pass(i)>0
        ind = find(fraction_in_rand(:,i)>cut.f_pass(i),1,'last');
        if ~isempty(ind)
            dG_cuts(i) = MFE_thresholds(ind);
            dG_cuts(i) = max([cut.MFE(i) dG_cuts(i)]);
        else
            dG_cuts(i) = cut.MFE(i);
        end
        cut.MFE(i) = dG_cuts(i);
    else
        dG_cuts(i) = cut.MFE(i);
    end
end


%% generating the final cut

% cut.MFE_cut = dG_cut;
bool_cuts = final_cut_20191020(cut,...
    loop_size_Us,n_bp_Us,MFE_Us,fraction_in_stem_Us,...
    n_hairpins_Us,distance_stem_U_Us);


%% plotting


load('categorical_colormap_v2.mat');
cat_col2 = cat_col2(5:end,:);


if plot_bool
    
    figure 
    n_rows = 2;
    n_cols = 4;
    frac_lims = [1E-3 1];
    
    subplot(n_rows,n_cols,1)
    hold on
    for j = 1:length(cut.bp_high)
        stairs(MFE_thresholds,fraction_in_U(:,j),'-','Color',cat_col2(j,:),'LineWidth',1.5);
        plot([dG_cuts(j) dG_cuts(j)],frac_lims,'--','Color',cat_col2(j,:))
    end
    xlabel('-\DeltaG_{hairpin} threshold')
    ylabel('Fraction in Us')
    default_plot(0);
    set(gca,'YScale','log','XLim',[min(MFE_thresholds) max(MFE_thresholds)],'YLim',frac_lims)
    
    
    subplot(n_rows,n_cols,5)
    hold on
    for j = 1:length(cut.bp_high)
        stairs(MFE_thresholds,fraction_in_rand(:,j),'-','Color',cat_col2(j,:),'LineWidth',1.5);
        plot([dG_cuts(j) dG_cuts(j)],frac_lims,'--','Color',cat_col2(j,:))
    end
    
%     stairs(MFE_thresholds,fraction_in_rand,'-k');
%     plot([dG_cut dG_cut],frac_lims,'--r','LineWidth',2)

    xlabel('-\DeltaG_{hairpin} threshold')
    ylabel('Fraction in Rand')
    default_plot(0);
    set(gca,'YScale','log','XLim',[min(MFE_thresholds) max(MFE_thresholds)],'YLim',frac_lims)
    
    
    subplot(n_rows,n_cols,[3 4 7 8])
    hold on
    for j = 1:length(cut.bp_high)
        stairs(fraction_in_rand(:,j),fraction_in_U(:,j),'-','Color',cat_col2(j,:),'LineWidth',1);
%         plot([dG_cut dG_cut],frac_lims,':','Color',cat_col2(j,:))
    end
%     stairs(fraction_in_rand,fraction_in_U,'-k')
%     plot(fraction_in_rand(ind),fraction_in_U(ind),'or','LineWidth',2)
    
    xlabel('Fraction in Rand')
    ylabel('Fraction in Us')
    default_plot(1)
    set(gca,'XLim',[1E-4 1],'YLim',[1E-4 1])
    
    set(gcf,'Position',[680   646   904   452]);
end

