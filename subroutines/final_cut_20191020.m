function bool_cut = final_cut_20191020(cuts,...
    loop_size_Us,n_bp_Us,MFE_Us,fraction_in_stem_Us,...
    n_hairpin_Us,distance_stem_U)

bool_cut = [];
for i = 1:length(cuts.bp_low)
    
    bool_cut{i} = (n_bp_Us>=cuts.bp_low(i)) & (n_bp_Us<=cuts.bp_high(i)) &...
        (loop_size_Us>=cuts.loop_low(i)) & (loop_size_Us<=cuts.loop_high(i)) &...
        distance_stem_U<=cuts.distance_stem_end_thresh(i) & ...
        (fraction_in_stem_Us>=cuts.frac_low(i)) & ...
        (MFE_Us<=-cuts.MFE(i)) & ...
        (n_hairpin_Us==1);
    
end

% for i = 1:length(n_bp_Us)
%     
%     MFE_f_stem_bool = false(length(n_bp_Us),1);
%     for i = 1:length(dG_cut)
%         MFE_f_stem_bool = MFE_f_stem_bool | ...
%             ( (MFE_Us<=-dG_cut(i)) & fraction_in_stem_Us>=fraction_thresholds(i));
%     end
%     
%     
%     bool_cut = MFE_f_stem_bool & (n_bp>=cut.bp_low) & (n_bp<=cut.bp_high) &...
%         (loop_size>=cut.loop_low) & (loop_size_Us<=cut.loop_high) &...
%         distance_stem_U<=cut.distance_stem_end_thresh & ...
%         (n_hairpin_Us==1);
% end


%   bool = (n_bp(i)>=cut.bp_low) & (n_bp(i)<=cut.bp_high) &...
%         (loop_size(i)>=cut.loop_low) & (loop_size<=cut.loop_high) &...
%         dist_stem_U(i)<=cut.distance_stem_end_thresh & ...
%         (fraction_in_stem>=cut.frac_low) & (fraction_in_stem<cut.frac_high) & ...
%         (MFE<=-cut.MFE_cut) & ...
%         (n_hairpin==1);