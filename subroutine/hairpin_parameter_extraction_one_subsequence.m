function hairpin_param = hairpin_parameter_extraction_one_subsequence(mountain_var_subsequence)


[pks,~] = findpeaks(mountain_var_subsequence);

if isempty(pks)
    hairpin_param.n_base_pairs = 0;
    hairpin_param.n_hairpins = 0;
    hairpin_param.loop_ind = zeros(length(mountain_var_subsequence),1);
    hairpin_param.stem_3pr_ind = zeros(length(mountain_var_subsequence),1);
    hairpin_param.stem_5pr_ind = zeros(length(mountain_var_subsequence),1);
    hairpin_param.bulges_3pr = [];
    hairpin_param.bulges_5pr = [];
    hairpin_param.mountain_var = mountain_var_subsequence;
    
else
    
    % find 3' most peak, corresponding to the 3' most hairpin. 
    value_3pr_most_peak = pks(end);
    
    % loop. Need to be careful about how the mountain variable is defined. 
    loop_var = mountain_var_subsequence==value_3pr_most_peak;
    right_most_loop = find(loop_var==1,1,'last');
    left_most_loop = find(loop_var==1,1,'first')+1;
    loop_var = zeros(length(mountain_var_subsequence),1);
    loop_var(left_most_loop:right_most_loop) = 1;
    
    % 3' stem. Again, need to be careful about how the mountain variable is
    % defined and adjust for it. 
    pos_ind = 1:length(mountain_var_subsequence);
    if size(pos_ind,1) ~= size(mountain_var_subsequence,1)
       pos_ind = pos_ind'; 
    end
    stem_3pr_ind = pos_ind>right_most_loop & mountain_var_subsequence>0;
    last_ind_stem_3pr = find(stem_3pr_ind==1,1,'last');
    stem_3pr_ind(last_ind_stem_3pr+1) = 1;
    if size(stem_3pr_ind,1)>1
        stem_3pr_ind = stem_3pr_ind';
    end
    stem_3pr = mountain_var_subsequence(stem_3pr_ind);
    
    % 5' stem
    pos_ind = find(mountain_var_subsequence==0 & pos_ind<left_most_loop,1,'last')+1:(left_most_loop-1);
    if isempty(pos_ind) && (mountain_var_subsequence(1)==1)
       pos_ind = 1:left_most_loop-1; 
    end
    stem_5pr_ind = false(length(mountain_var_subsequence),1)';
    stem_5pr_ind(pos_ind) = true;
    stem_5pr = mountain_var_subsequence(stem_5pr_ind);
    
    % output
    hairpin_param.n_base_pairs = value_3pr_most_peak;
    hairpin_param.n_hairpins = length(pks);
    hairpin_param.loop_ind = loop_var';
    hairpin_param.stem_3pr_ind = stem_3pr_ind;
    hairpin_param.stem_5pr_ind = stem_5pr_ind;
    hairpin_param.bulges_3pr = select_bulges_moutain_plots(stem_3pr);
    hairpin_param.bulges_5pr = select_bulges_moutain_plots(stem_5pr);
    hairpin_param.mountain_var = mountain_var_subsequence';
    
end






