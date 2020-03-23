function [ind_best_hairpin] = select_best_hairpin_20191017(cut,...
    MFE,n_hairpins,n_bp,loop_size,distance_stem_U,fraction_in_stem)

% number of sequences, and lengths folded per sequence
n_seq = size(MFE,1);
l_per_seq = size(MFE,2);

% ranking order for properties: 1 for low, 0 for high.
rank_order = [1 1 0 1];

% initialization
% bool_cut = NaN(n_seq,1);
ind_best_hairpin = NaN(n_seq,1);

for i = 1:n_seq
    
    inds = 1:l_per_seq;
    
    % the cut on terminator properties
%     bool = (n_bp(i,:)>=cut.bp_low) & (n_bp(i,:)<=cut.bp_high) &...
%         (loop_size(i,:)>=cut.loop_low) & (loop_size(i,:)<=cut.loop_high) &...
%         ((distance_stem_3pr(i,:)-consecutive_Us(i))<=cut.distance_stem_end_thresh) & ...
%         (fraction_in_stem(i,:)>=cut.frac_low) & (MFE(i,:)<=-cut.MFE_cut) & (n_hairpins(i,:)==1);
    
      bool = (n_bp(i,:)>=cut.bp_low) & (n_bp(i,:)<=cut.bp_high) &...
        (loop_size(i,:)>=cut.loop_low) & (loop_size(i,:)<=cut.loop_high) &...
        (distance_stem_U(i,:)<=cut.distance_stem_end_thresh) & ...
        (n_hairpins(i,:)==1);
    

    if sum(bool)==0     % no hairpin candidate passes the cut, still report the best hairpin within subset with a hairpin.
        bool_cut(i) = 0;
        inds = inds(n_hairpins(i,:)==1);
    else                % restrict the analysis to hairpins passing the cut if they do.
        bool_cut(i) = 1;
        inds = inds(bool);
    end
    
    if isempty(inds)
        ind_best_hairpin(i) = 1;    % if no hairpin for any, just report the first hairpin. 
    else
        
        % scoring properties
        prop = NaN(length(rank_order),length(inds));
        prop(1,:) = squeeze(MFE(i,inds));
        prop(2,:) = squeeze(MFE(i,inds)./n_bp(i,inds));
        prop(3,:) = squeeze(fraction_in_stem(i,inds));
        prop(4,:) = squeeze(distance_stem_U(i,inds));

        % ranking the sequences on these properties
        rank_prop = NaN(length(rank_order),length(inds));
        for j = 1:size(prop,1)
            [~,~,rank_prop(j,:)] = unique((-1)^(rank_order(j)+1)*prop(j,:));
        end
        
        % counting the top rank per sequence
        first_rank = NaN(length(inds),1);
        second_rank = NaN(length(inds),1);
        for j = 1:length(inds)
            first_rank(j) = sum(rank_prop(:,j)==1);
            second_rank(j) = sum(rank_prop(:,j)==2);
        end
        
        % reporting the first (shortest folded) index with the largest top number of top rank.
        sub_inds_first = first_rank==max(first_rank);
        inds_best = inds(sub_inds_first);
        
        % if there is a tie, report the hairpin with the maximum number of
        % second rank. If yet another tie, report lowest L hairpin. 
        if length(inds_best)>1
            max_second_rank = max(second_rank(sub_inds_first));
            ind_best_hairpin(i) = inds_best(find(second_rank(sub_inds_first)==max_second_rank,1,'first'));
        else
            ind_best_hairpin(i) = inds_best;
        end
        
    end
end
end
