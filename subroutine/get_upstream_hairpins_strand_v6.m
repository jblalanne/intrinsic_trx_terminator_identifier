function all_hairpin_param = ...
    get_upstream_hairpins_strand_v6(...
    genome, gene_names, start, stop, strand,...
    positions_U_f, positions_U_r, forward_bool,...
    folding_directory, lower_n, upper_n,lengths_folded)


%% concatenating the list of positions of poly-Us.

range = lower_n:upper_n;
consecutive_Us = [];
positions = [];

if forward_bool
    for i = 1:length(range)
        positions = [positions; positions_U_f{i}];
        consecutive_Us = [consecutive_Us; range(i)*ones(length(positions_U_f{i}),1)];
    end
else
    for i = 1:length(range)
        positions = [positions; positions_U_r{i}];
        consecutive_Us = [consecutive_Us; range(i)*ones(length(positions_U_r{i}),1)];
    end
end


%% folding the upstream RNA sequences

% folding parameters
constraint_bool_folding = 1;        % prevents the U tract from being folded

if forward_bool
    [MFE, ~, hairpin_param] = ...
        master_RNAfold_file_v2(genome,positions,lengths_folded,...
        constraint_bool_folding,[folding_directory '_f'],...
        consecutive_Us);
else
    positions_reverse_convention = length(genome)-positions+1;
    RC_genome = char(RCsequence(genome));
    [MFE, ~, hairpin_param] = ...
        master_RNAfold_file_v2(RC_genome, positions_reverse_convention,...
        lengths_folded,constraint_bool_folding,[folding_directory '_r'],...
        consecutive_Us);
end

% the output hairpin_param is an array that is length(positions) vs.
% length(lengths_folded) in size, containing information about the folded
% sequences.


%% getting the higher level hairpin parameters

n_bp = NaN(size(hairpin_param));
loop_size= NaN(size(hairpin_param));
distance_stem_3pr= NaN(size(hairpin_param));
fraction_in_stem= NaN(size(hairpin_param));
n_hairpins= NaN(size(hairpin_param));


for i = 1:length(positions)
    for j = 1:length(lengths_folded)
        [n_bp(i,j), loop_size(i,j), distance_stem_3pr(i,j), fraction_in_stem(i,j), n_hairpins(i,j)] = ...
            analyse_hairpin(hairpin_param(i,j).mountain_var);
    end
end


%% Report all hairpins


% [ind_best_hairpin, bool_cut] = select_best_hairpin(MFE,n_hairpins,n_bp,loop_size,distance_stem_3pr,fraction_in_stem,cut,consecutive_Us);

% saving the properties from the best hairping from folded sequences for
% downstream QC purposes.

% n_bp2 = NaN(length(positions),1);
% loop_size2 = NaN(length(positions),1);
% fraction_in_stem2 = NaN(length(positions),1);
% MFE2 = NaN(length(positions),1);
% n_hairpins2 = NaN(length(positions),1);

distance_stem_U = NaN(length(positions),1);
mountain_var = [];
gene_upstream = [];
gene_downstream = [];
upstream_sequence = [];

for i = 1:length(positions)
    %     n_bp2(i) = n_bp(i,ind_best_hairpin(i));
    %     loop_size2(i) = loop_size(i,ind_best_hairpin(i));
    %     fraction_in_stem2(i) = fraction_in_stem(i,ind_best_hairpin(i));
    %     MFE2(i) = MFE(i,ind_best_hairpin(i));
    %     n_hairpins2(i) = n_hairpins(i,ind_best_hairpin(i));

    for j = 1:length(lengths_folded)
        distance_stem_U(i,j) = distance_stem_3pr(i,j)-consecutive_Us(i);
        
        
        if forward_bool
            min_pos = max([1 (positions(i)-lengths_folded(j)+1)]);
            upstream_sequence{i,j} = genome(min_pos:positions(i));
        else
            max_pos = min([(positions(i)+lengths_folded(j)-1) length(genome)]);
            upstream_sequence{i,j} = RCsequence(genome(positions(i):max_pos));
        end
        
        mountain_var(i,j,:) = hairpin_param(i,j).mountain_var;

    end
    
    [gene_upstream{i}, gene_downstream{i}] = get_upstream_gene_v2(positions(i),forward_bool,...
        gene_names, strand, start, stop);
    
%     [closest_upstream_stop, distance_upstream_stop, distance_stem_stop] = ...
%         closest_upstream_stop_codon_v4(start, stop, strand,terminator_properties);
    
end

all_hairpin_param = struct([]);
all_hairpin_param(1).n_bp = n_bp;
all_hairpin_param(1).loop_size = loop_size;
all_hairpin_param(1).fraction_in_stem = fraction_in_stem;
all_hairpin_param(1).MFE = MFE;
all_hairpin_param(1).distance_stem_U = distance_stem_U;
all_hairpin_param(1).consecutive_Us = consecutive_Us;
all_hairpin_param(1).n_hairpins = n_hairpins;
all_hairpin_param(1).gene_upstream = gene_upstream;
all_hairpin_param(1).gene_downstream = gene_downstream;
all_hairpin_param(1).upstream_sequence = upstream_sequence;
% try
%     if ~isempty(hairpin_param)
all_hairpin_param(1).mountain_var = mountain_var;
%     else
%         all_hairpin_param(1).mountain_var = [];
%     end
% catch
%     bla
% end
all_hairpin_param(1).positions = positions;




%
%
% %% collect terminator properties (for those passing the cut!)
%
% terminator_properties = [];
% counter = 1;
%
% for i = 1:length(positions)
%
%     if bool_cut(i)
%         terminator_properties(counter).position = positions(i);
% %
% %         if forward_bool
% %             terminator_properties(counter).upstream_sequence = genome((positions(i)-max(lengths_folded)+1):positions(i));
% %         else
% %             terminator_properties(counter).upstream_sequence = RCsequence(genome(positions(i):(positions(i)+max(lengths_folded)-1)));
% %         end
%
%
% %         [gene_upstream, gene_downstream] = get_upstream_gene_v2(positions(i),forward_bool,...
% %             gene_names, strand, start, stop);
%         ind = ind_best_hairpin(i);
%         terminator_properties(counter).upstream_sequence = upstream_sequence{i};
%         terminator_properties(counter).folded_seq_id = i;
%         terminator_properties(counter).gene_upstream = gene_upstream{i};
%         terminator_properties(counter).gene_downstream = gene_downstream{i};
%         terminator_properties(counter).RNA_structure = hairpin_param(i,ind).mountain_var;
%         terminator_properties(counter).n_bp = n_bp(i,ind);
%         terminator_properties(counter).loop = loop_size(i,ind);
%         terminator_properties(counter).MFE = MFE(i,ind);
%         terminator_properties(counter).longest_U = consecutive_Us(i);
%         terminator_properties(counter).fraction_stem = fraction_in_stem(i,ind);
%         terminator_properties(counter).distance_stem_3pr = distance_stem_3pr(i,ind);
%         terminator_properties(counter).strand = forward_bool;
%         terminator_properties(counter).n_hairpins = n_hairpins(i,ind);
%
%         counter = counter+1;
%     end
%
% end
%
%
% %% obtain the distance from stem to stop for candidate terminators
%
% % moving forward only if certain candidates pass the cut.
% if ~isempty(terminator_properties)
%
%     [closest_upstream_stop, distance_upstream_stop, distance_stem_stop] = ...
%         closest_upstream_stop_codon_v4(start, stop, strand,terminator_properties);
%
%     % adding the information to the terminator struct.
%     for i = 1:length(terminator_properties)
%         terminator_properties(i).closest_upstream_stop = closest_upstream_stop(i);
%         terminator_properties(i).distance_stem_stop =  distance_stem_stop(i);
%         terminator_properties(i).distance_upstream_stop =  distance_upstream_stop(i);
%     end
%
%
%     %% identifying downstream terminators in multi-terminator instances
%
%     dist_thresh = 10;
%     [bool_multi_terminators, ~, close_sets] = exclude_multiple_terminators(terminator_properties,dist_thresh);
%
%     % for nearby terminators, select the best based on properties
%     best_close_sets = select_best_terminator_v2(close_sets,terminator_properties);
%
%     % update the bool_multi_terminators
%     for i = 1:length(close_sets)
%         bool_multi_terminators(best_close_sets(i)) = 1;
%         for j = 1:length(close_sets{i})
%             if close_sets{i}(j)~=best_close_sets(i)
%                 bool_multi_terminators(close_sets{i}(j)) = 0;
%             end
%         end
%     end
%
%     % add the information to the terminator struct.
%     for i = 1:length(terminator_properties)
%         terminator_properties(i).bool_multi_terminators = bool_multi_terminators(i);
%     end
%
%     % reporting the distance between stop codons and hairpins for final
%     % list.
%     final_stem_to_stop = distance_stem_stop((bool_multi_terminators==1) & ~isnan(distance_stem_stop));
%
% else
%     final_stem_to_stop = [];
% end
%

