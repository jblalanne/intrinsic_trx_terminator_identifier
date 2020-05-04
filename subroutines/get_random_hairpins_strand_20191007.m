function [terminator_properties, all_hairpin_param,...
    positions,MFE,n_hairpins,n_bp,loop_size,fraction_in_stem] = ...
    get_random_hairpins_strand_20191007(...
    genome,...
    rand_pos_f, rand_pos_r, forward_bool,...
    folding_directory,lengths_folded,cut)


%% concatenating the list of positions of poly-Us.

positions = [];
if forward_bool
    for i = 1:length(rand_pos_f)
        positions = [positions; rand_pos_f{i}];
    end
else
    for i = 1:length(rand_pos_r)
        positions = [positions; rand_pos_r{i}];
    end
end




%% folding the upstream RNA sequences

% folding parameters
constraint_bool_folding = 1;        % prevents the U tract from being folded


consecutive_Us = zeros(size(positions));
 
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


%% select the best hairpin candidate for next steps (includes the cut on the hairpin parameters).


[ind_best_hairpin, bool_cut] = select_best_hairpin_random(MFE,n_hairpins,n_bp,loop_size,distance_stem_3pr,fraction_in_stem,cut);

% saving the properties from the best hairping from folded sequences for
% downstream QC purposes.
n_bp2 = NaN(length(positions),1);
loop_size2 = NaN(length(positions),1);
fraction_in_stem2 = NaN(length(positions),1);
MFE2 = NaN(length(positions),1);
for i = 1:length(positions)
    n_bp2(i) = n_bp(i,ind_best_hairpin(i));
    loop_size2(i) = loop_size(i,ind_best_hairpin(i));
    fraction_in_stem2(i) = fraction_in_stem(i,ind_best_hairpin(i));
    MFE2(i) = MFE(i,ind_best_hairpin(i));
end

all_hairpin_param = struct([]);
all_hairpin_param(1).n_bp = n_bp2;
all_hairpin_param(1).loop_size = loop_size2;
all_hairpin_param(1).fraction_in_stem = fraction_in_stem2;
all_hairpin_param(1).MFE = MFE2;
all_hairpin_param(1).consecutive_Us = consecutive_Us;




%% collect terminator properties (for those passing the cut!)

terminator_properties = [];
counter = 1;

for i = 1:length(positions)
    
    if bool_cut(i)
        terminator_properties(counter).position = positions(i);
        
        if forward_bool
            terminator_properties(counter).upstream_sequence = genome((positions(i)-max(lengths_folded)+1):positions(i));
        else
            terminator_properties(counter).upstream_sequence = RCsequence(genome(positions(i):(positions(i)+max(lengths_folded)-1)));
        end
        
        ind = ind_best_hairpin(i);
        terminator_properties(counter).folded_seq_id = i;
        terminator_properties(counter).RNA_structure = hairpin_param(i,ind).mountain_var;
        terminator_properties(counter).n_bp = n_bp(i,ind);
        terminator_properties(counter).loop = loop_size(i,ind);
        terminator_properties(counter).MFE = MFE(i,ind);
        terminator_properties(counter).longest_U = consecutive_Us(i);
        terminator_properties(counter).fraction_stem = fraction_in_stem(i,ind);
        terminator_properties(counter).strand = forward_bool;
        terminator_properties(counter).n_hairpins = n_hairpins(i,ind);
        
        counter = counter+1;
    end
    
end


