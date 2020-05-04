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


% saving the properties from the best hairping from folded sequences for
% downstream QC purposes.


distance_stem_U = NaN(length(positions),1);
mountain_var = [];
gene_upstream = [];
gene_downstream = [];
upstream_sequence = [];

for i = 1:length(positions)


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
all_hairpin_param(1).mountain_var = mountain_var;
all_hairpin_param(1).positions = positions;





