function terminator_properties = identify_upstream_hairpins_v4(...
    genome, gene_names, start, stop, strand, ...
    positions_U_f,positions_U_r,...
    folding_directory,lower_n,upper_n,lengths_folded,...
    full_genome_descriptor,species,now_str,summary_plot_dir,...
    short_species_name,ind_species)


% housekeeping directory creation
mkdir(folding_directory);
cd(folding_directory);

% loop through the different sequences in the original fasta file
% and the different genomic region under consideration (stratification in 
% the position of putative terminators). 

terminator_properties = [];
terminator_properties_f = [];
terminator_properties_r = [];
stem_to_stop = [];
all_hairpin_param_f = [];
all_hairpin_param_r = [];

positions_f = [];
MFE_f = [];
n_hairpins_f = [];
n_bp_f = [];
loop_size_f = [];
distance_stem_3pr_f = [];
fraction_in_stem_f = [];
consecutive_Us_f = [];
positions_r = [];
MFE_r = [];
n_hairpins_r = [];
n_bp_r = [];
loop_size_r = [];
distance_stem_3pr_r = [];
fraction_in_stem_r = [];
consecutive_Us_r = [];


final_stem_to_stop_f = [];
final_stem_to_stop_r = [];


for i = 1:length(genome)
    for j = 1:length(positions_U_f{i})
        
        % forward strand
        forward_bool = 1;        
            all_hairpin_param_f{i}{j} = get_upstream_hairpins_strand_v6(...
            genome{i}, gene_names{i}, start{i}, stop{i}, strand{i},...
            positions_U_f{i}{j}, positions_U_r{i}{j}, forward_bool,...
            [folding_directory sprintf('/seq_%d_region_%d',i,j)],...
            lower_n, upper_n,lengths_folded);
        
        % reverse strand
        forward_bool = 0;
            all_hairpin_param_r{i}{j} = get_upstream_hairpins_strand_v6(...
            genome{i}, gene_names{i}, start{i}, stop{i}, strand{i},...
            positions_U_f{i}{j}, positions_U_r{i}{j}, forward_bool,...
            [folding_directory sprintf('/seq_%d_region_%d',i,j)],...
            lower_n, upper_n,lengths_folded);
        
    end
end


% saving final variables
save(sprintf('all_U_rich_upstream_RNA_fold_properties_%s_%s_%s.mat',species,short_species_name,now_str),...
    'all_hairpin_param_f','all_hairpin_param_r');


cd(folding_directory)
cd ..

