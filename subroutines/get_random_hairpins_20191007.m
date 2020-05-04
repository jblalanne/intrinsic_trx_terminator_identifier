function [terminator_properties, dG_cut]= get_random_hairpins_20191007(...
    genome, ...
    rand_pos_f,rand_pos_r,...
    folding_directory,lengths_folded,cut,...
    full_genome_descriptor,species,now_str,summary_plot_dir,...
    short_species_name,ind_species,f_pass_random)


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
    for j = 1:length(rand_pos_f{i})
        
        % forward strand
        forward_bool = 1; 
        
        [terminator_properties_f{i}{j}, all_hairpin_param_f{i}{j},...
            positions_f{i}{j},MFE_f{i}{j},n_hairpins_f{i}{j},n_bp_f{i}{j},...
            loop_size_f{i}{j},fraction_in_stem_f{i}{j}] = ...
            get_random_hairpins_strand_20191007(...
            genome{i}, ...
            rand_pos_f{i}{j}, rand_pos_r{i}{j}, forward_bool,...
            [folding_directory sprintf('/seq_%d_region_%d',i,j)],lengths_folded,cut);
        
        
        % reverse strand
        forward_bool = 0;
         [terminator_properties_r{i}{j}, all_hairpin_param_r{i}{j},...
            positions_r{i}{j},MFE_r{i}{j},n_hairpins_r{i}{j},n_bp_r{i}{j},...
            loop_size_r{i}{j},fraction_in_stem_r{i}{j}] = ...
            get_random_hairpins_strand_20191007(...
            genome{i}, ...
            rand_pos_r{i}{j}, rand_pos_r{i}{j}, forward_bool,...
            [folding_directory sprintf('/seq_%d_region_%d',i,j)],lengths_folded,cut);
        
        % combining important variables
        terminator_properties{i}{j} = [terminator_properties_f{i}{j} terminator_properties_r{i}{j}];
        
    end
end


% plot full summary combining all regions from genomic element >f_cut*(size
% of maximum element).
f_cut = 0;
[dG_cut,MFE,n_bp,loop_size,fraction_in_stem,n_hairpin] = generate_combined_summary_plot_random(species,short_species_name,genome,rand_pos_f,cut,...
    summary_plot_dir,f_cut,all_hairpin_param_f,all_hairpin_param_r,...
    terminator_properties_f,terminator_properties_r,now_str,ind_species,...
    n_hairpins_f,n_hairpins_r,f_pass_random);


save(sprintf('random_hairpins_terminators_%s_%s_%s.mat',species,short_species_name,now_str),'terminator_properties','stem_to_stop')
save(sprintf('random_position_RNA_fold_properties_s_%s_%s_%s.mat',species,short_species_name,now_str),...
    'MFE','n_bp','loop_size','fraction_in_stem','n_hairpin');



cd(folding_directory)
cd ..

