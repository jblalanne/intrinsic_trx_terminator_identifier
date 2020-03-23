function stem_to_stop_master_script_random_20191007(ind_species,species,...
    dir_fasta,files_fasta,summary_plot_dir,run_dir,...
    upper_n,cut,N_random,...
    f_pass_random)


%% get U rich positions

% nucleotide_oi = 2; % A:1, T:2, C:3, G: 4
% 

[genome, full_genome_descriptor,short_genome_descriptor] = get_genomes_fna(dir_fasta,files_fasta(ind_species).name);
fprintf('Done reading genome %s\n',files_fasta(ind_species).name);

%% folding random regions

% assumes that the first sequence element is the largest
L_genome = [];
for i =1:length(genome)
    L_genome(i) = length(genome{i});
end

rand_pos_f = [];
rand_pos_r = [];
for i = 1:length(genome)
    rand_pos_f{i}{1}{1} = randi(length(genome{i}),ceil(L_genome(i)/sum(L_genome)*N_random/2),1);
    rand_pos_r{i}{1}{1} = randi(length(genome{i}),ceil(L_genome(i)/sum(L_genome)*N_random/2),1);
end



%% fold upstream RANDOM region: get hairpin characteristics


% time at start
now_str = datestr(datetime('now'),30);

% lengths of folded upstream RNA regions for each hairpin.
lengths_folded = [40];

% getting the short species name
inds = regexp(full_genome_descriptor{1},' ');
short_species_name = strrep(full_genome_descriptor{1}((inds(1)+1):(inds(3)-1)),' ','_');
short_species_name = strrep(short_species_name,'''','');


% making the new directories.
cd(run_dir)
cwd = pwd;
folding_directory = [cwd '/' sprintf('%d_',ind_species) species{ind_species} '_' short_species_name '_folding_upstream_RANDOM_' now_str];

[putative_terminator_properties,dG_cut] = get_random_hairpins_20191007(...
    genome, ...
    rand_pos_f, rand_pos_r,...
    folding_directory,lengths_folded,cut,...
    full_genome_descriptor,species{ind_species},now_str,summary_plot_dir,...
    short_species_name,ind_species,f_pass_random);
fprintf('Done identifying putative terminators s %s\n',species{ind_species});


