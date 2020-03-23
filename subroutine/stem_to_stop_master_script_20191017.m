function stem_to_stop_master_script_20191017(ind_species,species,...
    dir_fasta,files_fasta,dir_gff,files_gff,summary_plot_dir,run_dir,GC_content_oi)



%% get U rich positions

if GC_content_oi>=0.6
    lower_n = 4;
else
    lower_n = 5;
end

if GC_content_oi<=0.35
    upper_n = 11;
else
    upper_n = 10;
end
nucleotide_oi = 2; % A:1, T:2, C:3, G: 4

% parameters restricting regions up for search for putative terminators
% relative to annotated ORFs. 
l_before_stop = 10;
l_after_stop = 200;
l_in_next_CDS = 30;

tic
[genome, gene_names,start,stop,strand, ...
    positions_U_f, positions_U_r,...
    full_genome_descriptor] = ...
    get_homonucleotide_stretches_v6(nucleotide_oi,...
    lower_n,upper_n,l_before_stop,l_after_stop,l_in_next_CDS,...
    dir_fasta,files_fasta(ind_species).name,dir_gff,files_gff(ind_species).name);


% output data structure: cell array of positions corresponding to poly-U
% sequences in the genome, stratified by genome sequence, region index,
% length of consecutive Us and position index. 
%
% positions_U_f{id_seq}{region_index}{id_length_U}(poly_U_index)
% positions_U_r{id_seq}{region_index}{id_length_U}(poly_U_index)
%
% id_seq: the index of the different sequence elements in the original
%         fasta file (e.g., chromosome, plasmids). Same index as genome.
% region_index: 1 corresponds to poly U downstream of head on genes, 2 to
%               downstrea of co-directional genes, with long distance to
%               the next gene, 3 to downstream of nearby co-directional
%               genes.
% id_length_U: index for the poly U tract length, spans lower_n to upper_n.
% poly_U_index: index of the positions identified
fprintf('Done identifying U rich regions %s\n',species{ind_species});

 
%% fold upstream region: get hairpin characteristics

% time at start
now_str = datestr(datetime('now'),30);

% lengths of folded upstream RNA regions for each hairpin.
lengths_folded = [30 35 40 45 50];


% getting the short species name
inds = regexp(full_genome_descriptor{1},' ');
short_species_name = strrep(full_genome_descriptor{1}((inds(1)+1):(inds(3)-1)),' ','_');
short_species_name = strrep(short_species_name,'''','');

% making the new directories.
cd(run_dir)
cwd = pwd;
folding_directory = [cwd '/' sprintf('%d_',ind_species) species{ind_species} '_' strrep(short_species_name,'''','') '_folding_upstream_polyUs_' now_str];
putative_terminator_properties = identify_upstream_hairpins_v4(...
    genome, gene_names, start, stop, strand, ...
    positions_U_f, positions_U_r,...
    folding_directory,lower_n,upper_n,lengths_folded,...
    full_genome_descriptor,species{ind_species},now_str,summary_plot_dir,...
    short_species_name,ind_species);
fprintf('Done identifying putative terminators s %s\n',species{ind_species});




