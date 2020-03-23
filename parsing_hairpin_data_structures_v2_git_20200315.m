

%% list of genome sequence and annotation files

cwd = pwd;
dir_gff = '/Users/jbl/Documents/MIT/research/D1/trx_trl_coupling/basic_bioinformatics/RefSeq_representative_reference_genomes/GFFs';
dir_fasta = '/Users/jbl/Documents/MIT/research/D1/trx_trl_coupling/basic_bioinformatics/RefSeq_representative_reference_genomes/genomes_fasta';

cd(dir_fasta);
files_fasta = dir('*.fna');
species = [];
display_name = [];
for i = 1:length(files_fasta)
    name = files_fasta(i).name;
    ind = find(name=='_');
    species{i} = name(1:(ind(2)-1));
    display_name{i} = species{i};
end
cd(cwd);

cd(dir_gff)
files_gff = dir('*.gff');
cd(cwd);



%%

files_rand = dir('/Users/jbl/Documents/MIT/research/D1/trx_trl_coupling/basic_bioinformatics/random_position_folding_Refseq_20191011/data_random_folding/*random_position*.mat');
files_U = dir('/Users/jbl/Documents/MIT/research/D1/trx_trl_coupling/basic_bioinformatics/U_rich_folded_Refseq_v4_20191021/*all_U_rich_*.mat');

% GC content
load('Refseq_GC_content_20191002.mat');
% genome and gene annotations

%%
% genomes and gene annotations
% load('genomes_gene_annotations_Refseq_20191017.mat');



%%

dir_final = '/Users/jbl/Documents/MIT/research/D1/trx_trl_coupling/basic_bioinformatics/final_stop_to_stem_pipeline_20191026';
cd(dir_final);


% loading appropriate files 

for ind_oi = 1:length(species)
    
now_str = datestr(datetime('now'),30);

tic

file_name_U = [files_U(ind_oi).folder '/' files_U(ind_oi).name];
file_name_rand = [files_rand(ind_oi).folder '/' files_rand(ind_oi).name];

ind1 = regexp(files_U(ind_oi).name,'GCF_');
ind2 = regexp(files_U(ind_oi).name,'_2019')-1;
species_name = files_U(ind_oi).name(ind1:ind2);
ind1 = regexp(species_name,'_');
GCF_name = species_name(1:(ind1(2)-1));
short_name = species_name((ind1(2)+1):end);

GC_content_oi = GC_content(ind_oi);

genome = genomes{ind_oi};
start = starts{ind_oi};
stop = stops{ind_oi};
strand = strands{ind_oi};



% % % % % % % % % % % % % % % 
% unpacking data structures %
% % % % % % % % % % % % % % % 

%loose cut parameters for the first pass (just to select the best hairpin).
geometrical_cut = [];
geometrical_cut.bp_low = 5;
geometrical_cut.bp_high = 15;
geometrical_cut.loop_low = 3;
geometrical_cut.loop_high = 8;
geometrical_cut.distance_stem_end_thresh = 0;   % MODIFICATION on 03/15/2020

% best_hairpins_struct
U_rich_hairpins=...
    unpack_data_structure_terminator_20191017(geometrical_cut,file_name_U);

% read random hairpins data
random_hairpins = load_random_hairpin_data_20191019(file_name_rand);


% % % % % % % % % % % % % % % % % % % %
% thresholding on hairpin parameters %
% % % % % % % % % % % % % % % % % % % %
MFE_thresholds = 0:0.1:35;

cuts = [];
% variable portion of the cut
cuts.f_pass = [1.5E-2 1E-2];
cuts.frac_low = [0.9 0.95];
cuts.MFE = [6.5 6.5];

% fixed portion of the cut
cuts.bp_low = 5*ones(size(cuts.f_pass));
cuts.bp_high = 15*ones(size(cuts.f_pass));
cuts.loop_low = 3*ones(size(cuts.f_pass));
cuts.loop_high = 8*ones(size(cuts.f_pass));
cuts.distance_stem_end_thresh = 0*ones(size(cuts.f_pass)); % MODIFICATION on 03/15/2020

% thresholding on hairpin parameters: use random hairpin properties to select MFE threhsold 
plot_bool = 0;
[bool_cuts,dG_cuts,cuts,fraction_in_U,fraction_in_rand] = get_final_cut_20191019(MFE_thresholds,...
   U_rich_hairpins,random_hairpins,cuts,plot_bool);


% % % % % % % % % % % % % % % % % % 
% obtaining stop-to-stem distance %
% % % % % % % % % % % % % % % % % % 
f_genome_cut = 0.05;
tic
stop_to_stem = get_stop_to_stem_20191025(...
    U_rich_hairpins,start,stop,strand,genome,bool_cuts,f_genome_cut);


% % % % % % % % % % % % % 
% plotting the results %
% % % % % % % % % % % % % 

h_fig = plotting_Urich_random_hairpins_20191024(GCF_name,short_name,...
    GC_content_oi,...
    U_rich_hairpins,random_hairpins,cuts,...
    bool_cuts,stop_to_stem,...
    MFE_thresholds,fraction_in_rand,fraction_in_U,dG_cuts);


% saving figure
saveas(h_fig,[files_U(ind_oi).name(1:4) '_' GCF_name '_' short_name '_stop_to_stem_plot_' now_str '.png']);
saveas(h_fig,[files_U(ind_oi).name(1:4) '_' GCF_name '_' short_name '_stop_to_stem_plot_' now_str '.fig']);
close(h_fig);

% % saving data files
save([files_U(ind_oi).name(1:4) '_' GCF_name '_' short_name '_stop_to_stem_variables_' now_str '.mat'],...
    'stop_to_stem','bool_cuts','dG_cuts','fraction_in_U',...
    'fraction_in_rand','U_rich_hairpins','random_hairpins');
end


