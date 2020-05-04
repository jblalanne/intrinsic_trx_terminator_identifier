

%% list of genome sequence and annotation files

cwd = pwd;

% the pipeline should be compatible with .fna and .gff files downloaded
% from NCBI for a species of interest. The directories below should be
% replaced with the full path to the directories containing these files.
dir_gff = '/Users/GFFs';
dir_fasta = '/Users/genomes_fasta';

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


%% result files from using terminator_identification_header_random_sequence_git_20200504.m
% directory below should be replaced with directory containing the output from terminator_identification_header_random_sequence_git_20200504.m
% an example file is included: example_output_variables/0001_random_position_RNA_fold_properties_GCF_000005825.2_Bacillus_pseudofirmus_20191007T165343.mat
files_rand = dir('/Users/data_random_folding/*random_position*.mat');

%% result filres from using terminator_identification_header_script_v4_git_20200504.m
% directory below should be replaced with directory containing the output from terminator_identification_header_script_v4_git_20200504.m
% an example file is included: example_output_variables/0001_all_U_rich_upstream_RNA_fold_properties_GCF_000005825.2_Bacillus_pseudofirmus_20191017T195724.mat
files_U = dir('/Users/U_rich_folded/*all_U_rich_*.mat');


%% parsing the candidate terminator structures
% performing thresholding on hairpin parameters based on result of folding random regions in genome.

% NEED TO RUN THE ANNOTATION PARSING SCRIPT FOR SPECIES PRIOR TO RUNNING THIS PORTION: see get_genome_annotations_GC_content_git_20200504.m


% this will be the directory where final results are saved
dir_final = '/Users/final_stop_to_stem_results';
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
f_genome_cut = 0.05;        % excluding terminators from sequence element smaller than f_genome_cut fraction of largest sequence element.
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


