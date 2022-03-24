
%%% get list of species

cwd = pwd;
% the pipeline should be compatible with .fna and .gff files downloaded
% from NCBI for a species of interest. The directories below should be
% replaced with the full path to the directories containing these files.
dir_gff = cwd + "/GFFs";
dir_fasta = cwd + "/genomes_fasta";

disp(dir_gff)

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


%% run the full script on all species

now_str = datestr(datetime('now'),30);
run_dir_name = ['trx_terminator_identification_' now_str];
mkdir(run_dir_name);
cd(run_dir_name);
run_dir = pwd;

summary_plot_dir_name = ['full_summary_plots_' now_str];
mkdir(summary_plot_dir_name);
cd(summary_plot_dir_name); 
summary_plot_dir = pwd; 
cd('..')

failed_species = [];

% U-tract length searched
upper_n = 11;

GC_thresh_low = 0.35;
GC_thresh_high = 0.65;

% cut parameters
cut = [];
cut.frac_low = 0.95;
cut.MFE_cut = 9;
cut.bp_low = 5;
cut.bp_high = 15;
cut.loop_low = 3;
cut.loop_high = 8;
cut.distance_stem_end_thresh = 2;

f_pass_random = 0.02;

% number of random region to fold
N_random = 10000;

failed_species = [];

for i = 1:length(species)
     
    try
        
        tic
        stem_to_stop_master_script_random_20191007(i,species,...
            dir_fasta,files_fasta,summary_plot_dir,...
            run_dir,upper_n,cut,N_random,...
            f_pass_random)
        time_species = toc;
        fprintf('Time =  %.2f min\n',time_species/60)
    catch e
        failed_species(end+1) = i;
        msg = "Species: " + i + " failed, the error was: " + e.message + newline;
        fprintf(2,'%s', msg)
    end
    pause(2);
end

cd(cwd);

