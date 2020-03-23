
%% get list of species

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


load('Refseq_GC_content_20191002.mat');


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

for i = 1:length(species)
    
    GC_content_oi = GC_content(i);
    try
        tic
        stem_to_stop_master_script_20191017(i,species,...
            dir_fasta,files_fasta,dir_gff,files_gff,summary_plot_dir,...
            run_dir,GC_content_oi)
        time_species = toc;
        fprintf('Time =  %.2f min\n',time_species/60)
    catch
        failed_species(end+1) = i;
    end
end

cd(cwd);

