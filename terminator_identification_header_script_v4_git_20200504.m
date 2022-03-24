
%% get list of species' GFF (annotation) and fasta (sequence) files

cwd = pwd;
% the pipeline should be compatible with .fna and .gff files downloaded
% from NCBI for a species of interest. The directories below should be
% replaced with the full path to the directories containing these files.
% Example files for B. pseudofirmus are included. 
dir_gff = cwd + "/GFFs";
dir_fasta = cwd + "/genomes_fasta";

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


%% GC content variable for species: 

sizes_genomes = zeros(2,length(species)*10);
tic
counter = 1;
for i = 1:length(species)
    [genome, full_genome_descriptor,short_genome_descriptor] = get_genomes_fna(dir_fasta,files_fasta(i).name);
    combined_genome = [];
    for j = 1:length(genome)
        combined_genome = [combined_genome genome{j}];
    end
    GC_content(i) = sum(combined_genome=='C' | combined_genome=='G')/length(combined_genome);
end

 

%% run the full script on all species


% creating directories for species. 
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

% looping through all species (list from the GFF and FASTA file list above)
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
        failed_species(end+1) = i;  % sometimes non-standard characters in genomes can lead script to fail.
    end
end

cd(cwd);

