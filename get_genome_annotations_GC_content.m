
%% get list of species

cwd = pwd;
% the pipeline should be compatible with .fna and .gff files downloaded
% from NCBI for a species of interest. The directories below should be
% replaced with the full path to the directories containing these files.
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



%% get genome annotations for species 
genomes = [];
full_genome_descriptors = []; 
gene_names =[];
starts = [];
stops = [];
strands = [];
 
tic
for i = 1:length(species)
    
    % get genome 
    [genomes{i}, full_genome_descriptors{i},short_genome_descriptors{i}] = get_genomes_fna(dir_fasta,files_fasta(i).name);
    fprintf('Done reading genome %s\n',files_fasta(i).name);
    
    % gene annotations
    [gene_names{i},starts{i},stops{i},strands{i}] =  ...
        get_genes_gff(dir_gff,files_gff(i).name,short_genome_descriptors{i});
    fprintf('Done parsing gene annotation %s\n',files_gff(i).name);
    
    if mod(i,ceil(length(species)/100))==0
       toc 
    end
 
end


%% get species GC content

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


