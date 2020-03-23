
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


