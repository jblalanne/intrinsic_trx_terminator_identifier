function [genome, gene_names,start,stop,strand, ...
    positions_homo_rich_f, positions_homo_rich_r,...
    full_genome_descriptor] = ...
    get_homonucleotide_stretches_v6(nucleotide_oi,...
    lower_n,upper_n,l_before_stop,l_after_stop,l_in_next_CDS,...
    dir_fasta,file_fasta,dir_gff,file_gff)


%% get genome sequences and gene annotation

[genome, full_genome_descriptor,short_genome_descriptor] = get_genomes_fna(dir_fasta,file_fasta);
fprintf('Done reading genome %s\n',file_fasta);

[gene_names,start,stop,strand] =  ...
    get_genes_gff(dir_gff,file_gff,short_genome_descriptor);
fprintf('Done parsing gene annotation %s\n',file_gff);


%% getting stretches of long homonucleotides

% key for nucleotide_oi:
% A : 1
% T : 2
% C : 3
% G : 4

% bool_region == 0: keep all
% bool_region == 1: keep CDS
% bool_region == 2: keep intergenic regions
% bool_region == 3: only downstream of genes with head genes downstream.
% bool_region == 4: downstream of genes with far co-directional genes downstream
% bool_region == 5: downstream of genes with close co-directional genes downstream

% need to break down the script in part for each genome element
% (chromosome, plasmid, etc.). 
positions_homo_rich_f = [];
positions_homo_rich_r = [];

bool_regions = [3 4 5];
for j = 1:length(bool_regions)
    for i = 1:length(genome)
        [positions_homo_rich_f{i}{j}, positions_homo_rich_r{i}{j}] = ...
            get_homonucleotide_subsequence_v3(genome{i},start{i},stop{i},strand{i},...
            bool_regions(j),nucleotide_oi,...
            lower_n,upper_n,l_before_stop,l_after_stop,l_in_next_CDS);
    end
end



% %% verification forward strand
% id_seq = 1;
% buff_u = 10;
% buff_d = 4;
% id_length_U = 3;
% id_U = 100;
% for i = 1:length(positions_homo_rich_f{id_seq}{id_length_U})
% disp(genome{id_seq}((positions_homo_rich_f{id_seq}{id_length_U}(i)-buff_u):(positions_homo_rich_f{id_seq}{id_length_U}(i)+buff_d)));
% end
% 
% 
% %% verification reverse strand
% id_seq = 2;
% buff_u = 4;
% buff_d = 10;
% id_length_U = 1;
% id_U = 100;
% for i = 1:length(positions_homo_rich_r{id_seq}{id_length_U})
% disp(genome{id_seq}((positions_homo_rich_r{id_seq}{id_length_U}(i)-buff_u):(positions_homo_rich_r{id_seq}{id_length_U}(i)+buff_d)));
% end