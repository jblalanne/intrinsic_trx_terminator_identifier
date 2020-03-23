function [gene_names,start,stop,strand] =  ...
    get_genes_gff(gff_dir,gff_file_name,short_genome_descriptor)

gene_names = [];
start = [];
stop = [];
strand = [];

cwd = pwd;
cd(gff_dir);

fid = fopen(gff_file_name);
gff_str_fmt = '%s %s %s %d %d %s %s %s %s';

    
L_gene_prealloc = 2E4;  % largest # of bacterial genes is ~10^4;
tline = fgetl(fid);
ind_seq = [];

tic 
while ~feof(fid) 
    
    
    % searching for type of containing genetic material (different sequence
    % regions appearing under different id in the fasta file, e.g.,
    % chromosomal DNA vs. plasmids). 
    ind = regexp(tline,'##sequence-region','ONCE');
    if ~isempty(ind)
        % matching the sequence descriptor
        ind_spaces = find(tline==' ');
        seq_descriptor = tline((ind_spaces(1)+1):(ind_spaces(2)-1));
        ind_seq = find(strcmp(seq_descriptor,short_genome_descriptor));
        
        
        % initialization & pre-allocation
        start{ind_seq} = NaN(L_gene_prealloc,1);
        stop{ind_seq} = NaN(L_gene_prealloc,1);
        strand{ind_seq} = NaN(L_gene_prealloc,1);
        gene_names{ind_seq} = cell(L_gene_prealloc,1);
        gene_counter = 1;
       
        
%         % skip the next line (##species ... etc.). 
%         tline=fgetl(fid);
%         tline=fgetl(fid);
    end    
    
    
    % reading the gene information line by line
    
    if ~isempty(ind_seq) && ~strcmp(tline(1),'#') %&& ~isempty(tline)
        % read line
        
        line_content = textscan(tline,gff_str_fmt,'Delimiter','\t');
        
        
        % identifying whether the line corresponds to a CDS
        ind_cds_bool = regexp(line_content{9}{1},'gene_biotype=protein_coding','ONCE');
        
        % need the additional protein coding flag, as non coding RNA (e.g.,
        % rRNA, tRNA, SRP_RNA) are also labelled as genes.
        if strcmp(line_content{3},'gene') && ~isempty(ind_cds_bool)
            
            % gene positions and strand
            start{ind_seq}(gene_counter) = line_content{4};
            stop{ind_seq}(gene_counter) = line_content{5};
            strand{ind_seq}(gene_counter) = (line_content{7}{1}=='+');
            
            % gene name
            ind_gene_name = regexp(line_content{9}{1},'Name=');
            ind_semicol = find(line_content{9}{1}==';' & (1:length(line_content{9}{1}))>ind_gene_name,1,'first');
            gene_names{ind_seq}{gene_counter} = line_content{9}{1}( (ind_gene_name+length('Name=')):(ind_semicol-1));
            
            gene_counter = gene_counter+1;
            
        end
    end
    
    tline=fgetl(fid);

end


% removing unallocated space
for i = 1:length(start)
    ind_unallocated = isnan(start{i});
    start{i}(ind_unallocated) = [];
    stop{i}(ind_unallocated) = [];
    strand{i}(ind_unallocated) = [];
    gene_names{i}(ind_unallocated) = [];
end

fclose(fid);    % housekeeping
cd(cwd);