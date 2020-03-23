function [gene_upstream, gene_downstream] = get_upstream_gene_v2(position,strand_term,...
    genes, strand, start, stop)


% get upstream and downstream genes

% if strcmp(species,'bsub')
%     [genes, strand, start, stop] = get_bsub_genes_v2();
% elseif strcmp(species,'ecoli')
%     [genes, strand, start, stop] = get_ecoli_genes_v2();
% end

genes_f = genes(strand==1);
genes_r = genes(strand==0);

% intragenic_bool = zeros(size(position));
gene_upstream = [];
gene_downstream = [];

if strand_term
    ind_up = find(stop(strand==1)<position,1,'last');
    if (length(ind_up)==1) && (ind_up < length(genes_f))
        gene_upstream = genes_f{ind_up};
        gene_downstream = genes_f{ind_up+1};
    else
        gene_upstream = '';
        gene_downstream = '';
    end
else
    ind_up = find(start(strand==0)>position,1,'first');
    
    if length(ind_up)==1 && (ind_up > 1)
        gene_upstream = genes_r{ind_up};
        gene_downstream = genes_r{ind_up-1};
    elseif length(ind_up)==1
        gene_upstream = genes_r{ind_up};
        gene_downstream = genes_r{end};
    else
        gene_upstream = '';
        gene_downstream = '';
    end

end