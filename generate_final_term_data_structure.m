
cwd = pwd;

% directory below should contains the output of script: parsing_hairpin_data_structures_v2_git_20200504.m
dir_oi = cwd + "/final_stop_to_stem_results/";
cd(dir_oi);

file_oi = dir(cwd + "/final_stop_to_stem_results/*stop_to_stem_variables*.mat");
% an example file is included: example_output_variables/0001_GCF_000005825.2_Bacillus_pseudofirmus_stop_to_stem_variables_20191026T163538.mat

stop_to_stem_AND = [];
stop_to_stem_OR = [];
f_pass_AND = [];
f_pass_OR = [];

%% run annotation reading and GC content script:
% see get_genome_annotations_GC_content_git_20200504.m

max_l = 0;
for i = 1:length(species)
    for j = 1:length(full_genome_descriptors{i})
       max_l = max([max_l length( full_genome_descriptors{i}{j})]);
    end
end

%% generate final "dereplicated" stop to stem distribution. 

putative_terminators = [];
for i = 1:length(species)
    putative_terminators(i).species = [];
    putative_terminators(i).genomic_el_id = [];
    putative_terminators(i).genomic_el_descriptor = [];
    putative_terminators(i).position = [];
    putative_terminators(i).strand = [];
    putative_terminators(i).gene_upstream = [];
    putative_terminators(i).upstream_sequence = [];
    putative_terminators(i).dot_bracket = [];
    putative_terminators(i).n_bp = [];
    putative_terminators(i).loop_size = [];
    putative_terminators(i).n_consecutive_Us = [];
    putative_terminators(i).MFE = [];
    putative_terminators(i).fraction_in_stem = [];
    putative_terminators(i).stop_to_stem = [];
%     putative_terminators(i).distance_stem_U = [];
end


stop_to_stem_OR_all = [];
L_mountain = 50;

MFE_thresh_4Us = 20;

tic
for i = 1:length(file_oi)
    load(file_oi(i).name)
    
    bool_OR = ~isnan(stop_to_stem{1}) | ~isnan(stop_to_stem{2});
    stop_to_stem_OR = NaN(length(bool_OR),1);
    for j = 1:length(bool_OR)
        if bool_OR(j)
            if ~isnan(stop_to_stem{1}(j))
                stop_to_stem_OR(j) = stop_to_stem{1}(j);
            elseif ~isnan(stop_to_stem{2}(j))
                stop_to_stem_OR(j) = stop_to_stem{2}(j);
            end
        end
    end
    
    up_seqs = {U_rich_hairpins.upstream_sequence};
    up_seqs = up_seqs{1};
    n_L = length(up_seqs);
    n_identical_seqs = NaN(n_L,1);

    temp_OR = [];
    
    bool_OR2 = bool_OR;
    bool_OR2_note = bool_OR;
    for j = 1:n_L
        
%         additional threshold based on C. crescentus and 4 U's:
        if U_rich_hairpins.consecutive_Us(j)==4 && bool_OR2(j)
            bool_OR2(j) = U_rich_hairpins.MFE(j)<-MFE_thresh_4Us;
            bool_OR2_note(j) = U_rich_hairpins.MFE(j)<-MFE_thresh_4Us;
        end

        % constrain the distance between stem and U tract (added: 03/14/2020)
        if U_rich_hairpins.distance_stem_U(j)>0
            bool_OR2(j) = 0;
        end

            
        if bool_OR2(j)
            
            inds = find(strcmp(up_seqs,up_seqs{j}));
            inds2 = find(stop_to_stem_OR(inds)==stop_to_stem_OR(j));
            n_identical_seqs(j) = length(inds2);
            temp_OR = [temp_OR stop_to_stem_OR(j)];
            
            % generating terminator structure
            putative_terminators(i).species{end+1} = species{i};
            putative_terminators(i).genomic_el_id{end+1} = U_rich_hairpins.genome_el_id(j);
            putative_terminators(i).genomic_el_descriptor{end+1} = full_genome_descriptors{i}{U_rich_hairpins.genome_el_id(j)};
            putative_terminators(i).position(end+1) = U_rich_hairpins.positions(j);
            putative_terminators(i).strand(end+1) = U_rich_hairpins.strand(j);
            putative_terminators(i).gene_upstream{end+1} = U_rich_hairpins.gene_upstream{j};
            putative_terminators(i).upstream_sequence{end+1} = U_rich_hairpins.upstream_sequence{j};
            range_mountain = (L_mountain-length(U_rich_hairpins.upstream_sequence{j})+1):L_mountain;
            dot_bracket_oi = mountain_var_to_dot_bracket(squeeze(U_rich_hairpins.mountain_var(j,range_mountain)));
            for k = 1:U_rich_hairpins.consecutive_Us(j)
                dot_bracket_oi(end-k+1) = 'x';
            end
            putative_terminators(i).dot_bracket{end+1} = dot_bracket_oi;
            putative_terminators(i).n_bp(end+1) = U_rich_hairpins.n_bp(j);
            putative_terminators(i).loop_size(end+1) = U_rich_hairpins.loop_size(j);
            putative_terminators(i).n_consecutive_Us(end+1) = U_rich_hairpins.consecutive_Us(j);
            putative_terminators(i).MFE(end+1) = U_rich_hairpins.MFE(j);
            putative_terminators(i).fraction_in_stem(end+1) = U_rich_hairpins.fraction_in_stem(j);
            putative_terminators(i).stop_to_stem(end+1) = stop_to_stem_OR(j);
           
            % remove other identical sequences to avoid double counting
            for k = 1:length(inds2)
                bool_OR2(inds(inds2(k))) = 0;
            end
        end
    end
    
    if sum(bool_OR2)==0
        putative_terminators(i).species{1} = species{i};
    end
    
    % total final number of terminator per species 
    n_term(i) = length(temp_OR);
    
    stop_to_stem_OR_all{i} = temp_OR;
    
    f_unique(i) = sum(n_identical_seqs==1)/sum(~isnan(n_identical_seqs));
    
    if mod(i,100)==0
       toc 
    end
end



%% printing final table

file_name = 'putative_terminators.txt';
fid = fopen(file_name,'w');

header = [];
header = [header 'Species\t'];
header = [header 'Genomic element #\t'];
header = [header 'Genomic element descriptor\t'];
header = [header 'Position\t'];
header = [header 'Strand\t'];
header = [header 'Upstream gene\t'];
header = [header 'Upstream sequence\t'];
header = [header 'RNA structure\t'];
header = [header 'Stem size\t'];
header = [header 'Loop size\t'];
header = [header 'MFE (kcal/mol)'];
header = [header 'Fraction paired in stem\t'];
header = [header 'Stop to stem distance\n'];
fprintf(fid,header);


l_descriptor = 0;

str_fmt = '%20s\t%3d\t%80s\t%10d\t%3d\t%20s\t%50s\t%50s\t%3d\t%3d\t%.1f\t%.2f\t%4d\n';

n_total_term = 0;
tic
for i = 1:length(putative_terminators)
    for j = 1:length(putative_terminators(i).position)
        fprintf(fid,str_fmt,...
            putative_terminators(i).species{j},...
            putative_terminators(i).genomic_el_id{j},...
            putative_terminators(i).genomic_el_descriptor{j},...
            putative_terminators(i).position(j),...
            putative_terminators(i).strand(j),...
            putative_terminators(i).gene_upstream{j},...
            putative_terminators(i).upstream_sequence{j},...
            putative_terminators(i).dot_bracket{j},...
            putative_terminators(i).n_bp(j),...
            putative_terminators(i).loop_size(j),...
            putative_terminators(i).MFE(j),...
            putative_terminators(i).fraction_in_stem(j),...
            putative_terminators(i).stop_to_stem(j));
        n_total_term = n_total_term+1;
    end
    
    if mod(i,100)==0
        toc
    end
end

fclose(fid);

