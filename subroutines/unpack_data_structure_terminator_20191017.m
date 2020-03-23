function terminator_properties =...
    unpack_data_structure_terminator_20191017(geometrical_cut,file_name)


load(file_name)

n_gen_el = length(all_hairpin_param_f);
n_regions = length(all_hairpin_param_f{1});


field_names = fieldnames(all_hairpin_param_f{1}{1});
field_names_oi = 1:12;

% initialize total fields
for l = 1:length(field_names_oi)
    str_oi = sprintf('%s_all = [];',...
        field_names{field_names_oi(l)});
    eval(str_oi);
end

fields_to_be_restricted = [1 2 3 4 5 7];
fields_to_be_restricted2 = [10];
fields_to_be_restricted3 = [11];
strand_str = {'f','r'};
var_str = 'all_hairpin_param';
new_var_str = 'terminator_properties';

% new fields
genome_el_id_all = [];
strand_all = [];



for j = 1:n_gen_el
    for k = 1:n_regions
        
        for m = 1:length(strand_str)
            
            % unpacking the variables
            for l = 1:length(field_names_oi)
                str_oi = sprintf('%s = %s_%s{%d}{%d}.%s;',...
                    field_names{field_names_oi(l)},var_str,strand_str{m},...
                    j,k,field_names{field_names_oi(l)});
                eval(str_oi);
            end
            
            ind_best_hairpin = select_best_hairpin_20191017(geometrical_cut,...
                MFE,n_hairpins,n_bp,loop_size,distance_stem_U,fraction_in_stem);
            

            % generating the new variable
            for l = 1:length(field_names_oi)
                
                if ismember(l,fields_to_be_restricted)
                    
                    str_oi = sprintf('%s2 =[];',field_names{field_names_oi(l)});
                    eval(str_oi);
                    
                    for n = 1:length(positions)
                        str_oi = sprintf('%s2(%d) = %s(%d,ind_best_hairpin(%d));',...
                            field_names{field_names_oi(l)},n,field_names{field_names_oi(l)},n,n);
                        eval(str_oi);
                    end
                    str_oi = sprintf('%s_all = [%s_all; %s2''];',...
                        field_names{field_names_oi(l)},field_names{field_names_oi(l)},...
                        field_names{field_names_oi(l)});
                    eval(str_oi);
                   
                elseif ismember(l,fields_to_be_restricted2)
                    
                    str_oi = sprintf('%s2 =[];',field_names{field_names_oi(l)});
                    eval(str_oi);
                    
                    for n = 1:length(positions)
                        str_oi = sprintf('%s2{%d} = %s{%d,ind_best_hairpin(%d)};',...
                            field_names{field_names_oi(l)},n,field_names{field_names_oi(l)},n,n);
                        try
                            eval(str_oi);
                        catch
                           bla 
                        end
                    end
                    str_oi = sprintf('%s_all = [%s_all; %s2''];',...
                        field_names{field_names_oi(l)},field_names{field_names_oi(l)},...
                        field_names{field_names_oi(l)});
                    eval(str_oi);
                    
                elseif ismember(l,fields_to_be_restricted3)
                    str_oi = sprintf('%s2 =[];',field_names{field_names_oi(l)});
                    eval(str_oi);
                    
                    for n = 1:length(positions)
                        str_oi = sprintf('%s2(%d,:) = squeeze(%s(%d,ind_best_hairpin(%d),:));',...
                            field_names{field_names_oi(l)},n,field_names{field_names_oi(l)},n,n);
                        eval(str_oi);
                    end
                    
                    str_oi = sprintf('%s_all = [%s_all; %s2];',...
                        field_names{field_names_oi(l)},field_names{field_names_oi(l)},...
                        field_names{field_names_oi(l)});
                    eval(str_oi);
                    
                else
                    
                    eval(sprintf('bool=size(%s,1)>1;',field_names{field_names_oi(l)}));
                    if bool
                        str_oi = sprintf('%s_all = [%s_all; %s];',...
                            field_names{field_names_oi(l)},field_names{field_names_oi(l)},...
                            field_names{field_names_oi(l)});
                    else
                        str_oi = sprintf('%s_all = [%s_all; %s''];',...
                            field_names{field_names_oi(l)},field_names{field_names_oi(l)},...
                            field_names{field_names_oi(l)});
                    end
                    try
                    eval(str_oi);
                    catch
                        bla
                    end
                    
                end
            end
            
            % need to add fields: strand, length folded index, genome element
            genome_el_id_all = [genome_el_id_all; j*ones(size(positions))];
            
            if strcmp(strand_str{m},'f')
                strand_all = [strand_all; ones(size(positions))];
            elseif strcmp(strand_str{m},'r')
                strand_all = [strand_all; zeros(size(positions))];
            end
            

        end
    end
end

% assembling the final data structure
for l = 1:length(field_names_oi)
    str_oi = sprintf('%s.%s=%s_all;',...
        new_var_str,field_names{field_names_oi(l)},...
        field_names{field_names_oi(l)});
    eval(str_oi);
end

terminator_properties.genome_el_id=genome_el_id_all;
terminator_properties.strand=strand_all;
% 
% 
% %%
% MFE_all = [];
% n_hairpin_all = [];
% n_bp_all = [];
% loop_size_all = [];
% distance_stem_U_all = [];
% fraction_in_stem_all = [];
% consecutive_U_all = [];
% strand_all = [];
% positions_all = [];
% length_folded_index_all = [];
% genome_el_id_all = [];
% upstream_gene_all = [];
% downstream_gene_all = [];
% upstream_sequence_all = [];
% mountain_var_all = [];
% 
% for s = {'f','r'}
%     for j = 1:n_gen_el
%         for k = 1:n_regions
%             
%             if strcmp(s,'f')
%                 positions = all_hairpin_param_f{j}{k}.positions;
%                 MFE = all_hairpin_param_f{j}{k}.MFE;
%                 n_hairpins =all_hairpin_param_f{j}{k}.n_hairpins;
%                 n_bp = all_hairpin_param_f{j}{k}.n_bp;
%                 loop_size = all_hairpin_param_f{j}{k}.loop_size;
%                 distance_stem_U = all_hairpin_param_f{j}{k}.distance_stem_U;
%                 fraction_in_stem = all_hairpin_param_f{j}{k}.fraction_in_stem;
%                 consecutive_Us = all_hairpin_param_f{j}{k}.consecutive_Us;
%                 upstream_gene = all_hairpin_param_f{j}{k}.gene_upstream;
%                 downstream_gene = all_hairpin_param_f{j}{k}.gene_upstream;
%                 upstream_sequence = all_hairpin_param_f{j}{k}.upstream_sequence;
%                 mountain_var = all_hairpin_param_f{j}{k}.mountain_var;
%             elseif strcmp(s,'r')
%                 positions = all_hairpin_param_r.positions{j}{k};
%                 MFE = all_hairpin_param_r.MFE{j}{k};
%                 n_hairpins =all_hairpin_param_r.n_hairpins{j}{k};
%                 n_bp = all_hairpin_param_r.n_bp{j}{k};
%                 loop_size = all_hairpin_param_r.loop_size{j}{k};
%                 distance_stem_U = all_hairpin_param_r.distance_stem_U{j}{k};
%                 fraction_in_stem = all_hairpin_param_r.fraction_in_stem{j}{k};
%                 consecutive_Us = all_hairpin_param_r.consecutive_Us{j}{k};
%             end
%             
%             % this is the length folded index
%             ind_best_hairpin = select_best_hairpin_20191017(geometrical_cut,...
%                 MFE,n_hairpins,n_bp,loop_size,distance_stem_U,fraction_in_stem,...
%                 consecutive_Us);
%             
%             % saving the properties from the best hairping from folded sequences for
%             % downstream QC purposes.
%             n_bp2 = NaN(length(positions),1);
%             loop_size2 = NaN(length(positions),1);
%             fraction_in_stem2 = NaN(length(positions),1);
%             MFE2 = NaN(length(positions),1);
%             distance_stem_U = NaN(length(positions),1);
%             n_hairpin2 = NaN(length(positions),1);
%             mountain_var2 = NaN(length(positions),1);
%             for i = 1:length(positions)
%                 n_bp2(i) = n_bp(i,ind_best_hairpin(i));
%                 loop_size2(i) = loop_size(i,ind_best_hairpin(i));
%                 fraction_in_stem2(i) = fraction_in_stem(i,ind_best_hairpin(i));
%                 MFE2(i) = MFE(i,ind_best_hairpin(i));
%                 distance_stem_U(i) = distance_stem_U(i,ind_best_hairpin(i))-consecutive_Us(i);
%                 n_hairpin2(i) = n_hairpins(i,ind_best_hairpin(i));
%             end
%             
%             length_folded_index_all = [length_folded_index_all; ind_best_hairpin];
%             consecutive_U_all = [consecutive_U_all; consecutive_Us];
%             MFE_all = [MFE_all; MFE2];
%             n_hairpin_all = [n_hairpin_all; n_hairpin2];
%             n_bp_all = [n_bp_all; n_bp2];
%             loop_size_all = [loop_size_all; loop_size2];
%             distance_stem_U_all = [distance_stem_U_all; distance_stem_U];
%             fraction_in_stem_all = [fraction_in_stem_all; fraction_in_stem2];
%             positions_all = [positions_all; positions];
%             upstream_gene_all = [upstream_gene_all; upstream_gene];
%             downstream_gene_all = [downstream_gene_all; downstream_gene];
%             upstream_sequence_all = [upstream_sequence_all; upstream_sequence];
%             mountain_var_all = [];
%             
%             
%             genome_el_id_all = [genome_el_id_all; j*ones(size(positions))];
%             
%             if strcmp(s,'f')
%                 strand_all = [strand_all; ones(size(positions))];
%             elseif strcmp(s,'r')
%                 strand_all = [strand_all; zeros(size(positions))];
%             end
%         end
%     end
% end