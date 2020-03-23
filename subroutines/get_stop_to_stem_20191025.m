function distance_stem_stop = get_stop_to_stem_20191025(U_rich_hairpins,...
    start,stop,strand,genome,bool_cuts,f_genome_cut)


% stop to stem distance for terminators passing the local criteria
% (bool_cuts). 

% for i = 1:length(bool_cuts)
[distance_stem_stop_all,closest_upstream_stop_all] = ...
    closest_upstream_stop_codon_20191022(start, stop, strand,U_rich_hairpins);
% end

distance_stem_stop = [];
closest_upstream_stop = [];

for i = 1:length(bool_cuts)
    distance_stem_stop{i} = NaN(size(bool_cuts{i}));
    closest_upstream_stop{i} = NaN(size(bool_cuts{i}));

    for j = 1:length(bool_cuts{i})
        if bool_cuts{i}(j)
            distance_stem_stop{i}(j) = distance_stem_stop_all(j);
            closest_upstream_stop{i}(j) = closest_upstream_stop_all(j);
        end
    end
end

% multiple terminators downstream of same gene
bool_multi_terminators = [];
for i = 1:length(bool_cuts)
    bool_multi_terminators{i} = exclude_multiple_terminators_v3(U_rich_hairpins,closest_upstream_stop{i});
end

% terminators on "minor" (using the size of the element as a proxy) genomic element (e.g., plasmid) excluded. 
max_genome_size = 0;
for i = 1:length(genome)
    max_genome_size = max([max_genome_size length(genome{i})]);
end
gen_id_ok = [];
for i = 1:length(genome)
    if length(genome{i})>=f_genome_cut*max_genome_size
        gen_id_ok(i) = 1;
    else
        gen_id_ok(i) = 0;
    end
end
bool_gen_id = [];
for i = 1:length(bool_cuts)
    for j = 1:length(bool_cuts{i})
        bool_gen_id{i}(j) = gen_id_ok(U_rich_hairpins.genome_el_id(j));
    end
end


% final boolean
final_bool = [];
for i = 1:length(bool_cuts)
    
    final_bool{i} = bool_cuts{i} & bool_multi_terminators{i} & bool_gen_id{i}';
    distance_stem_stop{i}(~final_bool{i})=NaN;
end






