function bool_multi_terminators = exclude_multiple_terminators_v3(terminator_properties,upstream_stop)

positions = [terminator_properties(:).positions];
strand = [terminator_properties(:).strand];
gen_el = [terminator_properties(:).genome_el_id];

connectivity_matrix = zeros(length(terminator_properties));
for i = 1:length(positions)
    for j = 1:length(positions)
        connectivity_matrix(i,j) = (strand(i)==strand(j)) && (gen_el(i)==gen_el(j)) && (upstream_stop(i)==upstream_stop(j)) && ~isnan(upstream_stop(i));
    end
end

% restrict to isolated terminators
bool_multi_terminators = zeros(length(terminator_properties.n_bp),1);
[~,sizes_components,members] = networkComponents(connectivity_matrix);
clusters = members(sizes_components==1);
for j = 1:length(clusters)
    bool_multi_terminators(clusters{j}(1)) = 1;
end


