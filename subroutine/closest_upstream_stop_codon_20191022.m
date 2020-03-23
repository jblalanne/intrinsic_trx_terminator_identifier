function [distance_stem_stop, closest_upstream_stop] = ...
    closest_upstream_stop_codon_20191022(start, stop, strand_gene,...
    terminator_summary)


positions = [terminator_summary(:).positions];
strand_term = [terminator_summary(:).strand];
n_hairpins = [terminator_summary(:).n_hairpins];


% initialization
closest_upstream_stop = NaN(length(positions),1);
distance_upstream_stop = NaN(length(positions),1);
distance_stem_stop = NaN(length(positions),1);


% looping through positions of the end of the U-tract.
for i = 1:length(positions)
    
    if n_hairpins(i)>0
        
        mountain_plot = terminator_summary.mountain_var(i,:);
        gen_el_id = terminator_summary.genome_el_id(i);
        
        % startifying gene annotation based on strand & genome sequence
        % element. 
        stop_forward = stop{gen_el_id}(strand_gene{gen_el_id}==1);
        start_reverse = start{gen_el_id}(strand_gene{gen_el_id}==0);
        
        % fine position details
        stop_forward = stop_forward-1;
        start_reverse = start_reverse+1;

        
        
        if strand_term(i)
            
            bool_upstream = stop_forward<positions(i);
            position_upstream = stop_forward(bool_upstream);
            
            if ~isempty(position_upstream)
                closest_upstream_stop(i) = position_upstream(end);
                distance_upstream_stop(i) = positions(i)-closest_upstream_stop(i);
                
                % distance between 5' stem and stop codon:
                stem_5pr = find(mountain_plot>0,1,'first');
                stem_5pr_end_position = positions(i)-(length(mountain_plot)-stem_5pr);
                distance_stem_stop(i) =  stem_5pr_end_position-closest_upstream_stop(i);
                
            end
            
        else
            
            bool_upstream = start_reverse>positions(i);
            position_upstream = start_reverse(bool_upstream);
            
            if ~isempty(position_upstream)
                closest_upstream_stop(i) = position_upstream(1);
                distance_upstream_stop(i) = -positions(i)+closest_upstream_stop(i);
                
                % distance between 5' stem and stop codon:
                stem_5pr = find(mountain_plot>0,1,'first');
                stem_5pr_end_position = positions(i)+(length(mountain_plot)-stem_5pr);
                distance_stem_stop(i) =  -stem_5pr_end_position+closest_upstream_stop(i);
            end
        end

    end
    
end
