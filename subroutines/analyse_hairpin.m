function [n_bp, loop_size, distance_stem_3pr, fraction_in_stem, n_hairpins] = analyse_hairpin(RNA_structure)


% number of hairpins (should be 1 for this analysis to make sense).
n_hairpins = length(findpeaks(RNA_structure));

% number of base pairs
n_bp = max(RNA_structure);

% size of loop
loop_size = sum(RNA_structure==max(RNA_structure))-1;

% distance stem from 3' end 
last_paired_base = find(RNA_structure>0,1,'last')+1;
if ~isempty(last_paired_base)
    distance_stem_3pr=length(RNA_structure)-last_paired_base;
else
    distance_stem_3pr=length(RNA_structure);
end


% fraction of bases in hairpin (assuming only one hairpin)
hairpin_start = find(RNA_structure>0,1,'first')-1;
hairpin_end = find(RNA_structure>0,1,'last')+1;
if max(RNA_structure)==0
    fraction_in_stem = 0;
else
    fraction_in_stem = (2*n_bp)/(hairpin_end-hairpin_start-loop_size);
end
