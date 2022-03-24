function [MFE, mountain_var, hairpin_param] = ...
    master_RNAfold_file_v2(genome,positions,lengths_folded,...
    constraint_bool,full_path_directory,consecutive_Us)

% current working directory
cwd = pwd;

% creating directory containing folding files
system(sprintf('mkdir %s',full_path_directory));
cd(full_path_directory)

% copying the batch_RNAfold script.
% The bash script batch_RNAfold.sh calls RNAfold to perform the RNA
% folding. RNAfold must be locally installed and on the path for the Matlab
% call of the function to work. The directory in the cp operation below
% should be replaced by the directory containing the batch_RNAfold.sh
% file.
system(sprintf('cp /Users/subroutines/batch_RNAfold.sh %s',...
    full_path_directory));

% generate subsequences to fold for each positions
write_batch_RNAfold_sequences(genome,lengths_folded,positions,...
    constraint_bool,consecutive_Us)

% any sequence to fold?
fas_files = dir('sequence*.fas');
if ~isempty(fas_files)
    
    % fold with RNAfold. Refer to RNAfold instruction if the Matlab call to
    % RNAfold is unsuccessful. 
    system('bash batch_RNAfold.sh');
    
    % read free energy data
    MFE = read_RNAfold_results(length(lengths_folded));
    
    % get mountain plot from dot-bracket representation
    [mountain_var,structure_char] = get_mountain_plot_v2(lengths_folded);
    
    % extract hairpin parameters for each folded upstream length
    hairpin_param(length(positions),length(lengths_folded)) = struct();
    hairpin_param(1,1).n_base_pairs = [];
    hairpin_param(1,1).n_hairpins = [];
    hairpin_param(1,1).loop_ind = [];
    hairpin_param(1,1).stem_3pr_ind = [];
    hairpin_param(1,1).stem_5pr_ind = [];
    hairpin_param(1,1).bulges_3pr = [];
    hairpin_param(1,1).bulges_5pr = [];
    hairpin_param(1,1).mountain_var = [];
    for i = 1:length(positions)
        for j = 1:length(lengths_folded)
            hairpin_param(i,j) = hairpin_parameter_extraction_one_subsequence(squeeze(mountain_var(i,j,:)));
        end
    end
    
else
    
    % no U-rich region to fold.
    MFE = [];
    mountain_var = [];
    hairpin_param = [];
    
end

% moving back to original directory
cd(cwd);

