function [dG_cut,MFE,n_bp,loop_size,fraction_in_stem,n_hairpin] =...
    generate_combined_summary_plot_random(species,short_species_name,genome,...
    rand_pos_f,cut,...
    summary_plot_dir,f_cut,all_hairpin_param_f,all_hairpin_param_r,...
    terminator_properties_f,terminator_properties_r,now_str,ind_species,...
    n_hairpins_f,n_hairpins_r,f_pass_random)


% moving directories
cwd = pwd;
cd(summary_plot_dir)


max_size_genome_element = 0;
for j = 1:length(genome)
    max_size_genome_element = max([max_size_genome_element length(genome{j})]);
end

%%
field_names = fieldnames(all_hairpin_param_f{1}{1});
strand_str = {'f','r'};
var_str = {'all_hairpin_param_'};

for m = 1:length(var_str)
    for k = 1:length(strand_str)
        eval(sprintf('%splot_%s{1}{1} = [];',var_str{m},strand_str{k}));
        for l = 1:length(field_names)
            
            eval(sprintf('%splot_%s{1}{1}.%s = [];',var_str{m},strand_str{k},field_names{l}));
            
            for i = 1:length(genome)
                if length(genome{i})>max_size_genome_element*f_cut
                    for j = 1:length(rand_pos_f{i})
                        
                        
                        eval(sprintf('%splot_%s{1}{1}.%s = [%splot_%s{1}{1}.%s %s%s{%d}{%d}.%s''];',...
                            var_str{m},strand_str{k},field_names{l},var_str{m},strand_str{k},...
                            field_names{l},var_str{m},strand_str{k},i,j,field_names{l}));
                    end
                end
            end
        end
    end
end

%% combining other variables

terminator_properties_plot_f{1}{1} = [];
terminator_properties_plot_r{1}{1} = [];

combined_genome = [];
n_hairpin_plot_f = [];
n_hairpin_plot_r = [];

for i = 1:length(genome)
    if length(genome{i})>max_size_genome_element*f_cut
        
        combined_genome = [combined_genome genome{i}];
        
        for j = 1:length(rand_pos_f{i})
            terminator_properties_plot_f{1}{1} = [terminator_properties_plot_f{1}{1} terminator_properties_f{i}{j}];
            terminator_properties_plot_r{1}{1} = [terminator_properties_plot_r{1}{1} terminator_properties_r{i}{j}];
            n_hairpin_plot_f = [n_hairpin_plot_f n_hairpins_f{i}{j}'];
            n_hairpin_plot_r = [n_hairpin_plot_r n_hairpins_r{i}{j}'];
        end
    end
end



% generating the summary plot
GC_content = sum(combined_genome=='C' | combined_genome=='G')/length(combined_genome);
[dG_cut,MFE,n_bp,loop_size,fraction_in_stem,n_hairpin] =...
    plot_hairpin_properties_summary_random(all_hairpin_param_plot_f{1}{1}, all_hairpin_param_plot_r{1}{1},...
    terminator_properties_plot_f{1}{1},terminator_properties_plot_r{1}{1}, cut,...
    species,short_species_name,GC_content,...
    now_str,ind_species,n_hairpin_plot_f,n_hairpin_plot_r,f_pass_random);


% returning to starting directory
cd(cwd);