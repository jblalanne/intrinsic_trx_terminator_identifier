function [mountain_var,structure_char] = get_mountain_plot_v2(lengths_folded)

files = dir('results_*');

mountain_var = zeros(length(files),length(lengths_folded),max(lengths_folded));
% index 1: sequence ID.
% index 2: folded length id. 
% index 3: positions.

structure_char = [];

for i = 1:length(files)
    
    fid = fopen(files(i).name);
    
    file_content = [];
    tline = fgetl(fid);
    
    counter = 1;
    counter2 = 1;
    while ischar(tline)
        if mod(counter,6)==3
            % getting the structure content
            ind_end = min(regexp(tline,' '));
            file_content{counter2} = tline(1:(ind_end-1));
            counter2=counter2+1;
        end
        tline = fgetl(fid);
        counter = counter+1;
    end
    
    % keeping only the nucleotides
    if length(file_content)==length(lengths_folded)
        for j = 1:length(lengths_folded)
            
            structure = file_content{j};     
           
            structure_num = zeros(lengths_folded(j),1);
            structure_num(structure=='(')=1;
            structure_num(structure==')')=-1;
            
            % padding with the appropriate number of zeros
            mountain_var(i,j,:) = [zeros(1,max(lengths_folded)-lengths_folded(j),1) cumsum(structure_num)'];
            structure_char{i,j} = structure;
        end
    end
    
    fclose(fid);
    
end

end


