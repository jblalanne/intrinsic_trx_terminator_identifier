function [genome, full_genome_descriptor,short_genome_descriptor] = get_genomes_fna(genome_dir,file_name)

% directory movement
cwd = pwd;
cd(genome_dir)

% opening file
fid = fopen(file_name);

% initialization
genome = [];
full_genome_descriptor = [];
short_genome_descriptor = [];
counter = 1;

% reading the file (assumes that first line if .fasta format (> etc.)). 
tline = fgetl(fid);
full_genome_descriptor{counter} = tline(2:end);
ind = find(tline==' ',1,'first');
short_genome_descriptor{counter} = tline(2:ind-1);


L_pre_allocation = 2E7;     % largest bacterial genome 13M
genome{1} = blanks(L_pre_allocation);
current_ind = 1; 

% reading until end of file
while ~feof(fid)
    tline = fgetl(fid);
    if ~isempty(tline)
        if strcmp(tline(1),'>')     % genome element (assuming marked by fasta like >).
            
            counter = counter+1;                        % updating the counter (for DNA elements: chromosome, plasmids, etc.).
            full_genome_descriptor{counter} = tline(2:end);       % get the element name
            ind = find(tline==' ',1,'first');
            short_genome_descriptor{counter} = tline(2:ind-1);
            
            genome{counter} = blanks(L_pre_allocation);  % initialization of genome variable
            current_ind = 1;                            % position index variable
            
        else    % sequence data
            current_range = current_ind:(current_ind+length(tline)-1);  % range in the genome with the current line entry
            genome{counter}(current_range) = tline;                     % filling in the genome sequence
            current_ind = current_ind + length(tline);                  % updating position index variable
            
        end
    end

end

% remove the unnecessary pre-allocated space
for i = 1:length(genome)
   genome{i}(genome{i}==' ') = []; 
end

fclose(fid);    % housekeeping
cd(cwd)
