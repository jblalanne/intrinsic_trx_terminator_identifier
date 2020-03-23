function MFE = read_RNAfold_results(n_seq)

% reading the corresponding output from RNAfold
files = dir('results_*');

MFE = zeros(length(files),n_seq);

for i = 1:length(files)
    
    fid = fopen(files(i).name);
    
    tline = fgetl(fid);
    counter = 1;
    counter2 = 1;
    while ischar(tline)
        tline = fgetl(fid);
        
        if mod(counter,6)==2
            MFE(i,counter2)=sscanf(tline,'%*s (%f)');
            counter2=counter2+1;
%         elseif mod(counter,6)==5
%             prob_MFE_struct(i,counter2)=sscanf(tline,'%*s %*s %*s %*s %*s %*s %f');
        end
        
        counter = counter+1;
    end
    
    fclose(fid);
end