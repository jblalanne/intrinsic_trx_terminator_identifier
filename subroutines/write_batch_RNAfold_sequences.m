function write_batch_RNAfold_sequences(genome,lengths_folded,positions,...
constraint_bool,consecutive_Us)

n_sequences = length(positions);

% generating one file per positions. each file with subsequences from the
% 3' end, previously chosen lengths upstream. 
for i = 1:n_sequences
    
    % creating file
    if i<10
        csvwrite(sprintf('sequence_upstream_000%d.fas',i),'');
        fid = fopen(sprintf('sequence_upstream_000%d.fas',i),'w');
    elseif i<100
        csvwrite(sprintf('sequence_upstream_00%d.fas',i),'');
        fid = fopen(sprintf('sequence_upstream_00%d.fas',i),'w');
    elseif i<1000
        csvwrite(sprintf('sequence_upstream_0%d.fas',i),'');
        fid = fopen(sprintf('sequence_upstream_0%d.fas',i),'w');
    else
        csvwrite(sprintf('sequence_upstream_%d.fas',i),'');
        fid = fopen(sprintf('sequence_upstream_%d.fas',i),'w');
    end
    
    
    % printing sequences from the 3' ends extending upstream. 
    for j = 1:length(lengths_folded)
               
        fprintf(fid,sprintf('>sequence %d (3'' end %d), folding L=%d upstream\n',i, positions(i),lengths_folded(j)));
        lower_pos = max([1 (positions(i)-lengths_folded(j)+1)]);
        fprintf(fid,[genome(lower_pos:positions(i)) '\n']);
        
        counter = 0;
        constraint = '';
        while counter<lengths_folded(j)
            if counter<consecutive_Us(i)
                if constraint_bool
                    constraint = ['x' constraint];
                else
                    constraint = ['.' constraint];
                end
            else
                constraint = ['.' constraint];
            end
            counter = counter+1;
        end
        fprintf(fid,[constraint '\n']);
    end
    fclose(fid);
end

end