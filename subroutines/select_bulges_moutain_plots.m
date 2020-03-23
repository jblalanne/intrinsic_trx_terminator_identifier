function bulges = select_bulges_moutain_plots(data)

diff_data = diff(~diff(data));

if isempty(diff_data)
    bulges = [];
else
    ind_start = 1+find(diff_data==1);
    ind_stop = 1+find(diff_data==-1);
    
%     try
    if ~isempty(ind_start)
        if ind_start(1)>ind_stop(1)
            if size(ind_start,1)>1
%             ind_start = [1; ind_start];
                ind_start = [1; ind_start];
            else
                ind_start = [1 ind_start];
            end
        end
    end
%     catch
%        bla 
%     end
    bulges = ind_stop-ind_start;
    
end
end

% need to be careful about the beginning and end
% contiguous = [];
% for i = 1:n_contig
%     contiguous{i} = ind_start(i):ind_stop;
% end



