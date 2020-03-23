function reverse_complement_sequence = RCsequence(sequence)

% complement
ind_A = strfind(sequence,'A');
ind_T = strfind(sequence,'T');
ind_C = strfind(sequence,'C');
ind_G = strfind(sequence,'G');
% ind_others = 1:length(sequence);
% ind_all = sort([ind_A ind_G ind_C ind_T]);
% ind_others(ind_all) = [];


complement_sequence = cell(1);

% length(sequence{1})
% length(ind_A{1})+length(ind_T{1})+length(ind_C{1})+length(ind_G{1})
% % what is not A, T, C or G?
% not_found = ones(length(sequence{1}),1);
% not_found(ind_A{1})=0;
% not_found(ind_T{1})=0;
% not_found(ind_C{1})=0;
% not_found(ind_G{1})=0;
% ind_not_found = find(not_found==1);


complement_sequence{1}(ind_A) = 'T';
complement_sequence{1}(ind_T) = 'A';
complement_sequence{1}(ind_C) = 'G';
complement_sequence{1}(ind_G) = 'C';

reverse_complement_sequence = char(cellstr(fliplr(char(complement_sequence))));

end