function [positions_homo_rich_f, positions_homo_rich_r] = ...
    get_homonucleotide_subsequence_v3(genome,start,stop,strand,...
    bool_region,nucleotide_oi,...
    lower_n,upper_n,l_before_stop,l_after_stop,l_in_next_CDS)



L = length(genome);
RC_genome = RCsequence_v2(genome);

%% go from nucleotides to numbers to facilitate downstream processing
gen_num_for = zeros(length(genome),1);
gen_num_rev = zeros(length(RC_genome),1);

for i = 1:length(genome)
    gen_num_for(i) = whichbase(genome(i));
    gen_num_rev(i) = whichbase(RC_genome(i));
end



%% restrict attention to region of interest

% bool_region == 0: keep all
% bool_region == 1: keep CDS
% bool_region == 2: keep intergenic regions
% bool_region == 3: only downstream of genes with head genes downstream.
% bool_region == 4: downstream of genes with far co-directional genes downstream
% bool_region == 5: downstream of genes with close co-directional genes downstream

L_genome = length(genome);

if bool_region==1
    
    start_forward = start(strand==1);
    stop_forward = stop(strand==1);
    
    % use reverse strand convention for the reverse strand.
    start_reverse = L-stop(strand==0)+1;
    stop_reverse = L-start(strand==0)+1;
    
    indices_CDS_forward = [];
    for i = 1:length(start_forward)
        indices_CDS_forward = [indices_CDS_forward start_forward(i):stop_forward(i)];
    end
    
    indices_CDS_reverse = [];
    for i = 1:length(start_reverse)
        indices_CDS_reverse = [start_reverse(i):stop_reverse(i) indices_CDS_reverse ];
    end
    
    
    % replace non coding sequences by 0
    genome_num_for_2 = zeros(length(genome),1);
    genome_num_for_2(indices_CDS_forward) = gen_num_for(indices_CDS_forward);
    
    genome_num_rev_2 = zeros(length(genome),1);
    genome_num_rev_2(indices_CDS_reverse) = gen_num_rev(indices_CDS_reverse);
    
    disp('Done keeping only CDS.');
elseif bool_region==2
    
    indices_forward = [];
    start_forward = start(strand==1);
    stop_forward = stop(strand==1);
    for i = 1:length(start_forward)-1
        ind_oi = (stop_forward(i)-l_before_stop):...
            min([(start_forward(i+1)+l_in_next_CDS) (stop_forward(i)+l_after_stop)]);
        indices_forward = [indices_forward ind_oi];
    end
    
    
    % use reverse strand convention for the reverse strand.
    indices_reverse = [];
    start_reverse = fliplr(L-stop(strand==0)+1);
    stop_reverse = fliplr(L-start(strand==0)+1);
    
    for i = 1:length(start_reverse)-1
        ind_oi = (stop_reverse(i)-l_before_stop):...
            min([(start_reverse(i+1)+l_in_next_CDS) (stop_reverse(i)+l_after_stop)]);
        indices_reverse = [indices_reverse ind_oi];
    end
    
    
    genome_num_for_2 = zeros(length(genome),1);
    genome_num_for_2(indices_forward) = gen_num_for(indices_forward);
    
    genome_num_rev_2 = zeros(length(genome),1);
    genome_num_rev_2(indices_reverse) = gen_num_rev(indices_reverse);
    
elseif bool_region==0
    genome_num_for_2 = gen_num_for;
    genome_num_rev_2 = gen_num_rev;
    
    
elseif bool_region == 3         % non co-directional genes
    
    indices_forward = [];
    for i = 1:(length(stop)-1)
        if (strand(i)==1) && (strand(i+1)==0)
            upper_value = min([(stop(i)+l_after_stop) L_genome]);
            lower_value = max([1 (stop(i)-l_before_stop)]);
            ind_oi = lower_value:upper_value;
            indices_forward = [indices_forward ind_oi];
        end
    end
    
    indices_reverse = [];
    start_reverse = (L-stop+1);       % reverse position convention
    stop_reverse = (L-start+1);
    for i = 2:length(start_reverse)
        if (strand(i)==0) && (strand(i-1)==1)
            upper_value = min([L_genome (stop_reverse(i)+l_after_stop)]);
            lower_value = max([1 (stop_reverse(i)-l_before_stop)]);
            ind_oi = lower_value:upper_value;
            indices_reverse = [indices_reverse ind_oi];
        end
    end
    
    
    genome_num_for_2 = zeros(length(genome),1);
    genome_num_for_2(indices_forward) = gen_num_for(indices_forward);
    
    genome_num_rev_2 = zeros(length(genome),1);
    genome_num_rev_2(indices_reverse) = gen_num_rev(indices_reverse);
    
elseif bool_region == 4         % co-directional, but downstream gene far ahead
    
    indices_forward = [];
    for i = 1:(length(stop)-1)
        if (strand(i)==1) && (strand(i+1)==1) && ((start(i+1)-stop(i))>l_after_stop)
            upper_value = min([(stop(i)+l_after_stop) L_genome]);
            lower_value = max([1 (stop(i)-l_before_stop)]);
            ind_oi = lower_value:upper_value;
            indices_forward = [indices_forward ind_oi];
        end
    end
    
    indices_reverse = [];
    start_reverse = (L-stop+1);       % reverse position convention
    stop_reverse = (L-start+1);
    for i = 2:length(start_reverse)
        if (strand(i)==0) && (strand(i-1)==0) && ((start_reverse(i-1)-stop_reverse(i))>l_after_stop)
            upper_value = min([L_genome (stop_reverse(i)+l_after_stop)]);
            lower_value = max([1 (stop_reverse(i)-l_before_stop)]);
            ind_oi = lower_value:upper_value;
            indices_reverse = [indices_reverse ind_oi];
        end
    end
    
    genome_num_for_2 = zeros(length(genome),1);
    genome_num_for_2(indices_forward) = gen_num_for(indices_forward);
    
    genome_num_rev_2 = zeros(length(genome),1);
    genome_num_rev_2(indices_reverse) = gen_num_rev(indices_reverse);
    
elseif bool_region == 5     % co-directional, but downstream gene close
    
    indices_forward = [];
    for i = 1:(length(stop)-1)
        if (strand(i)==1) && (strand(i+1)==1) && ((start(i+1)-stop(i))<l_after_stop)
            upper_value = min([(start(i+1)+l_in_next_CDS) L_genome]);
            lower_value = max([1 (stop(i)-l_before_stop)]);
            ind_oi = lower_value:upper_value;
            indices_forward = [indices_forward ind_oi];
        end
    end
    
    indices_reverse = [];
    start_reverse = (L-stop+1);       % reverse position convention
    stop_reverse = (L-start+1);
    for i = 2:length(start_reverse)
        if (strand(i)==0) && (strand(i-1)==0) && ((start_reverse(i-1)-stop_reverse(i))<l_after_stop)
            upper_value = min([(start_reverse(i-1)+l_in_next_CDS) L_genome]);
            lower_value = max([1 (stop_reverse(i)-l_before_stop)]);
            ind_oi = lower_value:upper_value;
            indices_reverse = [indices_reverse ind_oi];
        end
    end
    
    genome_num_for_2 = zeros(length(genome),1);
    genome_num_for_2(indices_forward) = gen_num_for(indices_forward);
    
    genome_num_rev_2 = zeros(length(genome),1);
    genome_num_rev_2(indices_reverse) = gen_num_rev(indices_reverse);
    
end

%% U rich regions


% 1 if base U, 0 otherwise (within the region of interest narrowed down
% above).
nucleotide_oi_2_for = zeros(length(genome),1);
nucleotide_oi_2_for(genome_num_for_2==nucleotide_oi)=1;

nucleotide_oi_2_rev = zeros(length(genome),1);
nucleotide_oi_2_rev(genome_num_rev_2==nucleotide_oi)=1;

% performing traveling sum:
% half_width = 4;
% summing_filter = filter_generator(length(bsub_genome),half_width,'reverse');
% U_sequence_CDS_for = int8((2*half_width)*convolve(U_CDS_for,summing_filter));

width = lower_n:upper_n;

positions_homo_rich_f = [];
positions_homo_rich_r = [];


for i = 1:length(width)
    summing_filter = zeros(length(genome),1);
    x1 = ceil(length(genome)/2)-ceil(width(i)/2);
    x2 = x1+width(i)-1;
    summing_filter(x1:x2) = 1;
    homo_sequence_for = int8(convolve(nucleotide_oi_2_for,summing_filter));
    homo_sequence_rev = int8(convolve(nucleotide_oi_2_rev,summing_filter));
    
    homo_sequence_for_0 = homo_sequence_for(2:end-1);
    homo_sequence_for_p1 = homo_sequence_for(1:end-2);
    homo_sequence_for_m1 = homo_sequence_for(3:end);
    position_f = find(homo_sequence_for_0==width(i) & ...
        homo_sequence_for_p1==(width(i)-1) & homo_sequence_for_m1==(width(i)-1))+width(i)-2;
    
    homo_sequence_rev_0 = homo_sequence_rev(2:end-1);
    homo_sequence_rev_p1 = homo_sequence_rev(1:end-2);
    homo_sequence_rev_m1 = homo_sequence_rev(3:end);
    position_r = (L-find(homo_sequence_rev_0==width(i) & ...
        homo_sequence_rev_p1==(width(i)-1) & homo_sequence_rev_m1==(width(i)-1))+1)-width(i)+2;
    
    
    % correcting for the way the filter was set-up so that the reported
    % position is at the very end of the U-tract.
    if (width(i)==5 || width(i)==4)
        position_f = position_f + 2;
        position_r = position_r - 2;
    elseif (width(i)==6) || (width(i)==7)
        position_f = position_f + 1;
        position_r = position_r - 1;
    elseif (width(i)==10) || (width(i)==11)
        position_f = position_f - 1;
        position_r = position_r + 1;
    elseif (width(i)==12) || (width(i)==13)
        position_f = position_f - 2;
        position_r = position_r + 2;
    end
    positions_homo_rich_f{i} = position_f;
    positions_homo_rich_r{i} = position_r;
    
end

end