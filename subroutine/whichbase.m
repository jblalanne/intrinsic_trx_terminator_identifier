function base_number = whichbase(base)
if strcmp(base,'A')
    base_number = 1;
elseif strcmp(base,'T')
    base_number = 2;
elseif strcmp(base,'C')
    base_number = 3;
elseif strcmp(base,'G')
    base_number = 4;
else   % dealing with N's
    base_number = NaN;
end
end