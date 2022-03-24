function dot_bracket = mountain_var_to_dot_bracket(mountain_var)

dot_bracket = [];
diff_mv = [0 diff(mountain_var)];

for i =1:length(diff_mv)
    if diff_mv(i)==0
        dot_bracket = [dot_bracket '.'];
    elseif diff_mv(i)==1
        dot_bracket = [dot_bracket '('];
    elseif diff_mv(i)==-1
        dot_bracket = [dot_bracket ')'];
    end
end

