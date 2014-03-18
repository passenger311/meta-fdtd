function var_inc = inc_type(var_inc)

% Create increment type strings

inc_str1    =   {'linear'};
inc_str2    =   {'log','logarithmic'};
inc_str3    =   {'custom'};

% Check increment type name against above examples

inc(1) = max(strcmpi(var_inc,inc_str1));
inc(2) = max(strcmpi(var_inc,inc_str2));
inc(3) = max(strcmpi(var_inc,inc_str3));

ind = find(inc == 1);

if size(ind) > 1 | isempty(ind) == 1 %#ok<OR2>
    fprintf(['Error: increment type not recognised, please see function '...
        ,'"var.m" for list of increment type names\n\n']);
    return
else 
    var_inc = ind;
end