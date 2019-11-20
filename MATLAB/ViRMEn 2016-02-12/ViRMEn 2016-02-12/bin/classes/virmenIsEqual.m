function x = virmenIsEqual(a,b)

exclude = {'name','parent','items'};

propsA = unique(properties(a));
propsB = unique(properties(b));
if ~isequal(propsA,propsB)
    x = false;
elseif iscell(a) && iscell(b)
    if ~isequal(size(a),size(b))
        x = false;
    else
        x = false(1,numel(a));
        for ndx = 1:length(x)
            x(ndx) = virmenIsEqual(a{ndx},b{ndx});
        end
        x = all(x);
    end
elseif isempty(propsA) && isempty(propsB)
    x = isequalwithequalnans(a,b);
elseif ischar(a) && ischar(b)
    x = strcmp(a,b);
else
    propsA = setdiff(propsA,exclude);
    propsB = setdiff(propsB,exclude);
    x = false(1,length(propsA));
    for ndx = 1:length(propsA)
        x(ndx) = virmenIsEqual(a.(propsA{ndx}),b.(propsB{ndx}));
    end
    x = all(x);
end