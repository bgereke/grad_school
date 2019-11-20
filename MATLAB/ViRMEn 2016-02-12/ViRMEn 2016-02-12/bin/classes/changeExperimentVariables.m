function changeExperimentVariables(src,evt) %#ok<INUSL>

obj = evt.AffectedObject;

% Determine which variables got updated
fldold = fieldnames(obj.backedUpVariables);
fldnew = fieldnames(obj.variables);
fldboth = intersect(fldold,fldnew);
upd = setdiff(union(fldold,fldnew),fldboth);
for ndx = 1:length(fldboth)
    if ~strcmp(obj.variables.(fldboth{ndx}),obj.backedUpVariables.(fldboth{ndx}))
        upd{end+1} = fldboth{ndx}; %#ok<AGROW>
    end
end

if isempty(upd)
    return
end

dsc = obj.descendants;
for ndx = 1:length(dsc)
    props = fieldnames(dsc{ndx}.symbolic);
    for p = 1:length(props)
        hasVar = false;
        if ischar(dsc{ndx}.symbolic.(props{p}))
            for v = 1:length(upd)
                if ~isempty(regexp([' ' dsc{ndx}.symbolic.(props{p}) ' '],['[^A-Za-z0-9_]' upd{v} '[^A-Za-z0-9_]'],'once'))
                    hasVar = true;
                end
            end
        else
            for v = 1:length(upd)
                if any(cellfun(@(x)~isempty(regexp([' ' x ' '],['[^A-Za-z0-9_]' upd{v} '[^A-Za-z0-9_]'],'once')),dsc{ndx}.symbolic.(props{p})))
                    hasVar = true;
                end
            end
        end
        if hasVar
            dsc{ndx}.(props{p}) = dsc{ndx}.symbolic.(props{p});
        end
    end
end