function fillAncestry(src,evt) %#ok<INUSL>

obj = evt.AffectedObject;
ch = obj.children;
for ndx = 1:length(ch)
    ch{ndx}.parent = obj;
end