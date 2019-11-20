function setNumericProperty(src,evt)

obj = evt.AffectedObject;
prop = src.Name;
val = obj.(prop);

if iscell(val) % cell array of strings
    obj.symbolic.(prop) = val;
elseif ischar(val) % symbolic string
    obj.symbolic.(prop) = val;
else % numeric array
    str = num2cell(val);
    str = cellfun(@num2str,str,'uniformoutput',false);
    obj.symbolic.(prop) = str;
end

str = obj.symbolic.(prop);
exper = obj.ancestor;
if iscell(str)
    val = zeros(size(str));
    for ndx = 1:numel(val)
        val(ndx) = symb2val(exper,str{ndx});
    end
else
    val = symb2val(exper,str);
end
obj.(prop) = val;


function virmenVal = symb2val(virmenExper,virmenString)

virmenIsBad = true;
while virmenIsBad
    virmenVars = fieldnames(virmenExper.variables);
    for virmenNdx = 1:length(virmenVars)
        eval([virmenVars{virmenNdx} '=eval(''' virmenExper.variables.(virmenVars{virmenNdx}) ''');']);
    end
    
    try
        virmenVal = eval(virmenString);
        virmenIsBad = false;
    catch virmenME
        if strcmp(virmenME.identifier,'MATLAB:UndefinedFunction')
            virmenF = strfind(virmenME.message,'''');
            virmenVarName = virmenME.message(virmenF(1)+1:virmenF(2)-1);
            virmenAnswer = inputdlg({['Enter a value for the new variable ''' virmenVarName '''.']},'New variable',1,{''});
            if isempty(virmenAnswer)
                virmenAnswer = {'NaN'};
            end
            virmenExper.variables.(virmenVarName) = virmenAnswer{1};
        else
            errordlg('Invalid expression.','Error');
            virmenVal = NaN;
            return
        end
    end
end