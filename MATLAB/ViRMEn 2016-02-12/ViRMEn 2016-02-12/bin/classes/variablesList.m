function virmenVariablesList = variablesList(virmenObj)

virmenVariablesList = struct;
virmenItems = virmenObj.descendants;
for virmenNdx = 1:length(virmenItems)
    virmenProps = fieldnames(virmenItems{virmenNdx}.symbolic);
    for virmenP = 1:length(virmenProps)
        virmenStr = virmenItems{virmenNdx}.symbolic.(virmenProps{virmenP});
        if ~iscell(virmenStr)
            virmenStr = {virmenStr};
        end
        for virmenS = 1:length(virmenStr)
            try
                eval([virmenStr{virmenS} ';']);
            catch virmenME
                virmenF = strfind(virmenME.message,'''');
                virmenVarName = virmenME.message(virmenF(1)+1:virmenF(2)-1);
                virmenVariablesList.(virmenVarName) = virmenItems{virmenNdx}.ancestor.variables.(virmenVarName);
            end
        end
    end
end