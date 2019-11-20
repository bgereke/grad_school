function backupVariableValues(src,evt) %#ok<INUSL>

obj = evt.AffectedObject;
obj.backedUpVariables = obj.variables;