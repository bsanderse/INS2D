% function options = accumulate_object(object,voi,options)
% accumulate object within options structure (used in accumulate structure)
% add 'object' as field to struct 'options' with values voi

for i=1:size(voi,1)
    if exist(voi{i,1},'var')
        % check if this variable is a variable in the workspace, and if so
        % store its value from the workspace into the struct (ignoring any
        % value given in the accumulate_structure file)
        if ~isempty(object)
            % options.object.variable = eval(variable)
            options.(sprintf(object)).(sprintf(voi{i,1}))=eval(voi{i,1});
        else
            % options.object = eval(variable)
            options.(sprintf(voi{i,1}))=eval(voi{i,1});
        end
        
    else
        % not yet existing in current workspace
        
        % check if length>0
        if ~isempty(voi{i,2})
        
            % create new variable in structure and give it the default
            % value as provided in accumulate_structure
            if ~isempty(object)
                options.(sprintf(object)).(sprintf(voi{i,1}))=voi{i,2};
            else
                options.(sprintf(voi{i,1}))=voi{i,2};
            end
            % put the variable also in the current workspace
            % (this can be dangerous if the same name is used in different
            % substructures)
%             evalc([sprintf(voi{i,1}) '=voi{i,2}']);
        end

    end
end