% accumulate object within CC structure (used in accumulate structure)

for i=1:size(voi,1)
    if exist(voi{i,1},'var')
        if length(object)>0
            options.(sprintf(object)).(sprintf(voi{i,1}))=eval(voi{i,1});
        else
            options.(sprintf(voi{i,1}))=eval(voi{i,1});
        end
        
    elseif length(voi{i,2})>0
        
        if length(object)>0
            options.(sprintf(object)).(sprintf(voi{i,1}))=voi{i,2};
        else
            options.(sprintf(voi{i,1}))=voi{i,2};
        end
        
        evalc([sprintf(voi{i,1}) '=voi{i,2}']);

    end
end