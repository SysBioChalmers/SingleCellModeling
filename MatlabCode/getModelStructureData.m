%Turns a cell array of models into a matrix with reaction presence
function msd = getModelStructureData(filename, baseModel)


x = load(filename);
models = x.models;

mRes = nan(length(baseModel.rxns),length(models));

%assemble into a matrix
for i = 1:length(models)
    mRes(:,i) = ismember(baseModel.rxns,models{i}.rxns);
end

msd = sum(mRes,2);

end