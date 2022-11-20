%Turns a cell array of models into a vector with pairwise jaccard values.
%Assumes an even number of models
function jacc = getBootstrapJaccardData(filename, baseModel)

x = load(filename);
models = x.models;

mRes = nan(length(baseModel.rxns),length(models));

%assemble into a matrix
for i = 1:length(models)
    mRes(:,i) = ismember(baseModel.rxns,models{i}.rxns);
end

jacc = nan(length(models)/2,1);
for i = 1:length(jacc)
    ind1 = i*2-1;
    ind2 = i*2;
    jacc(i) = jaccard(mRes(:,ind1), mRes(:,ind2));
end


end