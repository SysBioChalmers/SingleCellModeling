%Assembles structural model data from both the old and new tINIT
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old tINIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%


oldtINITGTExModels = {};

for i = 1:20
    disp(num2str(i))
    fn = ['../data/oldInitModelsGTEx/init_models-' num2str(i) '.mat'];
    x = load(fn);
    oldtINITGTExModels = [oldtINITGTExModels;x.init_models];
end

%length(oldtINITGTExModels) %265, ok

load('../data/prepDataHumanGEM.mat');

baseModel = prepDataHumanGEM.refModel;

%now build a matrix
allMod = oldtINITGTExModels;

compMat = false(length(baseModel.rxns), length(allMod));



for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,allMod{i}.rxns);
end

rng(1);%random seed
proj_coords = tsne(double(compMat.'),'Distance','hamming','NumDimensions',2,'Exaggeration',6,'Perplexity',10);

%read the gtex tissue data
T = readtable('../data/GTExInd/gtexIndSampTissues.txt', 'ReadVariableNames',false);
gtexTissues = table2cell(T(:,1));

%export to R
d = struct();
d.tsneX = proj_coords(:,1);
d.tsneY = proj_coords(:,2);
d.tissues = gtexTissues;
save('../data/tsneGTExOldtINIT.mat', 'd');

clear oldtINITGTExModels;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Also process the GTEx from the new model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gtexIndModels = {};

for i = 1:10
    disp(num2str(i))
    fn = ['../data/GTExInd/chunk_' num2str(i) '.mat'];
    x = load(fn);
    gtexIndModels = [gtexIndModels;x.models];
end

%length(gtexIndModels) %265, ok

load('../data/prepDataHumanGEM.mat');

baseModel = prepDataHumanGEM.refModel;

%now build a matrix
allMod = gtexIndModels;

compMat = false(length(baseModel.rxns), length(allMod));



for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,allMod{i}.rxns);
end

rng(1);%random seed
proj_coords = tsne(double(compMat.'),'Distance','hamming','NumDimensions',2,'Exaggeration',6,'Perplexity',10);

%read the gtex tissue data
T = readtable('../data/GTExInd/gtexIndSampTissues.txt', 'ReadVariableNames',false);
gtexTissues = table2cell(T(:,1));

%export to R
d = struct();
d.tsneX = proj_coords(:,1);
d.tsneY = proj_coords(:,2);
d.tissues = gtexTissues;
save('../data/tsneGTExNewtINIT.mat', 'd');

clear gtexIndModels;

