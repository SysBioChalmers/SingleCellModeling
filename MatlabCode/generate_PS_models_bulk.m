%Generates the models for the bulk data (the 8 T cell samples from the cell type profiles paper) with ftINIT. For Fig. 2A.
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

T = readtable('../data/PoolSize/scaledTMMMatrix.txt');


genes = T{:,1};
dataMat = table2array(T(:,2:end));
size(dataMat)%looks ok
length(genes)%looks ok
dataMat(1:10,1:8)


nModels = size(dataMat,2);

%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('../data/prepDataHumanGEM.mat');
%prepDataHumanGEM

model_indx = 1:nModels; %length(s.sampleIds); %run with 15 samples for now

paramsNewAlg = struct();
paramsNewAlg.TimeLimit = 120;
paramsNewAlg.MIPGap = 0.0004;

milpSkipMets.simpleMets.mets = {'H2O';'Pi';'PPi';'H+';'O2';'CO2';'Na+'};
milpSkipMets.simpleMets.compsToKeep = {'i'};


arrayData = struct(); 
arrayData.genes = genes;
arrayData.levels = dataMat;
arrayData.tissues = arrayfun(@(x) num2str(x), 1:nModels, 'UniformOutput', false);
arrayData.threshold = 1;

cd ../..

models = cell(nModels,1);

for i = 1:nModels
     if (isempty(models{i}))
         disp(['running model: ' num2str(i)])
         tic 
         mres = getINITModel9(prepDataHumanGEM,arrayData.tissues{i},[],[],arrayData,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,milpSkipMets,true,false,paramsNewAlg);
         toc
         mres.id = arrayData.tissues{i};
         models{i,1} = mres;
     end
end

% save results
save('../data/PoolSize/psBulkModels.mat','models');

