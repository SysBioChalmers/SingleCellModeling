%Generates the models for the L4 dataset using ftINIT
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'


T = readtable('../data/L4/L4_all_pooled_lung.txt');


genes = T{:,1};
dataMat = table2array(T(:,2:end));
%size(dataMat)%looks ok
%length(genes)%looks ok
%dataMat(1:10,1:10)

%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('../data/prepDataHumanGEM.mat');
%prepDataHumanGEM

paramsNewAlg = struct();
paramsNewAlg.TimeLimit = 120;
paramsNewAlg.MIPGap = 0.0004;

milpSkipMets.simpleMets.mets = {'H2O';'Pi';'PPi';'H+';'O2';'CO2';'Na+'};
milpSkipMets.simpleMets.compsToKeep = {'i'};


arrayData = struct(); 
arrayData.genes = genes;
arrayData.levels = dataMat;
arrayData.tissues = {'LT lung'};
arrayData.threshold = 1;
%CPM the data
arrayData.levels = arrayData.levels.*10^6./sum(arrayData.levels,1);
%sum(arrayData.levels,1)

models = cell(2,1);

tic %we run this without step 3
mres = getINITModel9(prepDataHumanGEM,arrayData.tissues{1},[],[],arrayData,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,milpSkipMets,true,false,paramsNewAlg);
toc
mres.id = arrayData.tissues{1};
models{1,1} = mres;


%Now spleen
T = readtable('../data/L4/L4_all_pooled_spleen.txt');


genes = T{:,1};
dataMat = table2array(T(:,2:end));

arrayData = struct(); 
arrayData.genes = genes;
arrayData.levels = dataMat;
arrayData.tissues = {'LT Spleen'};
arrayData.threshold = 1;
%CPM the data
arrayData.levels = arrayData.levels.*10^6./sum(arrayData.levels,1);


tic %we run this without step 3
mres = getINITModel9(prepDataHumanGEM,arrayData.tissues{1},[],[],arrayData,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,milpSkipMets,true,false,paramsNewAlg);
toc
mres.id = arrayData.tissues{1};
models{2,1} = mres;


% save results
save('../data/L4/L4Models.mat','models');

