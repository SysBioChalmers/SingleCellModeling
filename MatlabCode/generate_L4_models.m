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
mres = ftINIT(prepDataHumanGEM,arrayData.tissues{1},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);
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
mres = ftINIT(prepDataHumanGEM,arrayData.tissues{1},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);
toc
mres.id = arrayData.tissues{1};
models{2,1} = mres;


% save results
save('../data/L4/L4Models.mat','models');

