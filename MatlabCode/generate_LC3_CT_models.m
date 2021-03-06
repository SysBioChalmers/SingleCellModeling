%Generates the models from the full LC3 cell types, i.e., without bootstrapping. These are used in Fig 2 D-F
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'


T = readtable('../data/LC3/ModelDataLC3.txt');


genes = T{:,1};
dataMat = table2array(T(:,2:end));
%size(dataMat)%looks ok
%length(genes)%looks ok
%dataMat(1:10,1:10)


nModels = size(dataMat,2);

%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('../data/prepDataHumanGEM.mat');
%prepDataHumanGEM

model_indx = 1:nModels; 

arrayData = struct(); 
arrayData.genes = genes;
arrayData.levels = dataMat;
arrayData.tissues = arrayfun(@(x) num2str(x), 1:nModels, 'UniformOutput', false);
arrayData.threshold = 1;
%CPM the data
arrayData.levels = arrayData.levels.*10^6./sum(arrayData.levels,1);
%sum(arrayData.levels,1)

models = cell(nModels,1);

for i = 1:nModels
     if (isempty(models{i}))
         disp(['running model: ' num2str(i)])
         tic %we run this without step 3
         mres = ftINIT(prepDataHumanGEM,arrayData.tissues{i},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);
         toc
         mres.id = arrayData.tissues{i};
         models{i,1} = mres;
     end
end

% save results
save('../data/LC3/pooledLC3Models.mat','models');

