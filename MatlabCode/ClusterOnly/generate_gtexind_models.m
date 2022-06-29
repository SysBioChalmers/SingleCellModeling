%Generates models for GTEx in TPM, the new tINIT version
function generate_gtexind_models(chunk, threshold)

if nargin < 2
	threshold = 1;
end

%for debugging locally
%cd C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/data

% add paths (comment out if you run locally)
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');



T = readtable('GTExInd/gtexIndSamp.txt');


genes = T{:,1};
dataMat = table2array(T(:,3:end));
%size(dataMat)%looks ok
%length(genes)%looks ok
%dataMat(1:10,1:10)


nModels = size(dataMat,2);

%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('prepDataHumanGEMEns.mat');
%prepDataHumanGEM

model_indx = 1:nModels; %length(s.sampleIds); %run with 15 samples for now


arrayData = struct(); 
arrayData.genes = genes;
arrayData.levels = dataMat;
arrayData.tissues = arrayfun(@(x) num2str(x), 1:nModels, 'UniformOutput', false);
arrayData.threshold = threshold;
%CPM the data
arrayData.levels = arrayData.levels.*10^6./sum(arrayData.levels,1);
%sum(arrayData.levels,1)

%Now, handle the chunks
% define chunks into which data will be broken up for parallelization
nChunks = 10;
indx = round(linspace(0, numel(arrayData.tissues), nChunks+1))';
model_chunks = arrayfun(@(i) (indx(i)+1:indx(i+1))', (1:nChunks)', 'UniformOutput', false);
if strcmpi(chunk, 'test')
    model_indx = 1;
else
    model_indx = model_chunks{chunk};
end

% extract subset of array data for building models
arrayData.tissues = arrayData.tissues(model_indx);
arrayData.levels = arrayData.levels(:, model_indx);

nModels = length(model_indx);


models = cell(nModels,1);

for i = 1:nModels
     if (isempty(models{i}))
         disp(['running model: ' num2str(i)])
         tic %we run this without step 3
		 mres = ftINIT(prepDataHumanGEMEns,arrayData.tissues{i},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);
         toc
         mres.id = arrayData.tissues{i};
         models{i,1} = mres;
     end
end

% save results
if threshold == 1
	fn = ['GTExInd/chunk_' num2str(chunk) '.mat'];
else
	fn = ['GTExInd_th' num2str(threshold) '/chunk_' num2str(chunk) '.mat'];
end
save(fn,'models');

