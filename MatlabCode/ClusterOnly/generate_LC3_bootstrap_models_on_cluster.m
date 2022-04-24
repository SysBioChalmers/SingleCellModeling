%Generates all the bootstrap models for one celltype, specified by chunk, for the LC3 dataset. 
%Designed to be called from a cluster script.
function [] = generate_LC3_bootstrap_models_on_cluster(chunk)

sampleNames = {'N_Alveolar Mac'; ...
               'N_AT2'; ...
               'N_CD4+ Th'; ...
               'N_Cytotoxic CD8+ T'; ...
               'N_Monocytes'; ...
               'N_NK'; ...
               'T_CD4+ Th'; ...
               'T_CD8+_CD4+ Mixed Th'; ...
               'T_Exhausted CD8+ T'; ...
               'T_Follicular B cells'; ...
               'T_MAST'; ...
               'T_mo-Mac'; ...
               'T_Naive CD4+ T'; ...
               'T_Treg'; ...
               'T_tS1'; ...
               'T_tS2'};


inFilename = ['LC3Data5/Bootstrap_' sampleNames{chunk} '.txt'];
outFilename = ['LC3Data5/Bootstrap_' sampleNames{chunk} '_models.mat'];

% add paths
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');



T = readtable(inFilename);
genes = T{:,1};
dataMat = table2array(T(:,2:end));

nModels = size(dataMat,2);

%% run tINIT

load('prepDataHumanGEM.mat');

model_indx = 1:nModels; 

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
%CPM the data
arrayData.levels = arrayData.levels.*10^6./sum(arrayData.levels,1);
%sum(arrayData.levels,1)

models = cell(nModels,1);

for i = 1:nModels
    % First try to run tINIT with shorter time limit. If it fails, then
    % try again with a longer time limit (and then again).
     disp(['running model: ' num2str(i)])
     tic %we run this without step 3
     mres = getINITModel9(prepDataHumanGEM,arrayData.tissues{i},[],[],arrayData,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,milpSkipMets,true,false,paramsNewAlg);
     toc
     mres.id = arrayData.tissues{i};
     models{i,1} = mres;
end

% save results
save(outFilename,'models');

