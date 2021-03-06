%Runs task analysis on all bootstraps for one cell type (specificed by chunk) from the LC3 dataset
function [] = TaskAnalysisLC3OnCluster(chunk)

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


inFilename = ['LC3/Run5/Bootstrap_' sampleNames{chunk} '_models.mat'];
outFilename = ['LC3/Run5/Bootstrap_' sampleNames{chunk} '_taskResults.mat'];


% add paths
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi'); %only need to be done once, may cause threading problems!

x = load(inFilename);
models = x.models;


taskReports = cell(length(models),1);
allTasks = parseTaskList('metabolicTasks_Full.txt');

%parpool(8);
parfor i = 1:length(models)
    disp(['Running: ' num2str(i)])
    [taskReports{i}, ~]=checkTasks(closeModel(models{i}),[],true,false,false,allTasks);
end

save(outFilename, 'taskReports')

