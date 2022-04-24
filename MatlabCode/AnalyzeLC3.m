%This code assembles the LC3 bootstrap models' reaction contents and their task results into two files.

cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

allTasks = parseTaskList('../data/metabolicTasks_Full.txt');
filt = [allTasks.shouldFail] == 0;
filtTasks = allTasks(filt);

sampleIds5 = {'N_Alveolar Mac';...
    'N_AT2';...
    'N_CD4+ Th';...
    'N_Cytotoxic CD8+ T';...
    'N_Monocytes';...
    'N_NK';...
    'T_CD4+ Th';...
    'T_CD8+_CD4+ Mixed Th';...
    'T_Exhausted CD8+ T';...
    'T_Follicular B cells';...
    'T_MAST';...
    'T_mo-Mac';...
    'T_Naive CD4+ T';...
    'T_Treg';...
    'T_tS1';...
    'T_tS2'};

trs_5 = nan(sum(filt),length(sampleIds5));
for i = 1:length(sampleIds5)
    trs_5(:,i) = extractTaskResults(['../data/LC3/Run5/Bootstrap_' sampleIds5{i} '_taskResults.mat'], filt);
end

t_5 = struct();
t_5.data = trs_5;
t_5.tasks = {filtTasks.description}.';
t_5.sampleIds = sampleIds5;
save('../data/LC3/LC3_tasks_run_5.mat', 't_5')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load base model to get all possible rxns
x = load('../data/prepDataHumanGEM.mat');
baseModel = x.prepDataHumanGEM.refModel;
scmp_5 = nan(length(baseModel.rxns),length(sampleIds5));

%These take a while to run, a few minutes 
for i = 1:length(sampleIds5)
    scmp_5(:,i) = getModelStructureData(['../data/LC3/Run5/Bootstrap_' sampleIds5{i} '_models.mat'], baseModel);
end

%run tsne (not used)
proj_coords = tsne(double(scmp_5.'),'Distance','hamming','NumDimensions',2,'Exaggeration',20,'Perplexity',3);

s_5 = struct();
s_5.tsneX = proj_coords(:,1);
s_5.tsneY = proj_coords(:,2);
s_5.data = scmp_5;
s_5.rxns = baseModel.rxns;
s_5.sampleIds = sampleIds5;
s_5.subSys = baseModel.subSystems;
save('../data/LC3/LC3_rxns_run_5.mat', 's_5')

