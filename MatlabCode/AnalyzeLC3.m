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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investigate which pathway is used
% for odd-chain fatty acid synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'
models = load('../data/LC3/Run5/Bootstrap_T_tS2_models.mat').models;
m = models{1}; %all models can generate odd-chain fatty acids, so just grab one
[inExchRxns,inExchRxnsInd] = getExchangeRxns(m);
constructEquations(m, inExchRxns)
m.lb(inExchRxnsInd) = 0; %stops import, not export

%look at henicosanoic acid[c]
%inputs are according to the task in the task file, all output was allowed
m.lb(strcmp(m.rxns, 'MAR09048')) = -1000; %O2
m.lb(strcmp(m.rxns, 'MAR09034')) = -1000; %glucose

rxnsToAdd = {};
rxnsToAdd.rxns = {'isoleucine_imp', 'faexp'};
rxnsToAdd.equations = {'=> isoleucine[c]','henicosanoic acid[c] =>'};
rxnsToAdd.lb = [-1000,0];
rxnsToAdd.ub = [0,1000];
m2 = addRxns(m, rxnsToAdd, 3);
m2.c = double(strcmp(m2.rxns, 'faexp'));

%run simulation
res = solveLP(m2,1);
res %-52.5626 - looks good

listMetRxnsWithFluxes(m2, res, 'henicosanoic acid[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, 'heneicosanoyl-CoA[c]', false, 10^-3) %mostly from (2E)-heneicosenoyl-CoA
listMetRxnsWithFluxes(m2, res, '(2E)-heneicosenoyl-CoA[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '3-hydroxyheneicosanoyl-CoA[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '3-oxoheneicosanoyl-CoA[c]', false, 10^-3)%malonyl-CoA[c] + nonadecanoyl-CoA[c]
listMetRxnsWithFluxes(m2, res, 'malonyl-CoA[c]', false, 10^-3) %acetyl-CoA[c] - end of that branch
listMetRxnsWithFluxes(m2, res, 'nonadecanoyl-CoA[c]', false, 10^-3) %mostly from (2E)-nonadecenoyl-CoA[c]
listMetRxnsWithFluxes(m2, res, '(2E)-nonadecenoyl-CoA[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '3-hydroxynonadecanoyl-CoA[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '3-oxononadecanoyl-CoA[c]', false, 10^-3)%malonyl-CoA[c] + heptadecanoyl-CoA[c]
listMetRxnsWithFluxes(m2, res, 'heptadecanoyl-CoA[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, 'margaric acid[c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, 'heptadecanoyl-[ACP][c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '(2E)-heptadecenoyl-[ACP][c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '3-hydroxyheptadecanoyl-[ACP][c]', false, 10^-3)
listMetRxnsWithFluxes(m2, res, '3-oxoheptadecanoyl-[ACP][c]', false, 10^-3) %malonyl-[ACP][c] + pentadecanoyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'pentadecanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '(2E)-pentadecenoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-hydroxypentadecanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-oxopentadecanoyl-[ACP][c]', false, 10^-3) %malonyl-[ACP][c] + tridecanoyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'tridecanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '(2E)-tridecenoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-hydroxytridecanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-oxotridecanoyl-[ACP][c]', false, 10^-3) %malonyl-[ACP][c] + undecanoyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'undecanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '(2E)-undecenoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-hydroxyundecanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-oxoundecanoyl-[ACP][c]', false, 10^-3) % malonyl-[ACP][c] + nonanoyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'nonanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '(2E)-nonenoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-hydroxynonanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-oxononanoyl-[ACP][c]', false, 10^-3) %heptanoyl-[ACP][c] + malonyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'heptanoyl-[ACP]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '(2E)-heptenoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-hydroxyheptanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-oxoheptanoyl-[ACP][c]', false, 10^-3) %malonyl-[ACP][c] + pentanoyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'pentanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '(2E)-pentenoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-hydroxypentanoyl-[ACP][c]', false, 10^-3) 
listMetRxnsWithFluxes(m2, res, '3-oxopentanoyl-[ACP][c]', false, 10^-3) %malonyl-[ACP][c] + propanoyl-[ACP][c]
listMetRxnsWithFluxes(m2, res, 'propanoyl-[ACP]', false, 10^-3) %{'[ACP][c] + propanoyl-CoA[c] => CoA[c] + propanoyl-[ACP][c]'}
%Conclusion: the odd-chain fatty acids are generated from propanoyl-CoA and acetyl-CoA

