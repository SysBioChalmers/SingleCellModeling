%This code assembles the MCOR3 bootstrap models' reaction contents and their task results into two files.

cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

allTasks = parseTaskList('../data/metabolicTasks_Full.txt');
filt = [allTasks.shouldFail] == 0;
filtTasks = allTasks(filt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Task evaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampleIds5 = {'L2_3 IT';...
    'L5 IT Pld5';...
    'L5 IT S100b';...
    'L5 IT Tcap_1';...
    'L5 IT Tcap_2';...
    'L5 NP Slc17a8_1';...
    'L5 NP Slc17a8_2';...
    'L6 CT Cpa6_1';...
    'L6 CT Cpa6_2';...
    'L6 CT Nxph2 Pou3f2_2';...
    'L6 IT Sulf1_2';...
    'L6 IT Sulf1_3';...
    'L6 IT Sulf1_4';...
    'L6 NP Trh_1';...
    'Lamp5 Pdlim5_2';...
    'Lamp5 Slc35d3_1';...
    'Vip Chat'};

trs_5 = nan(sum(filt),length(sampleIds5));
for i = 1:length(sampleIds5)
    trs_5(:,i) = extractTaskResults(['../data/MCOR3/Run5/Bootstrap_' sampleIds5{i} '_taskResults.mat'], filt);
end

t_5 = struct();
t_5.data = trs_5;
t_5.tasks = {filtTasks.description}.';
t_5.sampleIds = sampleIds5;
save('../data/MCOR3/MCOR3_tasks_run_5.mat', 't_5')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rxns comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load base model to get all possible rxns
x = load('../data/prepDataMouseGEM.mat');
baseModel = x.prepDataMouseGEM.refModel;
scmp_5 = nan(length(baseModel.rxns),length(sampleIds5));

%These take a while to run, a few minutes 
for i = 1:length(sampleIds5)
    scmp_5(:,i) = getModelStructureData(['../data/MCOR3/Run5/Bootstrap_' sampleIds5{i} '_models.mat'], baseModel);
end

%run tsne
proj_coords = tsne(double(scmp_5.'),'Distance','hamming','NumDimensions',2,'Exaggeration',20,'Perplexity',3);


%filt1 = sum(scmp - scmp(:,1) ~= 0,2) > 0;%selects all rxns that differ across samples

s_5 = struct();
s_5.tsneX = proj_coords(:,1);
s_5.tsneY = proj_coords(:,2);
s_5.data = scmp_5;
s_5.rxns = baseModel.rxns;
s_5.sampleIds = sampleIds5;
s_5.subSys = baseModel.subSystems;
save('../data/MCOR3/MCOR3_rxns_run_5.mat', 's_5')

%also extract data for supplementary fig about bootstrap variation
bootstrapJaccs = nan(50, length(sampleIds5));
for i = 1:length(sampleIds5)
    disp(num2str(i))
    bootstrapJaccs(:,i) = getBootstrapJaccardData(['../data/MCOR3/Run5/Bootstrap_' sampleIds5{i} '_models.mat'], baseModel);
end

s = struct();
s.sampleIds = sampleIds5;
s.bootstrapJaccs = bootstrapJaccs;

save('../data/MCOR3/MCOR3_bootstrap_jacc.mat', 's')


