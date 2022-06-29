%Prep code for the old tINIT version
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode';

%load Human-GEM and tasks
load('C:/Work/MatlabCode/components/human-GEM/Human-GEM_1_12/Human-GEM/model/Human-GEM.mat')
%remove drug reactions from the model - we don't need them
model = removeDrugReactions(ihuman);
model = removeUnwantedReactions(model);
%remove all AA triplet rxns, i.e. rxns of the type
%2 H2O[c] + Tryptophanyl-Glycyl-Aspartate[c] <=> aspartate[c] + glycine[c] + tryptophan[c]
%This is only used for import of such compounds, and if you don't have that in your "medium",
%which we usually don't, these are pretty pointless.
AATriplets = getAATripletReactions(model,false);
%length(AATriplets)%735
model = removeReactions(model, AATriplets);
%Remove the duplicate complex I reaction with ROS - it only causes problems to include it
model = removeReactions(model, {'MAR13081'});


taskStruct = parseTaskList('C:/Work/MatlabCode/components/human-GEM/Human-GEM_1_12/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt');

model = addBoundaryMets(model);


% Run some preliminary steps that will allow us to skip some pre-processing
% steps in the tINIT algorithm, greatly reducing the overall run time.
[~,deletedDeadEndRxns] = simplifyModel(model,true,false,true,true,true);
cModel = removeReactions(model,deletedDeadEndRxns,false,true);
[taskReport, essentialRxnMat] = checkTasks(cModel,[],true,false,true,taskStruct);

% add pre-processing results to arrayData structure
arrayData = struct();
arrayData.deletedDeadEndRxns = deletedDeadEndRxns;
arrayData.taskReport = taskReport;
arrayData.essentialRxnMat = essentialRxnMat;
arrayData.threshold = 1;  % default value


save('../data/tINIT_inputs_human.mat','model','arrayData','taskStruct');


