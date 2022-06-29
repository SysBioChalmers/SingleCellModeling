function eGenes = getTaskEssentialGenes(models_file)
% Identify genes essential for different tasks in different models.
% This code is largely copied from the Human1 paper
%
% Input:
%
%   models_file  Path of .mat file containing a cell array of tINIT models.
%                If "models" is a number, then it will search the current
%                directory for a file named "init_models-#.mat".
%
%   task_file    Path of .mat file containing the metabolic tasks formatted
%                as a structure (i.e., generated using the parseTaskList
%                function).
%                Do NOT specify the path of a .xlsx file, as it is better
%                to avoid the use of java on the cluster (required for
%                parsing the Excel task file format).
%
%
%
% Output:
%
%   eGenes          results structure with the following fields:
%       taskList    list of metabolic tasks that were tested
%       tissues     list of tissues (model IDs) corresponding to each model
%       geneList    cell array of the list of genes from each model
%       essentialGenes   cell array with one entry per model, where each
%                        entry is a logical matrix with rows corresponding
%                        to genes (in geneList) and columns to tasks (in
%                        taskList). Entries in the matrix are true when a
%                        gene is essential for a task, and false otherwise.
%
%
% Usage:
%
%   eGenes = getTaskEssentialGenes(INIT_output, refModel, taskStruct);
%

% verify that the cluster node has sufficient CPUs available
x = parcluster('local');
if x.NumWorkers < 8
    error('Insufficient CPUs available (%u) for parallel processing!', x.NumWorkers); %for some reason, Vera often gives 8 cores, which is enough.
end

%for debugging:
%models_file = '../data/init_models_depmap15.mat';

% load models and metabolic tasks
if isnumeric(models_file)
    models_file = strcat('init_models-', num2str(models_file), '.mat');
end
x = load(models_file);
fns = fieldnames(x);
models = getfield(x,fns{1});

% specify paths for cluster use
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');

taskStruct = parseTaskList('metabolicTasks_Essential.txt');


% initialize variables
n = numel(models);
tissues = cell(n,1);
geneList = cell(n,1);
allEssentials = cell(n,1);

% Convert refModel genes from Ensembl IDs to gene abbreviations.
% This is done in a regular for-loop because it fails in parfor-loop.
for i = 1:length(models)
    [grRules,genes,rxnGeneMat] = translateGrRules(models{i}.grRules, 'Name', 'ENSG');
    models{i}.grRules = grRules;
    models{i}.genes = genes;
    models{i}.rxnGeneMat = rxnGeneMat;
    %also add boundary mets
    models{i} = closeModel(models{i});
end

% iterate through models in parallel for-loop
%parpool(8);  % assume 8 CPUs
parfor i = 1:length(models)
    % determine essential genes for each task
    [~,essentialGenes] = checkTasksGenes(models{i},[],false,false,true,taskStruct);
    
    % collect results
    tissues{i,1} = models{i}.id;
    geneList{i,1} = models{i}.genes;
    allEssentials{i,1} = essentialGenes;
    
end

% gather results into eGenes structure
eGenes = {};
eGenes.taskList = {taskStruct(:).description}';
eGenes.tissues = tissues;
eGenes.geneList = geneList;
eGenes.essentialGenes = allEssentials;

% save eGenes structure
out_filename = strcat(regexprep(models_file, '.mat', ''), '-eGenes.mat');
save(out_filename, 'eGenes');


