function [] = generate_GTEx_models_old_tINIT(chunk, restart)
% generates tINIT models from RNA-Seq profiles, using the old tINIT
%
% Note: This function was designed specifically for use on the cluster
%
% Input:
%
%   chunk     integer ranging from 1 to 10, or 'test' to run a test that
%             generates only one model.
%
%   restart   if TRUE, then re-start the model generation from a previous
%             incompleted run.
%             (Default = FALSE)
%

if nargin < 2
    restart = false;
end

% add paths
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');

%Note - I removed the first two lines in the file to simplify reading of the file
T = readtable('gtexIndSamp.txt');
genes = T{:,1};
dataMat = table2array(T(:,3:end));

sampleIds = arrayfun(@num2str, 1:size(dataMat,2), 'UniformOutput', false);



%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('tINIT_inputs_human.mat');

% define chunks into which data will be broken up for parallelization
nChunks = 20;
indx = round(linspace(0, size(dataMat,2), nChunks+1))';
model_chunks = arrayfun(@(i) (indx(i)+1:indx(i+1))', (1:nChunks)', 'UniformOutput', false);
if strcmpi(chunk, 'test')
    model_indx = 1;
else
    model_indx = model_chunks{chunk};
end

arrayData.genes = genes;
arrayData.tissues = sampleIds(model_indx);
arrayData.levels = dataMat(:,model_indx);


% initialize some variables
params = {};
init_models = {};
n_models = numel(arrayData.tissues);

% re-set to previous state if restart is requested
if ( restart )
    load(strcat('oldInitModelsGTEx/init_models-',num2str(chunk)));
    start_ind = numel(init_models) + 1;
else
    start_ind = 1;
end

for i = start_ind:n_models
    
    % First try to run tINIT with shorter time limit. If it fails, then
    % try again with a longer time limit (and then again).
    try
        params.TimeLimit = 1000;
        init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);
    catch
        try
            params.TimeLimit = 5000;
            init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);
        catch
            params.TimeLimit = 10000;
            init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);
        end
    end
    
    init_model.id = arrayData.tissues{i};
    init_models{i,1} = init_model;
    
    % save results
    if isnumeric(chunk)
        filename = strcat('oldInitModelsGTEx/init_models-',num2str(chunk));
    else
        filename = 'oldInitModelsGTEx/init_models-TEST';
    end
    save(filename,'init_models');
    
end



