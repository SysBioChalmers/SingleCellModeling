function [] = generate_DepMap_models(chunk, restart)
% generate tINIT models from RNA-Seq profiles
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

%% run tINIT


%we replace the arraydata with the one from depmap
load('DepMap/arrayDataDepMap.mat')

load('tINIT_inputs_human.mat');

subsetToUse = 1:15;

%replace the arrayData with the DepMap array data
arrayData.genes = arrayDataDepMap.genes;
arrayData.tissues = arrayDataDepMap.tissues(subsetToUse);
arrayData.levels = arrayDataDepMap.levels(:,subsetToUse);


% define chunks into which data will be broken up for parallelization
nChunks = 10;
indx = round(linspace(0, numel(arrayData.tissues), nChunks+1))';
model_chunks = arrayfun(@(i) (indx(i)+1:indx(i+1))', (1:nChunks)', 'UniformOutput', false);
if strcmpi(chunk, 'test')
    model_indx = 1;
else
    model_indx = model_chunks{chunk};
end

%get the subset of the chunk
arrayData.tissues = arrayData.tissues(model_indx);
arrayData.levels = arrayData.levels(:,model_indx);


% initialize some variables
params = {};
init_models_depmap = {};
n_models = numel(model_indx);

% re-set to previous state if restart is requested
if ( restart )
    load(strcat('DepMap/init_models_depmap15-',num2str(chunk)));
    start_ind = numel(init_models_depmap) + 1;
else
    start_ind = 1;
end

for i = start_ind:n_models
    
    % Run with a really long time limit to ensure a somewhat high quality model.
	params.TimeLimit = 10000;
	params.MIPGap = 0.0004;
	init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);

    
    init_model.id = arrayData.tissues{i};
    init_models_depmap{i,1} = init_model;
    
    % save results
    if isnumeric(chunk)
        filename = strcat('DepMap/init_models_depmap15-',num2str(chunk));
    else
        filename = 'init_models_depmap15-TEST';
    end
    save(filename,'init_models_depmap');
    
end

