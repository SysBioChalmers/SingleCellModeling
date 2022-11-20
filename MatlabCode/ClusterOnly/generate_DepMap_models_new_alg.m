function [] = generate_DepMap_models_new_alg(chunk, use1Plus1)
% generate tINIT models from RNA-Seq profiles
%
% Note: This function was designed specifically for use on the cluster
%
% Input:
%
%   chunk     integer ranging from 1 to 10, or 'test' to run a test that
%             generates only one model.
%
%   use1Plus1 if TRUE, run ftINIT with 1+1
%             (Default = FALSE)
%

if nargin < 2
    use1Plus1 = false;
end

% add paths
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');


%% run ftINIT

%we replace the arraydata with the one from depmap
load('DepMap/arrayDataDepMap.mat')
arrayData = arrayDataDepMap; %no need to add anything


% load inputs generated using "prepare_tINIT_inputs" function
load('prepDataHumanGEMEns.mat');

nChunks = 40;
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

depmap_models_newalg = cell(nModels, 1);

if use1Plus1
   INITSteps = '1+1';
else
   INITSteps = '1+0';
end   

for i = 1:nModels
     if (isempty(depmap_models_newalg{i}))
         disp(['running model: ' num2str(i)])
         tic 
         mres = ftINIT(prepDataHumanGEMEns,arrayData.tissues{i},[],[],arrayData,{},getHumanGEMINITSteps(INITSteps),false,true,[]);
         toc
         mres.id = arrayData.tissues{i};
         depmap_models_newalg{i,1} = mres;
     end
end

% save results
if isnumeric(chunk)
	if use1Plus1
		filename = strcat('DepMap/ftINIT2/depmap_models_newalg-',num2str(chunk));
	else
		filename = strcat('DepMap/ftINIT/depmap_models_newalg-',num2str(chunk));
	end
else
    filename = 'init_models_newalg-TEST';
end
save(filename,'depmap_models_newalg');

end
