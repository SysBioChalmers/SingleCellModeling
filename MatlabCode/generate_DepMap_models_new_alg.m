%Generates the models for the DepMap data with ftINIT
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

load('../data/arrayDataDepMap.mat')


%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('../data/prepDataHumanGEMEns.mat');
prepDataHumanGEMEns

model_indx = 1:15; %length(s.sampleIds); %run with 15 samples

n_models = numel(model_indx);

depmap_models_newalg = cell(n_models,1);

arrayData = arrayDataDepMap; %no need to add anything

for i = 1:n_models
    % First try to run tINIT with shorter time limit. If it fails, then
    % try again with a longer time limit (and then again).
     if (isempty(depmap_models_newalg{i}))
         disp(['running model: ' num2str(i)])
         tic %without step 3 for now
         mres = ftINIT(prepDataHumanGEMEns,arrayData.tissues{i},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);
         toc
         mres.id = arrayData.tissues{i};
         depmap_models_newalg{i,1} = mres;
     end
end


% save results
save('../data/depmap_models_newalg.mat','depmap_models_newalg');
