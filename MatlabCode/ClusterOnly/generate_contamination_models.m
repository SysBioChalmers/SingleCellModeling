%Generates models for GTEx in TPM, the new tINIT version
function generate_contamination_models(chunk)

%for debugging locally
%cd C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/data

% add paths (comment out if you run locally)
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');



load('contamination/PSPoolSizesCont.mat');

%Calculate cont_index and ps_index
cont_index = floor(1 + (chunk-1)/8);
ps_index = chunk - (cont_index-1)*8;

nModels = 120;


indOffset = (ps_index - 1)*nModels;


%% run tINIT

% load inputs generated using "prepare_tINIT_inputs" function
load('prepDataHumanGEM.mat');

arrayData = arrayDataCont{cont_index}; 

rxnContent = false(length(prepDataHumanGEM.refModel.rxns), nModels);

for i = 1:nModels
     disp(['running model: ' num2str(i)])
     mres = ftINIT(prepDataHumanGEM,arrayData.tissues{i + indOffset},[],[],arrayData,{},getHumanGEMINITSteps('1+0'),false,true,[]);
     rxnContent(:,i) = ismember(prepDataHumanGEM.refModel.rxns, mres.rxns);
end

% save results
fn = ['contamination/cont_models_' num2str(cont_index) '_' num2str(ps_index) '.mat'];
save(fn,'rxnContent');

