%Evaluates the speed of the old tINIT and ftINIT by running them for 10 GTEx samples
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'
%use gtex
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('../data/tINIT_inputs_human.mat');

%Note - I removed the first two lines in the file to simplify reading of the file
gtexPath = '../data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm_mod.gct';

L = importdata(gtexPath,'\t');
s.data = L.data;
s.sampleIds = L.textdata(1, 3:end);
s.genes = L.textdata(2:end, 1);
s.data = s.data*10^6./sum(s.data,1);%TPM

%we also need to get rid of the version in the genes, i.e. ENSG00000243485.5, the .5 should be removed
y = strlength('ENSG00000243485');
s.genes = cellfun(@(x) x(1:y), s.genes,'UniformOutput',false);

load('../data/prepDataHumanGEMEns.mat');

model_indx = 1:10;

arrayData.genes = s.genes;
arrayData.tissues = s.sampleIds(model_indx);
arrayData.levels = s.data(:,model_indx);
arrayData.threshold = 1;

n_models = numel(arrayData.tissues);

params = {};
init_models_depmap = {};
n_models = numel(model_indx);

params.TimeLimit = 10000;
params.MIPGap = 0.0030;

old_tinit_exec_times = nan(n_models, 1);

for i = 1:n_models
    disp(['Running old tINIT model: ' num2str(i)])
    try
        tic
        init_model = getINITModel2(model,arrayData.tissues{i},[],[],arrayData,[],true,[],true,true,taskStruct,params);
        old_tinit_exec_times(i) = toc;
    catch
        disp('Failed') % the results will be NA, maybe just treat them as 10000 in that case.
    end
end
%save
save('../data/old_tinit_exec_times.mat', 'old_tinit_exec_times');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Now with the new algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%Note - I removed the first two lines in the file to simplify reading of the file
gtexPath = '../data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm_mod.gct';

L = importdata(gtexPath,'\t');
s.data = L.data;
s.sampleIds = L.textdata(1, 3:end);
s.genes = L.textdata(2:end, 1);
s.data = s.data*10^6./sum(s.data,1);%TPM

%we also need to get rid of the version in the genes, i.e. ENSG00000243485.5, the .5 should be removed
y = strlength('ENSG00000243485');
s.genes = cellfun(@(x) x(1:y), s.genes,'UniformOutput',false);


load('../data/prepDataHumanGEMEns.mat');
prepDataHumanGEMEns

model_indx = 1:10;

arrayDataNewAlg = struct();
arrayDataNewAlg.genes = s.genes;
arrayDataNewAlg.tissues = s.sampleIds(model_indx);
arrayDataNewAlg.levels = s.data(:,model_indx);
arrayDataNewAlg.threshold = 1;

n_models = numel(model_indx);

%paramsNewAlg = struct();
%paramsNewAlg.TimeLimit = 120;
%paramsNewAlg.MIPGap = 0.0004;

%milpSkipMets.simpleMets.mets = {'H2O';'Pi';'PPi';'H+';'O2';'CO2';'Na+'};
%milpSkipMets.simpleMets.compsToKeep = {'i'};

gtex_models_newalg_tmp = cell(n_models,1);

new_tinit_exec_times = nan(n_models, 1);

%we run it without the second step, i.e., all transport rxns etc. without GPR are just added.
for i = 1:n_models
    disp(['running model: ' num2str(i)])
    tic 
    mres = ftINIT(prepDataHumanGEMEns,arrayDataNewAlg.tissues{i},[],[],arrayDataNewAlg,{},getHumanGEMINITSteps('1+0'),false,true,[]);
    mres.id = arrayDataNewAlg.tissues{i};
    gtex_models_newalg_tmp{i,1} = mres;
    new_tinit_exec_times(i) = toc();
end

save('../data/new_tinit_exec_times.mat', 'new_tinit_exec_times');


