%Generates all the bootstrap models for one celltype, specified by chunk, for the MCOR3 dataset. 
%Designed to be called from a cluster script.
function generate_MCOR3_bootstrap_models_run5(chunk)


sampleNames = {'L2_3 IT';...
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
    
sampleName = sampleNames{chunk};
infile = ['MCOR3Data5/Bootstrap_' sampleName '.txt'];
outfile = ['MCOR3Data5/Bootstrap_' sampleName '_models.mat'];


%for debugging
%cd 'C:/Work/MatlabCode/projects/HMASandbox/HMA_Sandbox/Single-cell Modeling/tINIT/clusterCode/code'%remove later

% add paths
addpath(genpath('../../components/RAVEN'));
addpath(genpath('../../components/COBRA'));
addpath(genpath('../../components/Human-GEM'));

%setRavenSolver('gurobi');

T = readtable(infile);


genes = T{:,1};
dataMat = table2array(T(:,2:end));

nModels = size(dataMat,2);

%% run tINIT

load('prepDataMouseGEM.mat');

paramsNewAlg = struct();
paramsNewAlg.TimeLimit = 120;
paramsNewAlg.MIPGap = 0.0004;

milpSkipMets.simpleMets.mets = {'H2O';'Pi';'PPi';'H+';'O2';'CO2';'Na+'};
milpSkipMets.simpleMets.compsToKeep = {'i'};

arrayData = struct(); 
arrayData.genes = genes;
arrayData.levels = dataMat;
arrayData.tissues = arrayfun(@(x) num2str(x), 1:nModels, 'UniformOutput', false);
arrayData.threshold = 1;
%CPM the data
arrayData.levels = arrayData.levels.*10^6./sum(arrayData.levels,1);
%sum(arrayData.levels,1)

models = cell(nModels,1);

for i = 1:nModels
     disp(['running model: ' num2str(i)])
     tic %we run this without step 3
     mres = getINITModel9(prepDataMouseGEM,arrayData.tissues{i},[],[],arrayData,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,milpSkipMets,true,false,paramsNewAlg);
     toc
     mres.id = arrayData.tissues{i};
     models{i,1} = mres;
end

% save results
save(outfile,'models');

