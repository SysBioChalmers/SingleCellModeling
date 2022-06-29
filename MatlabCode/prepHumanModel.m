%Runs prepInitModel for Human-GEM, including filtering of reactions.
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%load Human-GEM and tasks
%load('C:/Work/MatlabCode/components/human-GEM/Human-GEM_1_12/Human-GEM/model/Human-GEM.mat')
load('C:/Work/MatlabCode/components/human-GEM/Human-GEMftINIT/Human-GEM/model/Human-GEM.mat')

%first filter some duplicate reactions in the model to speed up the calculations and aviod randomness.
%let's not do this, it is difficult to track which were actually removed in the latest version of Human-GEM
m = removeUnwantedReactions(ihuman);

prepDataHumanGEM = prepHumanModelForftINIT(m, true);
save('../data/prepDataHumanGEM.mat','prepDataHumanGEM');

prepDataHumanGEMEns = prepHumanModelForftINIT(m, false);
save('../data/prepDataHumanGEMEns.mat','prepDataHumanGEMEns');


%x = prepDataHumanGEM.toIgnoreExch | prepDataHumanGEM.toIgnoreImportRxns | prepDataHumanGEM.toIgnoreSimpleTransp | prepDataHumanGEM.toIgnoreAdvTransp | prepDataHumanGEM.toIgnoreSpont | prepDataHumanGEM.toIgnoreS;
%sum(x)
%sum(prepDataHumanGEM.toIgnoreAllWithoutGPRs & ~x)%902 rxns extra if we skip all rxns without GPRs

%constructEquations(prepDataHumanGEM.refModel, prepDataHumanGEM.refModel.rxns(prepDataHumanGEM.toIgnoreAllWithoutGPRs & ~x))

%TEMP: Check if the unwanted reactions are removed in v 1.12

