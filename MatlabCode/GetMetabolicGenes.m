%Extracts the genes used in the model, which defines the metabolic genes
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'
load('C:\Work\MatlabCode\components\MouseGEM\Mouse-GEM\model\Mouse-GEM.mat')%version 1.2.0


%remove some kinases that only do signaling, i.e. add a phosphate group to a protein.
mouseGEMFilt = removeReactions(mouseGEM, {'MAR09577', 'MAR09578', 'MAR09579', 'MAR07617', 'MAR07618'});


[genes,rxnMapping] = getGenesFromGrRules(mouseGEMFilt.grRules);%takes a minute to run

%genes 
length(genes)%2961 - before filtering: 3513 - so those signaling reactions introduces ~550 genes!
metabolicGenesMouse = genes;
t = table(metabolicGenesMouse);
writetable(t, '../data/metabolicGenesMouse.txt', 'WriteVariableNames', false);


%Repeat for human as well
load('C:\Work\MatlabCode\components\human-GEM\Human-GEM\model\Human-GEM.mat')%version 1.10.0


%remove some kinases that only do signaling, i.e. add a phosphate group to a protein.
humanGEMFilt = removeReactions(ihuman, {'MAR09577', 'MAR09578', 'MAR09579', 'MAR07617', 'MAR07618'});

humanGEMFilt.grRules = translateGrRules(humanGEMFilt.grRules,'Name');
[metabolicGenesHuman,~] = getGenesFromGrRules(humanGEMFilt.grRules);
%genes
length(metabolicGenesHuman)%3060
t = table(metabolicGenesHuman);
writetable(t, 'data/metabolicGenesHuman.txt', 'WriteVariableNames', false);


