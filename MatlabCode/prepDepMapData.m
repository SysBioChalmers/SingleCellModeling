%Prepares the DepMap data, specifically by filtering out RNA-Seq samples for which no CRISPR data exists
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%% Load and prepare DepMap RNA-Seq data (cell lines)

% load RNA-Seq data from txt file
rna_data = readtable('../data/DepMap_tpm_ens.txt');

% load gene essentiality data (Achilles gene effect)
ach_data = readtable('F:/DepMap/Achilles_gene_effect.csv');
samples = ach_data.DepMap_ID;  % extract sample IDs

% filter RNA-Seq data to only include samples for which we have
% essentiality data
cellLineNames = rna_data.Properties.VariableNames;
%now replace '_' with '-'
cellLineNames = strrep(cellLineNames, '_', '-');

%{'Original column heading: 'ACH-001113''}

keep = ismember(cellLineNames, samples);

sum(keep) %891
sum(keep)/length(keep) % 65%, seems reasonable

% add RNA-Seq data to arrayData
arrayDataDepMap.genes = rna_data.gene;
arrayDataDepMap.tissues = cellLineNames(keep)';
arrayDataDepMap.levels = table2array(rna_data(:, keep));
arrayDataDepMap.threshold = 1;


% save tINIT inputs
save('../data/arrayDataDepMap.mat','arrayDataDepMap');











