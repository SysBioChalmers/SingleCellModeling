cd C:\Work\MatlabCode\projects\SingleCellModeling\SingleCellModeling\data
%first read the bootstrap sample names
data = fileread('MetAtlas/ClusterIds.txt');
sampleNames = strsplit(data,'\n');

load('prepDataHumanGEMEns.mat');

prepDataHumanGEMEns
totMat = nan(length(prepDataHumanGEMEns.refModel.rxns), length(sampleNames));

for i = 1:length(sampleNames)
    disp(i)
    inFilename = ['MetAtlas/BootstrapModels/Bootstrap_' sampleNames{i} '_models.mat'];
    x = load(inFilename);
    mat = x.rxnOnData;
    totMat(:,i) = sum(mat,2)/100;
end

%count number of tissues
nms = split(sampleNames, '_');
length(unique(nms(:,:,1))) %19

%Write text file
%first get the cell type names
T = readtable('MetAtlas/tissue_ct.txt');
[numTot,~] = size(T);%444
tableSampleNames = cell(numTot,1);
exportSampleNames = cell(numTot,1);
for i = 1:numTot
    tableSampleNames{i} = [T.tissueMD{i} '_' num2str(T.cluster(i))];
    exportSampleNames{i} = [T.tissueCT{i} '_' T.celltype{i} '_' num2str(T.cluster(i))];
end

[~, ia, ib] = intersect(sampleNames, tableSampleNames);
expNames = sampleNames.';%allocate right size
expNames(ia) = exportSampleNames(ib);
celltypes = sampleNames.';%allocate right size
celltypes(ia) = T.celltype(ib);
tissues = sampleNames.';%allocate right size
tissues(ia) = T.tissueCT(ib);

%test
%table(sampleNames.', expNames) %looks good


%now create a table with the names as columns, the reactions as the first row, and the matrix following
%as columns per tissue+ct.
doubles = repmat("double",length(sampleNames),1);
varTypes = ["string";doubles].';

varNames = ["id";string(expNames)].';
numRxns = length(prepDataHumanGEMEns.refModel.rxns);
expT = table('Size',[numRxns, length(sampleNames)+1],'VariableTypes',varTypes,'VariableNames',varNames);
expT(:,1) = prepDataHumanGEMEns.refModel.rxns;
expT(:,2:(length(sampleNames)+1)) = num2cell(totMat);

%Test that it looks ok:
%all(all(table2array(expT(1000:1010,100:110)) == totMat(1000:1010,99:109))) %OK
writetable(expT, 'MetAtlas/HPA_single-cell_reactions.tsv','FileType','text','Delimiter','\t');

%also run t-sne on the reactions included in the first step of ftINIT
mask = prepDataHumanGEMEns.toIgnoreExch | ...
       prepDataHumanGEMEns.toIgnoreImportRxns | ...
       prepDataHumanGEMEns.toIgnoreSimpleTransp | ...
       prepDataHumanGEMEns.toIgnoreAdvTransp | ...
       prepDataHumanGEMEns.toIgnoreSpont | ...
       prepDataHumanGEMEns.toIgnoreS | ...
       prepDataHumanGEMEns.toIgnoreCustomRxns;

matForTsne = totMat(~mask,:);
size(matForTsne )%7527         202
rng(1);%random seed
proj_coords = tsne(matForTsne.','NumDimensions',2,'Exaggeration',6,'Perplexity',10);
d = struct();
d.tsneX = proj_coords(:,1);
d.tsneY = proj_coords(:,2);
d.celltype = celltypes;
d.tissue = tissues;
save('../data/MetAtlas/metAtlasTsne.mat', 'd');

   