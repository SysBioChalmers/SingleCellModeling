%Generates the pool size data for Fig. 2A-C.
%Note that this code uses SingleCellToolbox (https://github.com/SysBioChalmers/SingleCellToolbox). 
%See https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0239495 for details
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'
%load('tINIT_preprocess_inputs')
%prep for scoring
%essentialReactions = any(preData.essentialRxnMat,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate pool size data sets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%assemble single-cell datasets
hca = DsHcaCB.get();
hcat = hca.cellSubset(hca.cellType == Celltype.TCellCD4Pos | hca.cellType == Celltype.TCellCD8Pos | hca.cellType == Celltype.TCellReg |  hca.cellType == Celltype.TCell );
hcab = hca.cellSubset(hca.cellType == Celltype.BCell);
[lcc,lch] = DsLC.get();
lcct = lcc.cellSubset(lcc.cellType == Celltype.TCellCD4Pos | lcc.cellType == Celltype.TCellCD8Pos | lcc.cellType == Celltype.TCellReg);
lccm = lcc.cellSubset(lcc.cellType == Celltype.Macrophage);

t68k = DsPbmc68k.get().cellSubset(DsPbmc68k.get().cellType == Celltype.TCellCD4Pos | DsPbmc68k.get().cellType == Celltype.TCellCD8Pos | DsPbmc68k.get().cellType == Celltype.TCellReg);
[~,~,cd8] = DsGSE112845.get();

% For checking contamination from other cell types
lcc4c = lcc.cellSubset((lcc.cellType == Celltype.Malignant) & strcmp(lcc.sampleIds,'4'));
lcc4t = lcc.cellSubset(lcc.cellType == Celltype.TCellCD4Pos | lcc.cellType == Celltype.TCellCD8Pos | lcc.cellType == Celltype.TCellReg & strcmp(lcc.sampleIds,'4'));


%keep this separate
melt = DsMel.get().cellSubset(DsMel.get().cellType == Celltype.TCell);

% Keep the size reasonable, otherwise we may run out of memory
hcatlim = hcat.cellSubset(1:25000);
hcablim = hcab.cellSubset(1:25000);

dss = { hcatlim, ... 
        hcablim, ...
        t68k, ...
        lcct, ...
        lccm, ...
        cd8 ...
      };

%export the data to be able to import it in R
Write10xMatrix('../data/exportToR/datasets/hcat', hcatlim);
Write10xMatrix('../data/exportToR/datasets/hcab', hcablim);
Write10xMatrix('../data/exportToR/datasets/t68k', t68k);
Write10xMatrix('../data/exportToR/datasets/lcct', lcct);
Write10xMatrix('../data/exportToR/datasets/lccm', lccm);
Write10xMatrix('../data/exportToR/datasets/tcd8', cd8);

%the smart-seq 2 data is different
if ~exist('../data/exportToR/datasets/melt', 'dir')
   mkdir('../data/exportToR/datasets/melt')
end
writematrix(melt.data, '../data/exportToR/datasets/melt/matrix.txt');
%write genes
fileID = fopen('../data/exportToR/datasets/melt/genes.txt', 'w');
for i = 1:length(melt.genes)
    fprintf(fileID, '%s\n', melt.genes{i});
end
fclose(fileID);



numPoints = 30;
numRep = 60;

%quick test
%dss = {hcatlim};
%numPoints = 5;
%numRep = 5;

%set random seed
rng(1);
  
psSamples = cell(1,numel(dss));
poolSizes = zeros(numel(dss),numPoints);
progbar = ProgrBar('Generating pool size samples');
for i = 1:numel(dss)
   [psSamples{i},poolSizes(i,:)] = genPoolSizeSamples(dss{i}, 10, 10000, numPoints, numRep, progbar.GetSubContext(1/numel(dss)));
end

progbar.Done();

%generate sm2 data, i.e. MelT
%set random seed
rng(1);
[psSamplesMelT,poolSizesMelT] = genPoolSizeSamples(melt, 10, 1000, numPoints, numRep);

%Generate pool size samples with cell contamination from other cell types
%set random seed
rng(1);
[psSamplesCont,poolSizesCont] = genPoolSizeSamplesWithContamination(lcc4t, lcc4c, 500, 5000, 8, numRep, [0 .02 .05 .1 .2]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate data for model generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%use two datasets only: hcatlim and t68k
numPoints = 10;
numRep = 20;
[psSamplesHcat,poolSizesHcat] = genPoolSizeSamples(hcatlim, 150, 10000, numPoints, numRep);
psSamplesHcat.writeToTextFile('../data/PoolSize/hcatPSSamples.txt');

[psSamplesT68k,poolSizesT68k] = genPoolSizeSamples(t68k, 150, 10000, numPoints, numRep);
psSamplesT68k.writeToTextFile('../data/PoolSize/t68kPSSamples.txt');
d = struct();
d.poolSizesHcat = poolSizesHcat;
d.poolSizesT68k = poolSizesT68k;
save('../data/PSPoolSizes.mat','d')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Measure variation across samples in scRNA-Seq, HCA_CB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unique(hcat.sampleIds) %8 samples, CB1 - CB8
s = Samples;
s.genes = hcat.genes;
s.sampleIds = strsplit(num2str(1:8));
s.data = nan(length(s.genes),8);
for i = 1:8
    subdat = hcat.data(:,strcmp(hcat.sampleIds,['CB' num2str(i)]));
    s.data(:,i) = sum(subdat,2);
    %TPM
    s.data(:,i) = s.data(:,i)*10^6/sum(s.data(:,i));
end

%test
s %looks reasonable


%also check number of cells per sample
cellsPerSample = nan(8,1);
for i = 1:8
    cellsPerSample(i) = sum(strcmp(hcat.sampleIds,['CB' num2str(i)]));
end
cellsPerSample %at least 10000 per patient, ok

%cellsPerSample =
%
%       10154
%       24384
%       30574
%       15605
%       21959
%       22371
%       24549
%       18016

%save for running tINIT
save('../data/PoolSize/pseudoBulkModelData.mat','s');
s.writeToTextFile('../data/PoolSize/pseudoBulkModelDataCPM.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate reaction scores (takes a few hours to run...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('C:/Work/MatlabCode/components/human-GEM/Human-GEM_1_12/Human-GEM/model/Human-GEM.mat')
[~,deletedDeadEndRxns] = simplifyModel(ihuman,true,false,true,true,true);
cModel = removeReactions(ihuman,deletedDeadEndRxns,true,true);
%convert genes
[cModel.grRules, cModel.genes, cModel.rxnGeneMat] = translateGrRules(cModel.grRules, 'Name');


%for the contaminated samples only
rng(1);  
progbar = ProgrBar('Generating contaminated reaction scores - (8h or so)');
rxnScoresCont = cell(1,size(psSamplesCont,2));
for i = 1:numel(psSamplesCont)
   tmpRxnScores = getRxnScoresFromSamples(psSamplesCont{i}, cModel, 1, progbar.GetSubContext(1/numel(psSamplesCont)));
   %change to 3D format
   rxnScoresCont{i} = rxnScoresTo3D(tmpRxnScores, poolSizesCont);
   disp(i)
end

progbar.Done();

save('../data/PSrxnScoresContaminated.mat', 'rxnScoresCont', 'poolSizesCont');
load('../data/PSrxnScoresContaminated.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calc Jaccard scores for contaminated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numContPoints = size(rxnScoresCont{1},2);
X = zeros(numContPoints,numel(rxnScoresCont));
Y = zeros(numContPoints,numel(rxnScoresCont));

progbar = ProgrBar('Generating data for pool size vs reaction scores Jaccard, contaminated');
fullVals = cell(numel(rxnScoresCont),1);
for i = 1:numel(rxnScoresCont)
   [vals, fullVals{i}] = calcJaccard(rxnScoresCont{i}, progbar.GetSubContext(1/numel(rxnScoresCont)));
   X(:,i) = poolSizesCont.';
   Y(:,i) = vals.';
end

resPoolSizeVsReactionScoresJaccardCont.X = X;
resPoolSizeVsReactionScoresJaccardCont.Y = Y;
save('../data/PSresPoolSizeVsReactionScoresJaccardCont.mat', 'resPoolSizeVsReactionScoresJaccardCont');

progbar.Done();

%save data for statistical test
save('../data/PSContStatTestData.mat', 'fullVals');


