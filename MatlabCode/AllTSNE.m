%This code generates all the data for fig 2 D-F and associated supplementary plots.
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

x = load('../data/LC3/pooledLC3Models.mat');
ctLCModels = x.models;
%length(ctLCModels)%16, ok

gtexIndModels = {};

for i = 1:10
    disp(num2str(i))
    fn = ['../data/GTExInd8/chunk_' num2str(i) '.mat'];
    x = load(fn);
    gtexIndModels = [gtexIndModels;x.models];
end

x = load('../data/L4/L4Models.mat');
L4Models = x.models;



load('../data/prepDataHumanGEM.mat');

baseModel = prepDataHumanGEM.refModel;

%now build a matrix
allMod = [gtexIndModels;L4Models;ctLCModels];

compMat = false(length(baseModel.rxns), length(gtexIndModels) + length(L4Models) + length(ctLCModels));



for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,allMod{i}.rxns);
end

rng(1);%random seed
proj_coords = tsne(double(compMat.'),'Distance','hamming','NumDimensions',2,'Exaggeration',6,'Perplexity',10);

%read the gtex tissue data
T = readtable('../data/GTExInd8/gtexIndSampTissues8.txt', 'ReadVariableNames',false);
gtexTissues = table2cell(T(:,1));

%export to R
d = struct();
d.tsneX = proj_coords(:,1);
d.tsneY = proj_coords(:,2);
d.tissues = [gtexTissues;'L4 Lung';'L4 Spleen';repmat({'scLC'}, length(ctLCModels),1)];
save('../data/allTSNE.mat', 'd');

clear gtexIndModels;

%first filter the compMat by removing the "always-on" rxns:
x = prepDataHumanGEM;
alwaysOn = x.toIgnoreExch | ...
           x.toIgnoreImportRxns | ...
           x.toIgnoreSimpleTransp | ...
           x.toIgnoreAdvTransp | ...
           x.toIgnoreSpont | ...
           x.toIgnoreS | ...
           x.toIgnoreCustomRxns;
baseRxns = baseModel.rxns(~alwaysOn);
baseRxns = baseRxns(~ismember(baseRxns,x.essentialRxns));
%length(baseRxns)%7327

%now filter the matrix
filt = ismember(baseModel.rxns, baseRxns);



compMat2 = compMat(filt,:);
%size(compMat2)%looks ok, 7327         283

%calc the average jaccard wihtin groups of gtex
%numTiss = 265/5;

numTiss = 53;
jaccWithinGtexGroups = nan(numTiss,1);
sampPerTiss = 8;
for i = 1:numTiss
    startInd = (i-1)*sampPerTiss + 1;
    endInd = i*sampPerTiss;
    jaccs = nan(sampPerTiss*(sampPerTiss-1)/2,1);
    currJac = 1;
    for j = 0:(sampPerTiss-2)
       for k = (j+1):(sampPerTiss-1)
           jaccs(currJac) = jaccard(compMat2(:,startInd + j), compMat2(:,startInd + k));
           currJac = currJac + 1;
       end
    end
    jaccWithinGtexGroups(i) = mean(jaccs);
end

jaccAcrossGtexGroups = nan((numTiss*sampPerTiss) *(numTiss*sampPerTiss-1)/2,1);
currJac = 1;
for i = 1:(numTiss*sampPerTiss-1)
    jaccs = nan((sampPerTiss-1)*sampPerTiss/2,1);
    for j = (i+1):(numTiss*sampPerTiss)
       jaccAcrossGtexGroups(currJac) = jaccard(compMat2(:,i), compMat2(:,j));
       currJac = currJac + 1;
    end
end

acrossTPM = jaccAcrossGtexGroups;
withinTPM = jaccWithinGtexGroups;

%also within/across lung

lungImmune = [1;3;5;6;7;8;9;10;11;12;13;14;16];
%length(lungImmune)
%lungEpithelial = c(2,4,15)


jaccWithinScs = nan(12*13/2,1);
currJac = 1;
indices = lungImmune + numTiss*sampPerTiss + 2;%2 for L4
for i = 1:(length(indices)-1)
    for j = (i+1):length(indices)
       jaccWithinScs(currJac) = jaccard(compMat2(:,indices(i)), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end
jaccWithinScTPM = mean(jaccWithinScs);

%compare sc and gtex
bloodInd = find(strcmp(gtexTissues, 'Whole Blood'));%1-8
jaccAcrossLungs = nan(sampPerTiss*length(indices),1);
currJac = 1;
for i = bloodInd(1):bloodInd(sampPerTiss)
    for j = 1:length(indices)
       jaccAcrossLungs(currJac) = jaccard(compMat2(:,i), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end
jaccAcrossLungTPM = mean(jaccAcrossLungs);

%compare sc and any gtex
jaccAcrossLungAnys = nan(numTiss*sampPerTiss*length(indices),1);
currJac = 1;
for i = 1:(numTiss*sampPerTiss)
    for j = 1:length(indices)
       jaccAcrossLungAnys(currJac) = jaccard(compMat2(:,i), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end
jaccAcrossLungAndAnyTPM = mean(jaccAcrossLungAnys);

d = struct();
d.acrossGTEx = acrossTPM;
d.withinGTEx = withinTPM;
d.withinScLung = jaccWithinScs;
d.scVsGTExLung = jaccAcrossLungs;
d.scVsGTExAny = jaccAcrossLungAnys;
d.compMat = compMat2;
save('../data/structCompTPM.mat', 'd');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now the TMM normalized samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmmModels = {};

for i = 1:10
    disp(num2str(i))
    fn = ['../data/TMMData/chunk_' num2str(i) '.mat'];
    x = load(fn);
    tmmModels = [tmmModels;x.models];
end

load('../data/prepDataHumanGEM.mat');

baseModel = prepDataHumanGEM.refModel;

%now build a matrix
compMat = false(length(baseModel.rxns), length(tmmModels));


for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,tmmModels{i}.rxns);
end

rng(1);%random seed
proj_coords_tmm = tsne(double(compMat.'),'Distance','hamming','NumDimensions',2,'Exaggeration',6,'Perplexity',10);

%export to R
d = struct();
d.tsneX = proj_coords_tmm(:,1);
d.tsneY = proj_coords_tmm(:,2);
save('../data/tmmTSNE.mat', 'd');

clear tmmModels;


%first filter the compMat by removing the "always-on" rxns:
x = prepDataHumanGEM;
alwaysOn = x.toIgnoreExch | ...
           x.toIgnoreImportRxns | ...
           x.toIgnoreSimpleTransp | ...
           x.toIgnoreAdvTransp | ...
           x.toIgnoreSpont | ...
           x.toIgnoreS | ...
           x.toIgnoreCustomRxns;
baseRxns = baseModel.rxns(~alwaysOn);
baseRxns = baseRxns(~ismember(baseRxns,x.essentialRxns));
%length(baseRxns)%7327

%now filter the matrix
filt = ismember(baseModel.rxns, baseRxns);



compMat2 = compMat(filt,:);
%size(compMat2)%looks ok 7340         442


%calc the average jaccard wihtin groups of gtex
numTiss = 53;
jaccWithinGtexGroups = nan(numTiss,1);
sampPerTiss = 8;
for i = 1:numTiss
    startInd = (i-1)*sampPerTiss + 1;
    endInd = i*sampPerTiss;
    jaccs = nan(sampPerTiss*(sampPerTiss-1)/2,1);
    currJac = 1;
    for j = 0:(sampPerTiss-2)
       for k = (j+1):(sampPerTiss-1)
           jaccs(currJac) = jaccard(compMat2(:,startInd + j), compMat2(:,startInd + k));
           currJac = currJac + 1;
       end
    end
    jaccWithinGtexGroups(i) = mean(jaccs);
end

jaccAcrossGtexGroups = nan((numTiss*sampPerTiss) *(numTiss*sampPerTiss-1)/2,1);
currJac = 1;
for i = 1:(numTiss*sampPerTiss-1)
    jaccs = nan((sampPerTiss-1)*sampPerTiss/2,1);
    for j = (i+1):(numTiss*sampPerTiss)
       jaccAcrossGtexGroups(currJac) = jaccard(compMat2(:,i), compMat2(:,j));
       currJac = currJac + 1;
    end
end

acrossTMM = jaccAcrossGtexGroups;
withinTMM = jaccWithinGtexGroups;

%also within/across lung

lungImmune = [1;3;5;6;7;8;9;10;11;12;13;14;16];
%length(lungImmune)
%lungEpithelial = c(2,4,15)


jaccWithinScs = nan(12*13/2,1);
currJac = 1;
indices = lungImmune + numTiss*sampPerTiss + 2;%2 for L4
for i = 1:(length(indices)-1)
    for j = (i+1):length(indices)
       jaccWithinScs(currJac) = jaccard(compMat2(:,indices(i)), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end

%compare sc and gtex
bloodInd = find(strcmp(gtexTissues, 'Whole Blood'));%1-8
jaccAcrossLungs = nan(sampPerTiss*length(indices),1);
currJac = 1;
for i = bloodInd(1):bloodInd(sampPerTiss)
    for j = 1:length(indices)
       jaccAcrossLungs(currJac) = jaccard(compMat2(:,i), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end

%compare sc and any gtex
jaccAcrossLungAnys = nan(numTiss*sampPerTiss*length(indices),1);
currJac = 1;
for i = 1:(numTiss*sampPerTiss)
    for j = 1:length(indices)
       jaccAcrossLungAnys(currJac) = jaccard(compMat2(:,i), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end

d = struct();
d.acrossGTEx = acrossTMM;
d.withinGTEx = withinTMM;
d.withinScLung = jaccWithinScs;
d.scVsGTExLung = jaccAcrossLungs;
d.scVsGTExAny = jaccAcrossLungAnys;
d.compMat = compMat2;
save('../data/structCompTMM.mat', 'd');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%And the quantile normalized samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qModels = {};

for i = 1:10
    disp(num2str(i))
    fn = ['../data/QuantileData/chunk_' num2str(i) '.mat'];
    x = load(fn);
    qModels = [qModels;x.models];
end

load('../data/prepDataHumanGEM.mat');

baseModel = prepDataHumanGEM.refModel;

%now build a matrix
compMat = false(length(baseModel.rxns), length(qModels));

for i = 1:size(compMat,2)
    compMat(:,i) = ismember(baseModel.rxns,qModels{i}.rxns);
end

rng(1);%random seed
proj_coords_q = tsne(double(compMat.'),'Distance','hamming','NumDimensions',2,'Exaggeration',6,'Perplexity',10);

%export to R
d = struct();
d.tsneX = proj_coords_q(:,1);
d.tsneY = proj_coords_q(:,2);
save('../data/quantileTSNE.mat', 'd');

clear qModels;

%%%%%%%%%%%%%%%
%calculate mean jaccard within gtex organs

%first filter the compMat by removing the "always-on" rxns:
x = prepDataHumanGEM;
alwaysOn = x.toIgnoreExch | ...
           x.toIgnoreImportRxns | ...
           x.toIgnoreSimpleTransp | ...
           x.toIgnoreAdvTransp | ...
           x.toIgnoreSpont | ...
           x.toIgnoreS | ...
           x.toIgnoreCustomRxns;
baseRxns = baseModel.rxns(~alwaysOn);
baseRxns = baseRxns(~ismember(baseRxns,x.essentialRxns));
%length(baseRxns)%7327

%now filter the matrix
filt = ismember(baseModel.rxns, baseRxns);



compMat2 = compMat(filt,:);
%size(compMat2)%looks ok 7340         442


%calc the average jaccard wihtin groups of gtex
numTiss = 53;
jaccWithinGtexGroups = nan(numTiss,1);
sampPerTiss = 8;
for i = 1:numTiss
    startInd = (i-1)*sampPerTiss + 1;
    endInd = i*sampPerTiss;
    jaccs = nan(sampPerTiss*(sampPerTiss-1)/2,1);
    currJac = 1;
    for j = 0:(sampPerTiss-2)
       for k = (j+1):(sampPerTiss-1)
           jaccs(currJac) = jaccard(compMat2(:,startInd + j), compMat2(:,startInd + k));
           currJac = currJac + 1;
       end
    end
    jaccWithinGtexGroups(i) = mean(jaccs);
end

jaccAcrossGtexGroups = nan((numTiss*sampPerTiss) *(numTiss*sampPerTiss-1)/2,1);
currJac = 1;
for i = 1:(numTiss*sampPerTiss-1)
    jaccs = nan((sampPerTiss-1)*sampPerTiss/2,1);
    for j = (i+1):(numTiss*sampPerTiss)
       jaccAcrossGtexGroups(currJac) = jaccard(compMat2(:,i), compMat2(:,j));
       currJac = currJac + 1;
    end
end

acrossQ = jaccAcrossGtexGroups;
withinQ = jaccWithinGtexGroups;

%also within/across lung

lungImmune = [1;3;5;6;7;8;9;10;11;12;13;14;16];
%length(lungImmune)
%lungEpithelial = c(2,4,15)


jaccWithinScs = nan(12*13/2,1);
currJac = 1;
indices = lungImmune + numTiss*sampPerTiss + 2;%2 for L4
for i = 1:(length(indices)-1)
    for j = (i+1):length(indices)
       jaccWithinScs(currJac) = jaccard(compMat2(:,indices(i)), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end

%compare sc and gtex
bloodInd = find(strcmp(gtexTissues, 'Whole Blood'));%1-8
jaccAcrossLungs = nan(sampPerTiss*length(indices),1);
currJac = 1;
for i = bloodInd(1):bloodInd(sampPerTiss)
    for j = 1:length(indices)
       jaccAcrossLungs(currJac) = jaccard(compMat2(:,i), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end

%compare sc and any gtex
jaccAcrossLungAnys = nan(numTiss*sampPerTiss*length(indices),1);
currJac = 1;
for i = 1:(numTiss*sampPerTiss)
    for j = 1:length(indices)
       jaccAcrossLungAnys(currJac) = jaccard(compMat2(:,i), compMat2(:,indices(j)));
       currJac = currJac + 1;
    end
end

d = struct();
d.acrossGTEx = acrossQ;
d.withinGTEx = withinQ;
d.withinScLung = jaccWithinScs;
d.scVsGTExLung = jaccAcrossLungs;
d.scVsGTExAny = jaccAcrossLungAnys;
d.compMat = compMat2;
save('../data/structCompQuantile.mat', 'd');
