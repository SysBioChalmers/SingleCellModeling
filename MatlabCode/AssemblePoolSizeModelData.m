%Assembles the data for Fig 2.1 A and some associated supplementary plots
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

hcatModels = {};

for i = 1:10
    disp(num2str(i));
    fn = ['../data/PoolSize/hcat/chunk_' num2str(i) '.mat'];
    x = load(fn);
    hcatModels = [hcatModels;x.models];
end

load('../data/prepDataHumanGEM.mat');

%assemble into a matrix

baseModel = prepDataHumanGEM.refModel;

%filter the reactions by removing the ones that are always on
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
length(baseRxns)%7327 - looks reasonable

hcatRes = nan(length(baseRxns), length(hcatModels));

for i = 1:length(hcatModels)
    if rem(i, 10) == 0
        disp(num2str(i));
    end
    hcatRes(:,i) = ismember(baseRxns,hcatModels{i}.rxns);
end

jaccardResHcat = nan(length(hcatModels)/2,1);
for i = 1:length(jaccardResHcat)
    ind1 = (i-1)*2 + 1; %9 for i=5, ok
    ind2 = i*2; %10 for i=5, ok
    jaccardResHcat(i) = jaccard(hcatRes(:,ind1), hcatRes(:,ind2));
end
numRep = 20;
numPoints = length(jaccardResHcat)/numRep;%10

totJaccardHcat = nan(numPoints,1);
for i = 1:numPoints
    startInd = (i-1)*numRep+1;%41 for i=3, ok (20 per point)
    endInd = i*numRep;%60 for i=3, ok
    totJaccardHcat(i) = mean(jaccardResHcat(startInd:endInd));
end

clear hcatModels,hcatRes;

%%%%%%%%%%%%%%%%%%%%%
%Now t68k:
%%%%%%%%%%%%%%%%%%%%%

%not tested, same as for hcat

t68kModels = {};

for i = 1:10
    disp(num2str(i));
    fn = ['../data/PoolSize/t68k/chunk_' num2str(i) '.mat'];
    x = load(fn);
    t68kModels = [t68kModels;x.models];
end

t68kRes = nan(length(baseRxns), length(t68kModels));

for i = 1:length(t68kModels)
    if rem(i, 10) == 0
        disp(num2str(i));
    end
    t68kRes(:,i) = ismember(baseRxns,t68kModels{i}.rxns);
end

jaccardRest68k = nan(length(t68kModels)/2,1);
for i = 1:length(jaccardRest68k)
    ind1 = (i-1)*2 + 1;
    ind2 = i*2;
    jaccardRest68k(i) = jaccard(t68kRes(:,ind1), t68kRes(:,ind2));
end
numRep = 20;
numPoints = length(jaccardRest68k)/numRep;

totJaccardt68k = nan(numPoints,1);
for i = 1:numPoints
    startInd = (i-1)*numRep+1;
    endInd = i*numRep;
    totJaccardt68k(i) = mean(jaccardRest68k(startInd:endInd));
end

clear t68kModels,t68kRes;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bulk data
%%%%%%%%%%%%%%%%%%%%%%%%%%

x = load('../data/PoolSize/psBulkModels.mat');
bulkModels = x.models;
%compare all pairs with each other
numPairs = length(bulkModels) * (length(bulkModels)-1)/2; %28


bulkRes = nan(length(baseRxns), length(bulkModels));

for i = 1:length(bulkModels)
    bulkRes(:,i) = ismember(baseRxns,bulkModels{i}.rxns);
end

jaccardRestBulk = nan(numPairs,1);
currRes = 1;
for i = 1:(length(bulkModels)-1)
    for j = (i+1):length(bulkModels)
        jaccardRestBulk(currRes) = jaccard(bulkRes(:,i),bulkRes(:,j));
        currRes = currRes + 1;
    end
end
currRes %29, as expected - so, 28 pairs

%no need to test further, the code above is pretty straightforward


clear bulkModels,bulkRes;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pseudobulk data
%%%%%%%%%%%%%%%%%%%%%%%%%%

%not tested much, very similar to bulk above

x = load('../data/PoolSize/psPseudoBulkModels.mat');
bulkModels = x.models;
%compare all pairs with each other
numPairs = length(bulkModels) * (length(bulkModels)-1)/2;%28, ok


bulkRes = nan(length(baseRxns), length(bulkModels));

for i = 1:length(bulkModels)
    bulkRes(:,i) = ismember(baseRxns,bulkModels{i}.rxns);
end

jaccardRestPseudoBulk = nan(numPairs,1);
currRes = 1;
for i = 1:(length(bulkModels)-1)
    for j = (i+1):length(bulkModels)
        jaccardRestPseudoBulk(currRes) = jaccard(bulkRes(:,i),bulkRes(:,j));
        currRes = currRes + 1;
    end
end
currRes %29, as expected


clear bulkModels,bulkRes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pseudobulk TMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%not tested much, very similar to bulk above

x = load('../data/PoolSize/psPseudoBulkModelsTMM.mat');
bulkModels = x.models;
%compare all pairs with each other
numPairs = length(bulkModels) * (length(bulkModels)-1)/2;%28, ok


bulkRes = nan(length(baseRxns), length(bulkModels));

for i = 1:length(bulkModels)
    bulkRes(:,i) = ismember(baseRxns,bulkModels{i}.rxns);
end

jaccardRestPseudoBulkTMM = nan(numPairs,1);
currRes = 1;
for i = 1:(length(bulkModels)-1)
    for j = (i+1):length(bulkModels)
        jaccardRestPseudoBulkTMM(currRes) = jaccard(bulkRes(:,i),bulkRes(:,j));
        currRes = currRes + 1;
    end
end
currRes %29, as expected
mean(jaccardRestPseudoBulkTMM)
clear bulkModels,bulkRes;


d = struct();
d.hcat = totJaccardHcat;
d.t68k = totJaccardt68k;
d.bulk = mean(jaccardRestBulk);
d.bulkAll = jaccardRestBulk;
d.pseudoBulk = mean(jaccardRestPseudoBulk)
d.pseudoBulkAll = jaccardRestPseudoBulk
d.pseudoBulkTMM = mean(jaccardRestPseudoBulkTMM)
d.pseudoBulkTMMAll = jaccardRestPseudoBulkTMM
save('../data/PoolSizeModelRes.mat', 'd');

