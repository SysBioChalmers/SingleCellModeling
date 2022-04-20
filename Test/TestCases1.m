
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TC 0101 - genPoolSizeSamples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ds = SCDataset;
ds.data = eye(50);
ds.genes = arrayfun(@num2str, 1:50, 'UniformOutput', false).';
ds.cellIds = arrayfun(@num2str, 1:50, 'UniformOutput', false);
ps = genPoolSizeSamples(ds, 5, 25, 2, 2);

%check that all are TPM/CPM
all(sum(ps.data,1) == 10^6) %ok

%check that there are 25 cells in the 8th pool, and that there are no reuse of cells
setdiff(unique(ps.data(:,8)), [0;10^6/25]) %empty, as it should be

%check that there are no overlaps between the pairs (only check col 7 and 8)
all(abs(ps.data(:,7) - ps.data(:,8)) == ps.data(:,7) + ps.data(:,8)) %ok

%all in all ok.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TC 0102 - genPoolSizeSamplesWithContamination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ds1 = SCDataset;
ds1.data = [eye(50);zeros(50)];
ds1.genes = arrayfun(@num2str, 1:100, 'UniformOutput', false).';
ds1.cellIds = arrayfun(@num2str, 1:50, 'UniformOutput', false);

ds2 = SCDataset;
ds2.data = [zeros(50);eye(50)];
ds2.genes = arrayfun(@num2str, 1:100, 'UniformOutput', false).';
ds2.cellIds = arrayfun(@num2str, 1:50, 'UniformOutput', false);

psTmp = genPoolSizeSamplesWithContamination(ds1, ds2, 5, 25, 2, 2, 0.20);
ps = psTmp{1};


sum1 = sum(ps.data(1:50,:),1);
sum2 = sum(ps.data(51:100,:),1);

%check that the first in each pair has all counts in the upper half
all(sum1(:,[1,3,5,7]) == 10^6) %ok

%check that 20% of the reads belong to the lower 50 genes for the second in each pair
all(sum1(:,[2,4,6,8]) == 800000) %ok
all(sum2(:,[2,4,6,8]) == 200000) %ok

%all in all ok.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define a model (copied from tINIT tests)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testModel = struct();
testModel.id = 'testModel';
testModel.rxns = {};
testModel.S=[];
testModel.rev=[];
testModel.mets = {'a';'a2';'b';'c';'d';'e';'e2';'f'};
testModel.metNames = {'a';'a';'b';'c';'d';'e';'e';'f'};
testModel.comps = {'s';'c'};
testModel.compNames = testModel.comps;
testModel.metComps = [1;2;2;2;2;2;1;2];
testModel.genes = {'G1';'G2';'G3';'G4';'G5';'G6';'G7';'G8';'G9';'G10'};
testModel.grRules = {};
testModel.rxnGeneMat = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10'};
rxnsToAdd.equations = {'=> a[s]';...
                       'a[s] <=> a[c]';...
                       'a[c] <=> b[c] + c[c]';...
                       'a[c] <=> 2 d[c]';...
                       'b[c] + c[c] => e[c]';...
                       '2 d[c] => e[c]';...
                       'e[c] => e[s]';...
                       'e[s] =>';...
                       'a[c] <=> f[c]';...
                       'f[c] <=> e[c]'};
rxnsToAdd.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};
testModel = addRxns(testModel,rxnsToAdd, 3);
testModel.c = [0;0;0;0;0;0;0;1;0;0];%optimize for output flux, if this is used, not sure
testModel.ub = repmat(1000,10,1);
testModel.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
testModel.rxnNames = testModel.rxns;
testModel.b = repmat(0,8,1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TC 0103 - tests both getRxnScoresFromSamples and getExprForRxnScore
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scores = [0.1, -1; ...
          1, 4; ...
          -1, 10; ...
          -3, -5; ...
          8, 4; ...
          10, 6; ...
          6, 7; ...
          2, -1; ...
          1, -4; ...
          -5, 6];

gexTmp =  getExprForRxnScore(scores,1);
gex = [gexTmp;[(10^6 - sum(gexTmp(:,1))), (10^6 - sum(gexTmp(:,2)))]];
sum(gex,1) %ok, TPM


expRes = scores;
expRes([1 2 8],:) = -2 ;

s = Samples();
s.data = gex;
s.genes = [testModel.genes; 'G11'];
s.sampleIds = {'S1','S2'};

rxnScores = getRxnScoresFromSamples(s, testModel, 1);

all(all(abs(expRes - rxnScores) < 10^-10))%ok, they are virtually the same

%all ok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TC 0104 - tests both rxnScoresTo3D and calcJaccard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scores = [0.1, -1, 1,   2, -1, -1, 2, 1; ...
          1,    4, -1,  1, 1,   1, -1, 1; ...
          -1,  10, 1,   1, -1,  1, -1, -1; ...
          -3,  -5, 1,  -1, -1, -1, 1, -1; ...
          8,    4, 2,   2, 1,  -1, 2, 1; ...
          10,   6, 1,  -1, 1,  -1, 1, 1; ...
          6,    7, -1, -1, -1,  1, 1, -1; ...
          2,   -1, 1,  -1, 1,   1, 1, 1; ...
          1,   -4, -1,  -1, -1,  1, -1, 1; ...
          -5,   6, 1,   1, 1,   1, 1, 1];

%allScores = [scores scores scores scores];

allCores3D = rxnScoresTo3D(scores, [1,2]);
all(size(allCores3D) == [10 2 4]) %ok, so 10 reactions, 2 pool sizes, 4 repetitions, which is really 2 repetitions of two pairs.

jacc = calcJaccard(allCores3D);
jaccPair1 = 4/9;
jaccPair2 = 4/8;
jaccPair3 = 3/8;
jaccPair4 = 5/9;

abs(mean([jaccPair1 jaccPair2]) - jacc(1)) < 10^-10 %ok
abs(mean([jaccPair3 jaccPair4]) - jacc(2)) < 10^-10 %ok

%so, since the Jaccard is right, they are most likely in the right order

%all ok


