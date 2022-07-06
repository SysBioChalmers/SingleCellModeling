%This code compares the output of the old and new tINIT versions on a small test model.
%The code is used for generating a supplementary figure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate a test model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode
testModelL = struct();
testModelL.id = 'testModel';
testModelL.rxns = {};
testModelL.S=[];
testModelL.rev=[];
testModelL.metNames = {'e1';'e2';'e3';'e4';'e5';'e6';'e7';'e8';'e9';'e1';'e2';'e3';'e4';'e5';'e6';'e7';'e8';'e9';'x1';'x2';'x3';'x4';'x5';'x6';'x7';'x8';'x9';'x10';'x11'};
testModelL.comps = {'s';'c'};
testModelL.compNames = testModelL.comps;
testModelL.metComps = [1;1;1;1;1;1;1;1;1;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2];
testModelL.mets = strcat(testModelL.metNames, testModelL.comps(testModelL.metComps));


testModelL.grRules = {};
testModelL.rxnGeneMat = [];

testModelL.genes = {'Ge1';'Ge2';'Ge4';'Ge5';'Ge7';'Ge9'; 'Gr1';'Gr2';'Gr3';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'Gr14';'Gr15'};


testModelL.ub = [];
testModelL.lb = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {  'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2';'E2b';'E3';'E4';'E5';'E6';'E7';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10';'R11';'R12';'R13';'R14';'R15'};
rxnsToAdd.grRules = {'';  '';  '';  '';  '';  '';  '';  '';  ''; 'Ge1';'Ge2';'';'';'Ge4';'Ge5';'';'Ge7';'';'Ge9'; 'Gr1';'Gr2';'Gr3';'';'Gr5';'Gr6';'Gr7';'Gr8';'Gr9';'Gr10';'Gr11';'Gr12';'';'Gr14';'Gr15'};
rxnsToAdd.equations = {'e1[s] <=>';...
                       'e2[s] <=>';...
                       'e3[s] <=>';...
                       'e4[s] <=>';...
                       'e5[s] <=>';...
                       'e6[s] <=>';...
                       'e7[s] <=>';...
                       'e8[s] <=>';...
                       'e9[s] <=>';...
                       'e1[s] <=> e1[c]';...
                       'e2[s] <=> e2[c]';...
                       'e2[s] <=> e2[c]';... %b variant
                       'e3[s] <=> e3[c]';...
                       'e4[s] <=> e4[c]';...
                       'e5[s] <=> e5[c]';...
                       'e6[s] <=> e6[c]';...
                       'e7[s] <=> e7[c]';...
                       'e8[s] <=> e8[c]';...
                       'e9[s] <=> e9[c]';...
                       'e1[c] + e2[c] <=> x1[c]';... %R1
                       'e1[c] + e3[c] => x2[c] + x3[c]';... %R2
                       'e4[c] + x3[c] => x4[c] + x5[c]';... %R3
                       'e5[c] + e6[c] + x4[c] => 2 x2[c] + x6[c]';... %R4
                       'x1[c] + x2[c] <=> x7[c] + 2 x8[c]';... %R5
                       'x2[c] + x8[c] => x3[c] + x9[c]';... %R6
                       'x4[c] <=> x9[c]';... %R7
                       'x5[c] <=> x9[c]';... %R8
                       'x6[c] <=> x10[c]';... %R9
                       'x6[c] <=> x11[c]';... %R10
                       'x10[c] + 2 x11[c] => e7[c]';... %R11
                       'x9[c] + x10[c] <=> e8[c]';... %R12
                       'x7[c] + x8[c] + x9[c] => e9[c]';... %R13
                       'x6[c] => x9[c]';... %R14
                       'x3[c] => x9[c]'... %R15
                       };
testModelL = addRxns(testModelL,rxnsToAdd, 3);
testModelL.c = [0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];%optimize for output flux, if this is used, not sure
testModelL.rxnNames = testModelL.rxns;
testModelL.b = repmat(0,length(testModelL.mets),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define test model gene scores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testModelLGeneScores = [3; -1;   8;    6;    -5;    5;     4;    5;    2;    3;    6;    1;    3;    1;    -3;    1;     3;      1;    2];

testParams = struct();

arrayDataL = struct();
arrayDataL.genes = testModelL.genes;
arrayDataL.tissues = {'t1'};
arrayDataL.levels = getExprForRxnScore(testModelLGeneScores,1);
arrayDataL.threshold = 1;

%Run prep data
prepDataL = prepINITModel(testModelL, [], {}, false, {}, 's');

mres = ftINIT(prepDataL,arrayDataL.tissues{1},[],[],arrayDataL,[],getINITSteps(),true,true,testParams);
mres2 = ftINIT(prepDataL,arrayDataL.tissues{1},[],[],arrayDataL,[],getINITSteps([], 'full'),true,true,testParams);

%expResult = {  'S1';'S2';'S3';'S4';'S5';'S6';'S7';'S8';'S9';'E1';'E2';'E3';'E4';'E5';'E6';'E8';'E9';'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R12';'R13';'R14';'R15'};


%run the old tINIT version (in Human-GEM)
paramsL2 = struct();
paramsL2.TimeLimit = 1000;
testModelL2 = closeModel(testModelL);
init_modelOrig = getINITModel2(testModelL2,arrayDataL.tissues{1},[],[],arrayDataL,[],true,[],true,true,[],paramsL2);
%run a modified version of the old tINIT (in Human-GEM), where secretion of metabolites and fluxes in both directions for reversible reactions are not allowed
%This is "ground truth"
init_modelOrigNoSecrOneDirOnly = getINITModel2RunFull(testModelL2,arrayDataL.tissues{1},[],[],arrayDataL,[],true,[],true,true,[],paramsL2);

table(testModelL.rxns, ismember(testModelL.rxns, mres.rxns), ismember(testModelL.rxns, mres2.rxns), ismember(testModelL.rxns, init_modelOrig.rxns), ismember(testModelL.rxns, init_modelOrigNoSecrOneDirOnly.rxns))
