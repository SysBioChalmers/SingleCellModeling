cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define test models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
rxnsToAdd.grRules = {'';'';'G3';'G4';'G5';'G6';'G7';'';'G9';'G10'};;
testModel = addRxns(testModel,rxnsToAdd, 3);
testModel.c = [0;0;0;0;0;0;0;1;0;0];%optimize for output flux, if this is used, not sure
testModel.ub = repmat(1000,10,1);
testModel.lb = [0;-1000;-1000;-1000;0;0;0;0;-1000;-1000];
testModel.rxnNames = testModel.rxns;
testModel.b = repmat(0,8,1);

testRxnScores = [-2;-2;-1;7;0.5;0.5;-1;-2;-1;1.5];

detectedMets = {};
params.TimeLimit = 10;

testModelTasks = struct();
testModelTasks.id = 'Gen e[s] from a[s]';
testModelTasks.description = 'Gen e[s] from a[s]';
testModelTasks.shouldFail = false;
testModelTasks.printFluxes = false;
testModelTasks.comments = '';
testModelTasks.inputs = {'a[s]'};
testModelTasks.LBin = 0;
testModelTasks.UBin = inf;
testModelTasks.outputs = {'e[s]'};
testModelTasks.LBout = 1;
testModelTasks.UBout = 1;
testModelTasks.equations = {};
testModelTasks.LBequ = [];
testModelTasks.UBequ = [];
testModelTasks.changed = {};
testModelTasks.LBrxn = {};
testModelTasks.UBrxn = {};



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testModel2 = struct();
testModel2.id = 'testModel2';
testModel2.rxns = {};
testModel2.S=[];
testModel2.rev=[];
testModel2.mets = {'a';'b'};
testModel2.metNames = {'a';'b'};
testModel2.comps = {'s'};
testModel2.compNames = testModel2.comps;
testModel2.metComps = [1;1];
testModel2.genes = {'G1';'G2';'G3';'G4'};
testModel2.grRules = {};
testModel2.rxnGeneMat = [];

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'R3';'R4'};
rxnsToAdd.equations = {'a[s] <=>';...
                       'a[s] => b[s]';...
                       'a[s] <=> b[s]';...
                       'b[s] =>'};
rxnsToAdd.grRules = testModel2.genes;
testModel2 = addRxns(testModel2,rxnsToAdd, 3,true,true);
testModel2.c = [0;0;0;1];%optimize for output flux, if this is used, not sure
testModel2.ub = repmat(1000,4,1);
testModel2.lb = [-1000;0;-1000;0];
testModel2.rxnNames = testModel2.rxns;
testModel2.b = zeros(2,1);
%testRxnScores2 = [-1.1;-1;8;-1];
%testRxnScores2 = [-1;-1;8;-1];
testRxnScores2 = [-1.1;-1;8;-1];
testRxnScores2b = [-1.1;8;-1;-1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testModel3 = testModel2;

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R5';'R6';'R7'};
rxnsToAdd.equations = {'a[s] <=> 0.5 c[s] + 0.5 d[s]';...
                       'c[s] <=> b[s]';...
                       'd[s] <=> b[s]'};
rxnsToAdd.grRules = {'G5';'G6';'G7'};
testModel3 = addRxns(testModel3,rxnsToAdd, 3, [], true, true);
testRxnScores3 = [-1;-1;2;-1;0.5;-2;1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testModel4 = testModel2;

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R5';'R6';'R7';'R8';'R9';'R10';'R11'};
rxnsToAdd.equations = {'5 a[s] <=> 5 d[s]';...
                       'e[s] <=> d[s]';
                       'f[s] + g[s] <=> e[s]';
                       'b[s] <=> f[s]';
                       'h[s] <=> g[s]';
                       'h[s] =>';
                       'e[s] => g[s]'};
rxnsToAdd.grRules = {'G5';'G6';'G7';'G8';'G9';'G10';'G11'};
testModel4 = addRxns(testModel4,rxnsToAdd, 3, [], true, true);
testRxnScores4 = [-1;-1;2;-1;0.5;-2;1;1.3;-0.5;-0.4;8];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testModel5 - for testing functions for separating the problem into subgroups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

testModel5 = testModel2;

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R5';'R6';'R7';'R8';'R9';'R10';'R11';'R12'};
rxnsToAdd.equations = {'5 c[s] <=> 5 d[s]';...
                       'e[s] <=> f[s]';
                       'f[s] + g[s] <=> h[s]';
                       'c[s] <=>';
                       'd[s] <=>';
                       'e[s] <=>';
                       'f[s] <=>';
                       'g[s] <=>'};
rxnsToAdd.grRules = {'G5';'G6';'G7';'G8';'G9';'G10';'G11';'G12'};
testModel5 = addRxns(testModel5,rxnsToAdd, 3, [], true, true);
testRxnScores5 = [-1;-1;2;-1;0.5;-2;1;1.3;-0.5;-0.4;8;8];
constructEquations(testModel5)


testParams = struct();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test model 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testModel5 = testModel2;
testModel5 = removeReactions(testModel5, testModel5.rxns);
testModel5.comps = {'c'};

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R8';'R9';'R10'};
rxnsToAdd.equations = {'a[c] => b[c]';...
                       'a[c] <=> b[c]';...
                       'a[c] <=> b[c]';...
                       'a[c] <=> b[c]';...
                       'b[c] <=> c[c]';...
                       'b[c] <=> c[c]';...
                       'b[c] <=> c[c]';...
                       'b[c] <=> c[c]';...
                       'c[c] <=>';...
                       'a[c] <=>'};
rxnsToAdd.grRules = {'G1';'G2';'G3';'G4';'G5';'G6';'G7';'G8';'G9';'G10'};
testModel5 = addRxns(testModel5,rxnsToAdd, 3, [], true, true);
testRxnScores5 = [1;-1;1;1;1;1;1;1;0;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test model 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testModel6 = testModel2;
testModel6 = removeReactions(testModel6, testModel6.rxns);
testModel6.comps = {'c'};

rxnsToAdd = struct();
rxnsToAdd.rxns = {'R1';'R2';'R3';'R4';'R5'};
rxnsToAdd.equations = {'a[c] => b[c]';...
                       'a[c] <=> b[c]';...
                       'a[c] <=> b[c]';...
                       'a[c] <=>';...
                       'b[c] <=>'};
rxnsToAdd.grRules = {'G1';'G2';'G3';'G4';'G5'};
testModel6 = addRxns(testModel6,rxnsToAdd, 3, [], true, true);
testRxnScores6 = [1;-1;1;0;0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cases - first tests of the full function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0001: testModel without tasks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prepDataTest1 = prepINITModel(testModel, {}, {}, false);
%check some things in the prepData
%1. We expect 3 rxns in origRxnsToZero:
all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch | prepDataTest1.toIgnoreImportRxns) , {'R1';'R2';'R8'})) %ok
%note that R7 should not be there, since it has a GPR.

arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(testRxnScores);
arrayData1.threshold = 1;
tst1ResModel1 = getINITModel9(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[], true,false,testParams);

%We expect R1, R2, and R8 to be added since they have no GPRs and are exch/simple transport rxns
%R7 however will not be added since it has a GPR
all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R8';'R9';'R10'})) %ok

%also test spontaneous
prepDataTest1 = prepINITModel(testModel, {}, {'R7';'R10'}, false);
all(strcmp(prepDataTest1.refModel.rxns(prepDataTest1.toIgnoreExch | prepDataTest1.toIgnoreImportRxns | prepDataTest1.toIgnoreSpont), {'R1';'R2';'R7';'R8';'R10'})) %ok
arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(testRxnScores);
arrayData1.threshold = 1;
tst1ResModel1 = getINITModel9(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[], true,false,testParams);
%the model should now change to include the "correct" path and 
%skip R9/R10:
all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8'})) %ok


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0002: Create a task that wants to generate e[s] from a[s] for testModel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prepDataTest1 = prepINITModel(testModel, testModelTasks, {}, false);
%We now expect to R2 and R7 to be essential. Note that R1 and R8 are not essential,
%the exchange rxns are not used when checking tasks.
%This is a bit complicated to check, because the essential rxns are expressed
%as rxn ids of the linearly merged model. We expect those to be called R1 and R7:
all(strcmp(prepDataTest1.essentialRxns,{'R1';'R7'})) %ok

arrayData1.genes = testModel.genes;
arrayData1.tissues = {'a'};
arrayData1.levels = getExprForRxnScore(testRxnScores);
arrayData1.threshold = 1;
tst1ResModel1 = getINITModel9(prepDataTest1,arrayData1.tissues{1},[],[],arrayData1,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[], true,false,testParams);
%Since both R2 and R7 are now essential, we expect all rxns to be on except R3 and 
%R5 (which have a negative total score and are not needed for the task)
all(strcmp(tst1ResModel1.rxns,{'R1';'R2';'R4';'R6';'R7';'R8';'R9';'R10'})) %ok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0003: The second step - gapfilling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First generate a model with gaps. We can use testModel. First we add 
%boundary mets. Then we remove the exchange reactions 
%(which is required for filling gaps) and create a gap by removing the R7 reaction.
mTempRef = addBoundaryMets(testModel);
mTempRef = removeReactions(mTempRef, {'R1';'R8'});
mTemp = removeReactions(mTempRef, {'R7'});
mTemp.id = 'tmp';
tmpRxnScores = testRxnScores([2;3;4;5;6;7;9;10]);
%now check that R7 is added back
[outModel,addedRxnMat] = fitTasksOpt(mTemp,mTempRef,[],true,min(tmpRxnScores,-0.1),testModelTasks);
strcmp(mTempRef.rxns(addedRxnMat),'R7')%ok

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0004: MergeLinear and groupRxnScores
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first testModel
[reducedModel,origRxnIds,groupIds,reversedRxns]=mergeLinear(testModel, {});
%We expect mergeLinear to merge {R1,R2}, {R3,R5}, {R4,R6}, {R7,R8}, {R9,R10}
all(groupIds == [1;1;2;3;2;3;4;4;5;5])%ok
%we expect R1, R3, R4, R7 to be irreversible, R9 to be reversible
all(reducedModel.rev == [0;0;0;0;1]) %ok
all(reducedModel.lb == [0;0;0;0;-1000])

newRxnScores=groupRxnScores(reducedModel, testRxnScores, origRxnIds, groupIds, ismember(origRxnIds, {'R1';'R2';'R8'}));
all(newRxnScores == [0;-0.5;7.5;-1;0.5])%ok


%then testModel4
[reducedModel,origRxnIds,groupIds,reversedRxns]=mergeLinear(testModel4, {});
%we expect {R5,R6},{R7,R8}, and{R9,R10} to be merged
all(groupIds == [0;0;0;0;1;1;2;2;3;3;0])%ok
%check reversibility
all(reducedModel.rev == [1;0;1;0;1;1;0;0])%ok
%check that some reactions have flipped direction when turned to irrev
strcmp(constructEquations(reducedModel, 'R9'),'g[s] => ')%ok
all(find(reversedRxns) == [6;9])%ok

%constructEquations(testModel4)
%constructEquations(reducedModel)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0006: reverseRxns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%R1 = '=> a[s]';...
%R3 = 'a[c] <=> b[c] + c[c]';...

tmpModel = reverseRxns(testModel, {'R1';'R3'});
res = constructEquations(tmpModel, {'R1';'R3'});
expRes = {'a[s] => ';'b[c] + c[c] <=> a[c]'}
all(strcmp(res,expRes)) %ok


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%T0007: rescaleModelForINIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
miniModel = struct();
miniModel.S = [1,1000;-1,-1000];
res = rescaleModelForINIT(miniModel,10);
expRes = [1,10;-1,-10];
all(all(res.S == expRes)) %ok


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% T0050 - Build a more complex model with transport etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is not a formal test case that can be easily run, but it has been run to confirm that the method works well on a larger model
% It requires manual modification of the code in Human-GEM, which makes it a bit impractical to run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
testModelLGeneScores = [3; -1;   8;    6;    -5;    5;     4;    5;    2;    3;    6;    1;    3;    1;    -3;    1;     3;      1;    2];


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


arrayDataL = struct();
arrayDataL.genes = testModelL.genes;
arrayDataL.tissues = {'t1'};
cd 'C:\Work\MatlabCode\projects\HMASandbox\HMA_Sandbox\Single-cell Modeling\tINIT'
arrayDataL.levels = getExprForRxnScore(testModelLGeneScores,1);
arrayDataL.threshold = 1;

%Run prep data
prepDataL = prepINITModel(testModelL, [], {}, false, {});

paramsL = struct();
paramsL.TimeLimit = 120;
paramsL.MIPGap = 0.0004;

mres = getINITModel9(prepDataL,arrayDataL.tissues{1},[],[],arrayDataL,[],[1;1;1;1;1;1;1;0],[1;1;1;1;1;1;1;0],true,true,[],true,false,paramsL);
mres2 = getINITModel9(prepDataL,arrayDataL.tissues{1},[],[],arrayDataL,[],[1;1;1;1;1;1;1;0],[1;0;0;0;1;1;1;0],true,true,[],true,false,paramsL);

%run the old tINIT version
paramsL2 = struct();
paramsL2.TimeLimit = 1000;
testModelL2 = addBoundaryMets(testModelL);
init_modelOrig = getINITModel2(testModelL2,arrayDataL.tissues{1},[],[],arrayDataL,[],true,[],true,true,[],paramsL2);

%in this call, I have modified the code - the possibility to turn off met secretion + don't allow flux in both directions is not possible.
%The following line, around line 337, is changed
%from:
%[~, deletedRxnsInINIT, metProduction] = runINIT(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,true,false,params);
%to:
%[~, deletedRxnsInINIT, metProduction] = runINIT(simplifyModel(cModel),rxnScores,metabolomicsData,essentialRxnsForTasks,0,false,true,params);
init_modelOrigNoSecrOneDirOnly = getINITModel2(testModelL2,arrayDataL.tissues{1},[],[],arrayDataL,[],true,[],true,true,[],paramsL2);

%The models init_modelOrigNoSecrOneDirOnly and mres2 are very similar, (only one exch rxn differ, which is expected) 
%init_modelOrig is quite different, with a lot of gaps, and worse. So, the conclusion is that the new version does a pretty good job.

