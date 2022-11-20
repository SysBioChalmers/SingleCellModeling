%Assembles the data for Fig 2C and some associated supplementary plots
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

numPairs = 60;
numSeries = 5;
numPS = 8;

fullVals = cell(numSeries,1);

jaccardData = nan(numSeries,numPS);
for series = 1:numSeries
    fullValsSer = nan(numPS,numPairs);
    for ps = 1:numPS
        fn = ['../data/contamination/cont_models_' num2str(series) '_' num2str(ps) '.mat'];
        modRes = load(fn).rxnContent;
        jaccVals = nan(numPairs,1);
        for i = 1:numPairs
            ind1 = i*2 - 1;
            ind2 = ind1 + 1;
            jaccVals(i) = jaccard(modRes(:,ind1),modRes(:,ind2));
        end
        jaccardData(series,ps) = mean(jaccVals);
        fullValsSer(ps,:) = jaccVals.';
    end
    fullVals{series} = fullValsSer;
end

%read the pool sizes
load('../data/contamination/PSPoolSizesCont.mat');

s = struct();
s.figData = jaccardData;
s.seriesNames = [0 .02 .05 .1 .2];
s.poolSizes = poolSizesCont;
s.fullVals = fullVals;
save('../data/contamination/PSContResults.mat','s')

