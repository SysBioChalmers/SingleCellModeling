%Gets reaction scores for each sample in a samples object
function rxnScores = getRxnScoresFromSamples(s, model, tINITThreshold, progrBarCtxt)

if nargin < 4
    progrBarCtxt = [];
end

if nargin < 3
    tINITThreshold = 1;
end

s = TPM(s);


progbar = ProgrBar('getRxnScoresFromSamples', progrBarCtxt);

rxnScores = zeros(size(model.rxns,1), numel(s.sampleIds));

%prepare array data
arrayData = struct();
arrayData.genes = s.genes;
arrayData.tissues = {'1'};
arrayData.celltypes = {'1'};
arrayData.threshold = tINITThreshold;
progbar.Progress(0);
for i = 1:numel(s.sampleIds)
    arrayData.levels = s.data(:, i);
    [rxnScores(:,i), ~, ~, ~] = scoreComplexModel(model, [], arrayData, '1', '1');
    progbar.Progress(i/numel(s.sampleIds));
end

progbar.Done();

end



