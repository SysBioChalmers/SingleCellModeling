function meanRes = calcJaccard(rxnScores, progrBarCtxt)

if nargin < 2
    progrBarCtxt = [];
end

progbar = ProgrBar('calcJaccard', progrBarCtxt);

points = size(rxnScores,2);
pairs = floor(size(rxnScores,3)/2);

res = zeros(points, pairs);

for point = 1:points
    for pair = 1:pairs
        repInd = 2*pair-1;
        rxnScoresA = rxnScores(:, point, repInd);
        rxnScoresB = rxnScores(:, point, repInd+1);
        jc = jaccard(rxnScoresA >= 0, rxnScoresB >= 0);
        res(point, pair) = jc;
    end
    progbar.Progress(point/points);
end

meanRes = mean(res,2);

progbar.Done;

end
