function rxnScores3D = rxnScoresTo3D(rxnScores, xParams)
% rxnScores3D = rxnScoresTo3D(rxnScores, xParams)
% Reshapes a 2D matrix with pairs and repetitions into a 3D matrix
% rxnScores   the 2D matrix
% xParams     pool sizes
% rxnScores3D The 3D matrix

rxns = size(rxnScores,1);
points = length(xParams);
repetitions = size(rxnScores,2)/points;
if round(repetitions) ~= repetitions
   error('rxnScoresTo3D: The rxnScores and xParams does not match'); 
end

rxnScores3D = zeros(rxns, points, repetitions);

for point = 1:points
    for rep = 1:repetitions
        index = (point-1)*repetitions + rep;
        rxnScores3D(:,point, rep) = rxnScores(:,index);
    end
end

end



