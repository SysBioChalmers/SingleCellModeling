function rxnScores3D = rxnScoresTo3D(rxnScores, xParams)
%x params here is for example pool sizes or noise magnitude

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



