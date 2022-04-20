function expr=getExprForRxnScore(scores, threshold)
%this function returns the expression that will give the scores sent in when sent to scoreComplexModel
%this is useful for writing test code
if nargin < 2
    threshold = 1;
end

%This is how the score is calculated: 5*log(expression./threshold)
%expression = threshold*10.^(scores/5)

expr = threshold*exp(scores/5);

end
