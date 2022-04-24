%Extracts the task results into a simple matrix
function taskResults = extractTaskResults(filename, filter)
x = load(filename);
tr = x.taskReports;

%assemble into a matrix
mRes = NaN(length(tr{1}.ok), length(tr));
for i = 1:length(tr)
    mRes(:,i) = tr{i}.ok;
end

taskResults = sum(mRes,2);

taskResults = taskResults(filter);

end
