%Generates random collections of cells that are pooled to pseudo-bulk samples.
function [samples, poolSizes] = genPoolSizeSamples(ds, minSizeSet, maxSizeSet, numPoints, repetitions, progrBarCtxt)
%This function generates paired samples next to each other. If we have 2 repetitions and three points, the samples will 
%be distributed as: p1r1A p1r1B p1r2A p1r2B p2r1A p2r1B p2r2A p2r2B p3r1A p3r1B p3r2A p3r2B 
%where A and B can be compared for each repetition (r) within each point (p)
if nargin < 6
    progrBarCtxt = [];
end

ds_red = TPM(ds);
totset = ds_red.data;
totNumCells = size(totset,2);


progbar = ProgrBar(strcat('getPoolSizeSamples:',ds.name), progrBarCtxt);

steps = numPoints;
maxSize = min(maxSizeSet, floor(totNumCells/2));

%spread the x:es out evenly in log space between min and max
poolSizes = round(10 .^ (linspace(log10(minSizeSet),log10(maxSize),steps)));

samples = Samples();
samples.genes = ds.genes;
samples.data = zeros(numel(ds.genes), numPoints*repetitions*2);
samples = samples.fillEmpties();
totNumCellsVec = 1:totNumCells;

for sz = 1:steps
    numSamp = poolSizes(sz);
    for repetition = 1:repetitions
        indexBase = (sz-1)*repetitions*2;
        aInd = repetition*2 - 1 + indexBase;
        bInd = aInd + 1;
        as = randsample(totNumCells, numSamp);
        lgc = false(totNumCells,1);
        lgc(as) = true;
        indLeft = totNumCellsVec(~lgc);
        bs = indLeft(randsample(totNumCells - length(as), numSamp));
        
        samples.data(:,aInd) = mean(totset(:,as),2);
        samples.data(:,bInd) = mean(totset(:,bs),2);
        progbar.Progress((sz-1 + repetition/repetitions)/steps);
    end
    progbar.Progress(sz/steps);
end

progbar.Done();

end
