%Same as genPoolSizeSamples, but draws from two populations, a certain fraction from each
function [samplesList, poolSizes] = genPoolSizeSamplesWithContamination(ds, dsCont, minSizeSet, maxSizeSet, numPoints, repetitions, contFrac, progrBarCtxt)

if nargin < 8
    progrBarCtxt = [];
end

%check that the genes are synchronized
if (length(ds.genes) ~= length(dsCont.genes)) | ~all(strcmp(ds.genes,dsCont.genes))
    error('The genes between ds and dsCont must match!')
end

ds_red = TPM(ds);
totset = ds_red.data;

dsc_red = TPM(dsCont);
totsetC = dsc_red.data;

totNumCells = size(totset,2);
totNumCellsC = size(totsetC,2);


progbar = ProgrBar(strcat('getPoolSizeSamples:',ds.name), progrBarCtxt);

steps = numPoints;
maxSize = min(maxSizeSet, floor(min(totNumCells/2, totNumCellsC/max(contFrac))));

%spread the x:es out evenly in log space between min and max
poolSizes = round(10 .^ (linspace(log10(minSizeSet),log10(maxSize),steps)));

samplesList = cell(length(contFrac),1);

contSteps = length(contFrac);
totNumCellsVec = 1:totNumCells;

for cont = 1:contSteps
    samplesList{cont} = Samples();
    samplesList{cont}.genes = ds.genes;
    samplesList{cont}.data = zeros(numel(ds.genes), numPoints*repetitions*2);
    samplesList{cont} = samplesList{cont}.fillEmpties();
    
    for sz = 1:steps
        numSamp = poolSizes(sz);
        numSampOrig = round(numSamp*(1-contFrac(cont)));
        numSampCont = numSamp - numSampOrig;
        for repetition = 1:repetitions
            indexBase = (sz-1)*repetitions*2;
            aInd = repetition*2 - 1 + indexBase;
            bInd = aInd + 1;
            as = randsample(totNumCells, numSamp);
            lgc = false(totNumCells,1);
            lgc(as) = true;
            indLeft = totNumCellsVec(~lgc);
            
            bsOrig = indLeft(randsample(totNumCells - length(as), numSampOrig));
            if numSampCont > 0
                bsCont = randsample(totNumCellsC, numSampCont);
            end
            

            samplesList{cont}.data(:,aInd) = mean(totset(:,as),2);
            if numSampCont > 0
                samplesList{cont}.data(:,bInd) = (mean(totset(:,bsOrig),2) .* numSampOrig + mean(totsetC(:,bsCont),2) .* numSampCont)./numSamp;
            else
                samplesList{cont}.data(:,bInd) = (mean(totset(:,bsOrig),2) .* numSampOrig)./numSamp;
            end
            progbar.Progress((sz-1 + repetition/repetitions)/steps/contSteps + (cont-1)/contSteps);
        end
        progbar.Progress(sz/steps/contSteps + (cont-1)/contSteps);
    end
    progbar.Progress(cont/contSteps);
end
progbar.Done();

end
