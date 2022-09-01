%This function is copied from the TME modeling project
function res = listMetRxnsWithFluxes(model, sol, met, consumption, lowerLimit)
if nargin < 5
    lowerLimit = 0;
end

%first check if the met is specific for a certain compartment
startIndex = regexp(met,'\[[a-zA-Z_0-9\-]+\]');
startIndex = startIndex(length(startIndex)); %the compound can look like this, so pick the last one: 'heptadecanoyl-[ACP][c]'
if ~isempty(startIndex)
    compStr = extractBetween(met, startIndex+1, strlength(met) - 1);
    comp = find(strcmp(model.comps, compStr));
    if isempty(comp) %if the compound looks like this, i.e., the [ACP] is not a compartment, list for all compartments: 'heptadecanoyl-[ACP]'
        mets = strcmp(model.metNames,met);
    else
        met = extractBetween(met, 1, startIndex-1);
        mets = strcmp(model.metNames,met) & (model.metComps == comp);
    end
else
    mets = strcmp(model.metNames,met);
end
SMultFlux = model.S(mets,:) .* sol.x.';
if (consumption)
    rxns = sum(SMultFlux < -lowerLimit, 1) > 0;
    SMultFluxNeg = SMultFlux(:,rxns);
    SMultFluxNeg(SMultFluxNeg > 0) = 0;
    fluxes = sum(-SMultFluxNeg,1);
else
    rxns = sum(SMultFlux > lowerLimit, 1) > 0;
    SMultFluxPos = SMultFlux(:,rxns);
    SMultFluxPos(SMultFluxPos < 0) = 0;
    fluxes = sum(SMultFluxPos,1);
end
tbl = table(model.rxns(rxns), fluxes.', constructEquations(model, model.rxns(rxns)));
tbl.Properties.VariableNames = {'Rxn', 'Metabolite Flux', 'Formula'};
tbl
res = tbl;

end
