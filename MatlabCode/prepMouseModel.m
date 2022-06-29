%Runs prepInitModel for Mouse-GEM, including filtering of reactions.
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%load Mouse-GEM and tasks
load('C:/Work/MatlabCode/components/MouseGEM/Mouse-GEM_1_3/Mouse-GEM/model/Mouse-GEM.mat')

%first filter some duplicate reactions in the model to speed up the calculations and aviod randomness.
m = removeUnwantedReactions(mouseGEM);%We use an older version of mouse-GEM, so this operation should be safe here
%for some reason the b vector has two columns in the mouse model - change that to 1
m.b = zeros(length(m.b),1);

%Remove a lot of unnecessary fields - the models become too large and cannot be saved otherwise
m = rmfield(m,'metBiGGID');
m = rmfield(m,'metChEBIID');
m = rmfield(m,'metEHMNID');
m = rmfield(m,'metHMDBID');
m = rmfield(m,'metHMR2ID');
m = rmfield(m,'metHepatoNET1ID');
m = rmfield(m,'metLipidMapsID');
m = rmfield(m,'metMetaNetXID');
m = rmfield(m,'metPubChemID');
m = rmfield(m,'metRecon3DID');
m = rmfield(m,'metKEGGID');
m = rmfield(m,'metRetired');

m = rmfield(m,'rxnBiGGID');
m = rmfield(m,'rxnEHMNID');
m = rmfield(m,'rxnHMR2ID');
m = rmfield(m,'rxnHepatoNET1ID');
m = rmfield(m,'rxnKEGGID');
m = rmfield(m,'rxnMetaNetXID');
m = rmfield(m,'rxnREACTOMEID');
m = rmfield(m,'rxnRatconID');
m = rmfield(m,'rxnRecon3DID');
m = rmfield(m,'rxnRetired');
m = rmfield(m,'rxnRheaID');
m = rmfield(m,'rxnRheaMasterID');
m = rmfield(m,'rxnTCDBID');

m = rmfield(m,'rxnMiriams');
m = rmfield(m,'metMiriams');
m = rmfield(m,'inchis');
m = rmfield(m,'rxnReferences');
m = rmfield(m,'rxnNotes');



%This is a little bit tricky - this will use the spontaneous reactions from Human-GEM, which will also work.
prepDataMouseGEM = prepHumanModelForftINIT(m, false);
save('../data/prepDataMouseGEM.mat','prepDataMouseGEM');
