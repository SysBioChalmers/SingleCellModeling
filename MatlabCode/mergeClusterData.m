%just assembles the depmap15 chunk model files from the old tINIT run on the cluster into one file
cd 'C:/Work/MatlabCode/projects/SingleCellModeling/SingleCellModeling/MatlabCode'

%% DepMap15
%%%%%%%%%%%%%%
depmap15_init_models = {};
for i = 1:10
    fn = ['../data/DepMap15FromCluster/init_models_depmap15-' num2str(i) '.mat'];
    x = load(fn);
    depmap15_init_models = [depmap15_init_models;x.init_models_depmap];
end
save('../data/init_models_depmap15.mat','depmap15_init_models');





