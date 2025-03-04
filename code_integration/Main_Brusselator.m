% test model

clear all;
close all;

path(path,'..');

feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

modelpar = struct;
modelpar.model = 'Brusselator';
modelpar.init = 'random';
modelpar.sets = 'pos';
modelpar.xdim = 2;
modelpar.threshold = 0.7;

featpar = struct;
featpar.feature = 'num';
featpar.alpha = 1;
featpar.avoid_steady = 0;  % Add this line


contpar = struct;
contpar.step_size = 0.5;
contpar.max_step_arc_length = 0.5;
contpar.min_step_arc_length = contpar.max_step_arc_length/2;
contpar.L = 8;
contpar.max_sim = 10;
contpar.stopping = 1e-3;  % Add this line

reproduce = 'steady';

switch reproduce
    case 'pos'
        start.point = [3.5,6.5];
        start.normal = [0,1];
        modelpar.sets = 'pos';
        modelpar.threshold = 0.7;


    case 'neg'

        start.point = [4.2,6.5];
        start.normal = [0,1];
        modelpar.sets = 'neg';
        modelpar.threshold = 0.3;

    case 'steady'
        start.point = [4.5,6.5];
        start.normal = [0,1];
        featpar.feature = 'steady';
        

end



cd ..
[p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
cd ReproduceCurves