% test model

clear all;
close all;

modelpar = struct;

modelpar.model = 'Barkley';
modelpar.xdim = 2;
modelpar.threshold = 0.7;


featpar = struct;
featpar.feature = 'retract';
modelpar.sets = 'pos';
featpar.N = 1;
featpar.alpha = 1;
featpar.avoid_steady = 0;  % Or set to 1 if needed  


reproduce = 'retract';

contpar = struct;
contpar.step_size = 0.01;
contpar.max_step_arc_length = 0.01;
contpar.min_step_arc_length = contpar.max_step_arc_length/2;
contpar.L = .5;
contpar.max_sim = 1;
contpar.stopping = 0;  % Add this to set the stopping threshold

switch reproduce
    case 'retract'
        % Find retract
        featpar.feature = 'retract';
        start.point = [0.375,0.015];
        start.normal = [0.01,0.01];
        modelpar.sets = 'pos';
    case 'meander1'
        start.point = [0.40,0.02];%+0.02;
        start.normal = [0.01,0.01];
        featpar.feature = 'meander';
        modelpar.sets = 'tip_points';
    case 'meander2'
        start.point = [0.66,0.01];%+0.02;
        start.normal = [0.01,0.01];
        featpar.feature = 'meander';
        modelpar.sets = 'tip_points';

    case 'drift'
        start.point = [0.45,0.005];
        start.normal = [0.01,0.01];
        contpar.max_step_arc_length = 0.005*1.5;
        contpar.min_step_arc_length = 0.005*1.5;
        contpar.step_size = 0.005*1.5;
        featpar.feature = 'drift';
        modelpar.sets = 'tip_points';




end









feature_handle = @feature_evaluation;
obj_handle = @objective_evaluation;

cd ..
[p_history,counts,L_history,p_history_all,metric_history_all] = continuation(contpar, featpar, modelpar, feature_handle, obj_handle, start )
cd ReproduceCurves