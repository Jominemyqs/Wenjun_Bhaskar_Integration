%========================================================
% File: run_parameter_sweep.m
% Description:
%   Runs Bhaskar's agent-based model at 25 parameter sets 
%   to (hopefully) capture transitions in patterns.
%   JBB, JOO vary in [0,0.05,0.1,0.15,0.2], 
%   while JBO=0.10 is held constant.
%   Then moves output to separate folders and plots final snapshot.
%========================================================

clear all; close all;

% 1) Define the sets for JBB and JOO:
JBBvals = [0.00, 0.05, 0.10, 0.15, 0.20];
JOOvals = [0.00, 0.05, 0.10, 0.15, 0.20];
JBO_val = 0.10;  % fix blue-orange adhesion

% 2) Other simulation parameters:
pol_val        = 0.005;
cell_pop_prop  = [0.4, 0.6];  % 40% orange, 60% blue
n              = 200;         % total # of cells

% 3) Filenames / output directories:
mat_dir         = 'Output_mat';
base_output_dir = 'ParamSweep_Results';

% Ensure those directories exist
if ~exist(mat_dir, 'dir'), mkdir(mat_dir); end
if ~exist(base_output_dir, 'dir'), mkdir(base_output_dir); end

% We'll store 'n' in a .mat so that coculture_model() sees it
matfile_name = fullfile(mat_dir, 'ParamSweep_1.mat');
save(matfile_name, 'n');  % minimal .mat containing n

% 4) Loop over the 5×5 = 25 combos
for iBB = 1:length(JBBvals)
    for iOO = 1:length(JOOvals)
        
        % Unpack the JBB, JOO, JBO from these loops
        JBB_val = JBBvals(iBB);
        JOO_val = JOOvals(iOO);
        % JBO_val is already set above (0.10)
        
        % Prepare struct for model()
        modelpar.model          = 'Bhaskar';
        modelpar.RR            = JOO_val;  % orange-orange
        modelpar.RG            = JBO_val;  % blue-orange
        modelpar.GG            = JBB_val;  % blue-blue
        modelpar.pol_val       = pol_val;
        modelpar.cell_pop_prop = cell_pop_prop;
        
        % Write them to the same ParamSweep_1.mat
        RR = modelpar.RR; 
        RG = modelpar.RG; 
        GG = modelpar.GG; 
        save(matfile_name, 'RR','RG','GG','pol_val','cell_pop_prop','-append');
        
        % Make a nice folder name, e.g. JBB0.10_JOO0.00_JBO0.10
        folder_name = sprintf('JBB%.2f_JOO%.2f_JBO%.2f', JBB_val, JOO_val, JBO_val);
        output_folder = fullfile(base_output_dir, folder_name);
        if ~exist(output_folder, 'dir'), mkdir(output_folder); end
        
        % 5) Run simulation for this param set
        disp(['Running final-snapshot simulation for ', folder_name]);
        model(modelpar);
        
        % 6) Move the newly-created ParamSweep_1_Output => named folder
        if exist('ParamSweep_1_Output','dir')
            movefile('ParamSweep_1_Output', output_folder, 'f');
        else
            warning('ParamSweep_1_Output not found after run. Skipping.');
            continue;
        end
        
        % 7) Quick verification: read positions at final time
        pos_file  = fullfile(output_folder, 'ParamSweep_1_Output', 'Pos_0500000.dat');
        type_file = fullfile(output_folder, 'ParamSweep_1_Output', 'Types_0500000.dat');
        
        if ~exist(pos_file,'file') || ~exist(type_file,'file')
            warning('Missing final position/type files for %s', folder_name);
            continue;
        end
        
        positions = readmatrix(pos_file);
        X = real(positions(:));
        Y = imag(positions(:));
        
        types = readmatrix(type_file);
        if length(types) ~= length(X)
           warning('Dimension mismatch at final snapshot for %s', folder_name);
           continue;
        end
        
        % 8) Make a scatter plot
        figure('visible','off');
        scatter(X, Y, 10, types, 'filled');
        axis equal
        title(sprintf('Final Snapshot: JBB=%.2f, JOO=%.2f, JBO=%.2f', ...
                       JBB_val, JOO_val, JBO_val));
        saveas(gcf, fullfile(output_folder, 'final_pattern.png'));
        close(gcf);
        
        disp(['Done with param set: ', folder_name]);
    end
end

disp('All final snapshots are now generated (25 total).');
