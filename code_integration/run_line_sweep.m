%========================================================
% File: run_line_sweep.m
% Purpose:
%   Sweep along a 1‑D line in (JBB,JOO,JBO) space:
%      • JBB fixed  = 0.05
%      • JBO fixed  = 0.10
%      • JOO varies = 0 : h : JOOmax
%   This is Step‑1 of Wenjun’s “Initialization” paragraph:
%   generate snapshots we will later feed into the Wasserstein test.
%========================================================

clear; close all; clc;

% -------- line parameters ---------------------------------------------
JBB_val   = 0.05;          % fixed blue‑blue
JBO_val   = 0.10;          % fixed blue‑orange
h         = 0.005;         % step‑size along JOO   (choose 0.002‑0.005)
JOO_max   = 0.05;          % end of the line sweep
JOO_vals  = 0 : h : JOO_max;   % vector of JOO points
% -----------------------------------------------------------------------

% other model constants (unchanged)
pol_val        = 0.005;
cell_pop_prop  = [0.4, 0.6];
n              = 200;

% folders
mat_dir         = 'Output_mat';
base_output_dir = 'LineSweep_JBB0p05';     % keep separate from 5×5 grid

if ~exist(mat_dir,'dir'),         mkdir(mat_dir); end
if ~exist(base_output_dir,'dir'), mkdir(base_output_dir); end

% minimal .mat so coculture_model sees total population n
save(fullfile(mat_dir,'ParamSweep_1.mat'),'n');

% -----------------------------------------------------------------------
for k = 1:length(JOO_vals)
    
    JOO_val = JOO_vals(k);
    
    % build modelpar struct
    modelpar.model          = 'Bhaskar';
    modelpar.RR             = JOO_val;   % orange‑orange
    modelpar.RG             = JBO_val;   % blue‑orange
    modelpar.GG             = JBB_val;   % blue‑blue
    modelpar.pol_val        = pol_val;
    modelpar.cell_pop_prop  = cell_pop_prop;
    
    % (optional) save current parameters in the same .mat file
    save(fullfile(mat_dir,'ParamSweep_1.mat'),...
         'JOO_val','JBO_val','JBB_val','pol_val','cell_pop_prop','-append');
    
    % folder name
    folder_name   = sprintf('JBB%.2f_JOO%.3f_JBO%.2f',JBB_val,JOO_val,JBO_val);
    output_folder = fullfile(base_output_dir,folder_name);
    if ~exist(output_folder,'dir'), mkdir(output_folder); end
    
    % run the simulation
    fprintf('Running %s ...\n',folder_name);
    model(modelpar);                       % <‑‑ Bhaskar simulation
    
    % move output
    if exist('ParamSweep_1_Output','dir')
        movefile('ParamSweep_1_Output',output_folder,'f');
    else
        warning('No ParamSweep_1_Output for %s',folder_name);
        continue;
    end
    
    % (optional) quick thumbnail
    pos_file  = fullfile(output_folder,'ParamSweep_1_Output','Pos_0500000.dat');
    type_file = fullfile(output_folder,'ParamSweep_1_Output','Types_0500000.dat');
    if exist(pos_file,'file') && exist(type_file,'file')
        % --- tiny scatter just for sanity check
        fid  = fopen(pos_file,'r'); raw = fscanf(fid,'%c'); fclose(fid);
        toks = strsplit(raw,',');  pos = cellfun(@str2double,toks);
        X = real(pos);  Y = imag(pos);
        types = load(type_file);
        figure('visible','off');
        scatter(X,Y,6,types,'filled'); axis equal off
        saveas(gcf, fullfile(output_folder,'thumb.png'));
        close(gcf);
    end
end

disp('Line sweep finished.');
