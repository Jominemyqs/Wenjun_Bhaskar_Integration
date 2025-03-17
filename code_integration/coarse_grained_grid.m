%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: coarse_grained_grid.m.m
% Description:
%   Loads final-snapshot data from the Bhaskar model parameter sweep,
%   then plots a mosaic (coarse-grained grid) of subplots showing 
%   each param set’s final pattern. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; close all; clc;

%--- 1) Define your param-set subfolders
set_folders = {
    'JBB0.05_JOO0.05_JBO0.05',...
    'JBB0.23_JOO0.05_JBO0.05',...
    'JBB0.23_JOO0.23_JBO0.13'
};

%--- 2) Define the base directory containing those subfolders.
base_dir = 'ParamSweep_Results';  
% e.g. "code_integration/ParamSweep_Results" if needed

%--- 3) Create a figure
figure('Name','Coarse-Grained Diagram',...
       'Units','normalized','Position',[0.1 0.1 0.85 0.45],...
       'Color','w');

%--- 4) Loop over the subfolders
n = length(set_folders);
for idx = 1:n
    
    subplot(1, n, idx);
    
    % Path to final .dat files
    this_folder   = set_folders{idx};
    output_folder = fullfile(base_dir, this_folder, 'ParamSweep_1_Output');
    pos_file      = fullfile(output_folder, 'Pos_5000000.dat');
    type_file     = fullfile(output_folder, 'Types_5000000.dat');
    
    % If either file is missing, skip or show empty subplot
    if ~exist(pos_file,'file') || ~exist(type_file,'file')
        title(sprintf('Missing data:\n%s', this_folder),'Interpreter','none');
        axis off
        continue;
    end
    
    %--- 4a) Manually parse Pos_5000000.dat, since it's comma-delimited complex
    fid = fopen(pos_file,'r');
    rawText = fscanf(fid, '%c');    % read entire file into one string
    fclose(fid);
    
    % Split by commas => e.g. "0.58+0.69i"
    tokens = strsplit(rawText, ',');
    N = numel(tokens);
    
    % Convert each token to complex
    posData = zeros(N,1);
    for kToken = 1:N
        posData(kToken) = str2double(tokens{kToken});
    end
    
    % If you had Nx2 real instead, you'd do: posData = load(pos_file);
    % but here we have Nx1 complex in comma-delimited text.
    
    %--- 4b) Load types (assuming normal numeric row-by-row)
    type_data = load(type_file);   % Nx1 ints (0 or 1)
    
    %--- 4c) Extract X, Y
    X = real(posData);
    Y = imag(posData);
    
    % If type_data doesn't match length, skip
    if length(type_data) ~= length(X)
        title(sprintf('Mismatch #positions vs #types:\n%s', this_folder), ...
              'Interpreter','none');
        axis off
        continue;
    end
    
    %--- 4d) Make a scatter plot
    scatter(X, Y, 7, type_data, 'filled');  % color-coded by type
    axis equal;
    
    % Adjust domain if needed. Example here: [-10, 10] in both directions:
    xlim([-10, 10]); ylim([-10, 10]);
    
    title(strrep(this_folder,'_','\_'), 'FontSize',11, 'Interpreter','none');
end

%--- 5) Overall title:
sgtitle('Coarse-Grained Bifurcation Diagram','FontSize',14);