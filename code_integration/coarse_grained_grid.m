%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: coarse_grained_grid.m
% Description:
%   Loads final-snapshot data from the Bhaskar model parameter sweep,
%   then plots a 5×5 mosaic (coarse-grained grid) of subplots showing 
%   each final pattern for JBB, JOO ∈ {0.00, 0.05, 0.10, 0.15, 0.20},
%   with JBO=0.10 held constant.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%--- 1) Define the parameter sets for JBB and JOO:
JBBvals = [0.00, 0.05, 0.10, 0.15, 0.20];
JOOvals = [0.00, 0.05, 0.10, 0.15, 0.20];
JBO_fixed = 0.10;  % We'll hold JBO=0.10

%--- 2) Define the base directory containing subfolders:
base_dir = 'ParamSweep_Results';

%--- 3) Create a figure sized to have 5×5 subplots
figure('Name','Coarse-Grained Diagram',...
       'Units','normalized','Position',[0.1 0.1 0.80 0.80],...
       'Color','w');

nB = length(JBBvals);
nO = length(JOOvals);

%--- 4) Nested loop over JBB, JOO to produce 5×5 = 25 subplots
plotIndex = 1;
for iB = 1:nB
    for iO = 1:nO
        
        JBB_val = JBBvals(iB);
        JOO_val = JOOvals(iO);
        
        % Build folder name, e.g.: 'JBB0.00_JOO0.10_JBO0.10'
        folderName = sprintf('JBB%.2f_JOO%.2f_JBO%.2f',...
                             JBB_val, JOO_val, JBO_fixed);
        output_folder = fullfile(base_dir, folderName, 'ParamSweep_1_Output');
        
        %--- 4a) Make a subplot for this (JBB,JOO)
        subplot(nB, nO, plotIndex);
        plotIndex = plotIndex + 1;
        
        %--- 4b) Construct paths to final .dat files
        pos_file  = fullfile(output_folder, 'Pos_0500000.dat');
        type_file = fullfile(output_folder, 'Types_0500000.dat');
        
        % If missing, skip or show empty subplot
        if ~isfile(pos_file) || ~isfile(type_file)
            title(sprintf('Missing data:\n%s', folderName),'Interpreter','none','FontSize',8);
            axis off
            continue;
        end
        
        %--- 4c) Load the final position data (comma-delimited single-column complex)
        fid = fopen(pos_file,'r');
        rawText = fscanf(fid, '%c');
        fclose(fid);
        
        tokens = strsplit(rawText, ',');
        N = numel(tokens);
        posData = zeros(N,1);
        for kToken = 1:N
            posData(kToken) = str2double(tokens{kToken});
        end
        
        %--- 4d) Load type data
        typeData = load(type_file);  % Nx1 (0 or 1)
        
        if length(typeData) ~= length(posData)
            title(sprintf('Mismatch #positions vs #types:\n%s', folderName), ...
                'Interpreter','none','FontSize',8);
            axis off
            continue;
        end
        
        % Convert complex to x,y
        X = real(posData);
        Y = imag(posData);
        
        %--- 4e) Scatter
        scatter(X, Y, 7, typeData, 'filled');
        axis equal; 
        xlim([-10,10]); ylim([-10,10]);  % adapt if domain is bigger
        title(folderName, 'Interpreter','none','FontSize',8);
        
    end
end

%--- 5) Overall label
sgtitle('Coarse-Grained Grid: JBO=0.10, JBB & JOO in [0,0.2]','FontSize',14);
