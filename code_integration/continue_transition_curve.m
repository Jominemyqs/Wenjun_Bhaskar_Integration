%===========================================================
% File : continue_transition_curve.m (updated version)
% Goal : Continue the transition curve and checkpoint progress.
%===========================================================

function continue_transition_curve()
    %% Setup
    baseDir   = 'LineSweep_JBB0p05';
    posTag    = 'JBO0.10';
    alphaRad  = 1.5;
    nbins     = 30;
    h         = 0.002;
    s         = 0.005;       % arclength step
    maxSteps  = 200;         % safety cap

    figFile   = fullfile(baseDir, 'transition_curve.png');
    saveFile  = fullfile(baseDir, 'pts.mat');
    
    % Set a total cap on number of points to generate
    maxPts = 200;


    %% Step 0: Try to load pts.mat
    if isfile(saveFile)
        load(saveFile, 'pts');
        fprintf('Loaded existing pts with %d points.\n', size(pts,1));
    else
        %% Scan folders to reconstruct pts
        folders = dir(fullfile(baseDir, 'JBB*_JOO*_JBO*.10'));
        pts = [];
        for k = 1:length(folders)
            tokens = regexp(folders(k).name, ...
                'JBB([\d\.]+)_JOO([\d\.]+)_JBO', 'tokens');
            if ~isempty(tokens)
                nums = str2double(tokens{1});
                pts(end+1,:) = nums; %#ok<AGROW>
            end
        end
        pts = sortrows(unique(pts, 'rows'), [1 2]);
        fprintf('Scanned %d folders, recovered %d unique points.\n', ...
            length(folders), size(pts,1));

        % If still empty, fallback to hardcoded points
        if size(pts,1) < 2
            fprintf('Starting from scratch with initial two points.\n');
            pts = [ 0.05, 0.01016 ;
                    0.05, 0.03621 ];
        end
        save(saveFile, 'pts');
        % Immediately plot and save the recovered points
        figure(10); clf;
        plot(pts(:,1), pts(:,2), 'ko-');
        xlabel('J_{BB}'); ylabel('J_{OO}');
        title('Transition Curve (J_{BB}, J_{OO})');
        saveas(gcf, figFile);
        drawnow;

    end

    %% Start continuation
    while size(pts, 1) < maxPts
    stepCount = size(pts, 1);  % dynamic counter
    fprintf('Step %d:\n', stepCount);

        fprintf('Step %d:\n', stepCount);

        % Compute secant
        t = pts(end,:) - pts(end-1,:);
        tangent = t / norm(t);

        % Predictor
        pred = pts(end,:) + s * tangent;

        % Corrector
        normal = [-tangent(2), tangent(1)];
        gfun = @(z) bif_g(pred + z*normal, h, posTag, baseDir, alphaRad, nbins);
        [zopt, ~] = fminbnd(@(z) -gfun(z), -0.05, 0.05);
        corr = pred + zopt * normal;

        % Append, save, update
        pts(end+1,:) = corr;
        save(saveFile, 'pts');

        % Plot update
        figure(10); clf;
        plot(pts(:,1), pts(:,2), 'ko-');
        xlabel('J_{BB}'); ylabel('J_{OO}');
        title('Transition Curve (J_{BB}, J_{OO})');
        saveas(gcf, figFile);
        drawnow;

        % Safety: stop if out of bounds
        if any(corr < 0) || any(corr > 0.1)
            fprintf('Out of bounds — terminating.\n');
            break;
        end
    end
end


%======================  END MAIN =======================================


%-----------------------------------------------------------------------
% bif_g  – computes Wasserstein distance feature between (p+h) and (p-h)
%-----------------------------------------------------------------------
function g = bif_g(p, h, posTag, baseDir, alphaRad, nbins)

p_plus  = p + h * [1 0];
p_minus = p - h * [1 0];

% ensure both simulations exist
runIfNeeded(p_plus(1), p_plus(2), 0.10, baseDir, posTag);
runIfNeeded(p_minus(1), p_minus(2), 0.10, baseDir, posTag);

% compute histograms
h1 = getRoundnessHist(p_plus(1), p_plus(2), 0.10, baseDir, alphaRad, nbins, posTag);
h2 = getRoundnessHist(p_minus(1), p_minus(2), 0.10, baseDir, alphaRad, nbins, posTag);

% compute Wasserstein distance
g = sum(abs(cumsum(h1) - cumsum(h2)));
end

%-----------------------------------------------------------------------
% runIfNeeded  – run simulation if Pos_ file is missing
%-----------------------------------------------------------------------
function runIfNeeded(JBB, JOO, JBO, baseDir, posTag)

folder  = sprintf('JBB%.2f_JOO%.5f_JBO%.2f', JBB, JOO, JBO);
outDir  = fullfile(baseDir, folder, 'ParamSweep_1_Output');
posFile = fullfile(outDir, ['Pos_' posTag '.dat']);

if exist(posFile, 'file')
    return; % already done
end

fprintf('  -> running simulation for  JBB=%.5f  JOO=%.5f\n', JBB, JOO);

modelpar              = struct();
modelpar.model        = 'Bhaskar';
modelpar.sets         = 'Bhaskar';
modelpar.GG           = JBB;
modelpar.RR           = JOO;
modelpar.RG           = JBO;
modelpar.pol_val      = 0.005;
modelpar.cell_pop_prop= [0.60 0.40];
modelpar.iterSave     = str2double(posTag);

model(modelpar);

fresh = 'ParamSweep_1_Output';
caseDir = fullfile(baseDir, folder);

if ~exist(caseDir, 'dir')
    mkdir(caseDir);
end

if exist(outDir, 'dir')
    rmdir(outDir, 's');
end

if ~exist(fresh, 'dir')
    error('Solver did not create folder "%s".', fresh);
end

movefile(fresh, caseDir);
end

%-----------------------------------------------------------------------
% getRoundnessHist  – roundness-score histogram for orange cells
%-----------------------------------------------------------------------
function histProb = getRoundnessHist(JBB, JOO, JBO, baseDir, alphaRad, nbins, posTag)
folder  = sprintf('JBB%.2f_JOO%.5f_JBO%.2f', JBB, JOO, JBO);
dataDir = fullfile(baseDir, folder, 'ParamSweep_1_Output');

posFile = fullfile(dataDir, ['Pos_' posTag '.dat']);
typFile = fullfile(dataDir, ['Types_' posTag '.dat']);

if ~exist(posFile, 'file') || ~exist(typFile, 'file')
    error('files not found in %s', dataDir);
end

txt = fileread(posFile);  toks = strsplit(txt, ',');
posC = str2double(toks).';  X = real(posC);  Y = imag(posC);
types = load(typFile);

orange = (types == 1);
X = X(orange); Y = Y(orange);

if numel(X) < 3
    histProb = zeros(1, nbins);
    return;
end

shp = alphaShape(X, Y, alphaRad);
nr  = numRegions(shp);
roundness = zeros(nr,1);

for k = 1:nr
    A = area(shp,k);
    P = perimeter(shp,k);
    roundness(k) = min(4*pi*A / P^2, 1);
end

edges = linspace(0,1,nbins+1);
histProb = histcounts(roundness, edges, 'Normalization', 'probability');
end