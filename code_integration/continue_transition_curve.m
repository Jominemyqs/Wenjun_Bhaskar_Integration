function continue_transition_curve()
    %% Setup
    baseDir   = 'LineSweep_JBB0p05';
    posTag    = 'JBO0.10';
    alphaRad  = 1.5;
    nbins     = 30;
    h         = 0.002;
    s         = 0.005;       % arclength step
    maxPts    = 4;           % limit to 2 steps for debug

    figFile   = fullfile(baseDir, 'transition_curve.png');
    saveFile  = fullfile(baseDir, 'pts.mat');

    %% Step 0: Try to load or create pts
    if isfile(saveFile)
        load(saveFile, 'pts');
        fprintf('Loaded existing pts with %d points.\n', size(pts,1));
    else
        fprintf('Starting from scratch with initial two points.\n');
        pts = [ 0.05, 0.01016 ;
                0.05, 0.03621 ];
        save(saveFile, 'pts');
    end

    %% Plot current state
    figure(10); clf;
    plot(pts(:,1), pts(:,2), 'ko-');
    xlabel('J_{BB}'); ylabel('J_{OO}');
    title('Transition Curve (J_{BB}, J_{OO})');
    saveas(gcf, figFile);
    drawnow;

    %% Start continuation
    while size(pts, 1) < maxPts
        stepCount = size(pts, 1);
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

        % Append, save, plot
        pts(end+1,:) = corr;
        save(saveFile, 'pts');

       figure(10); clf;
set(gcf, 'Color', 'w');       % Set figure background to white
set(gca, 'Color', 'w');       % Set axes background to white

plot(pts(:,1), pts(:,2), 'ko-');

xlabel('J_{BB}'); ylabel('J_{OO}');
title('Transition Curve (J_{BB}, J_{OO})');

saveas(gcf, figFile);
drawnow;

        % Stop if out of bounds
        if any(corr < 0) || any(corr > 0.1)
            fprintf('Out of bounds — terminating.\n');
            break;
        end
    end
end

%=========================================================================

function g = bif_g(p, h, posTag, baseDir, alphaRad, nbins)
    p_plus  = p + h * [1 0];
    p_minus = p - h * [1 0];

    runIfNeeded(p_plus(1), p_plus(2), 0.10, baseDir, posTag);
    runIfNeeded(p_minus(1), p_minus(2), 0.10, baseDir, posTag);

    h1 = getRoundnessHist(p_plus(1), p_plus(2), 0.10, baseDir, alphaRad, nbins, posTag);
    h2 = getRoundnessHist(p_minus(1), p_minus(2), 0.10, baseDir, alphaRad, nbins, posTag);

    g = sum(abs(cumsum(h1) - cumsum(h2)));
end

%=========================================================================

function runIfNeeded(JBB, JOO, JBO, baseDir, posTag)
    folder  = sprintf('JBB%.2f_JOO%.5f_JBO%.2f', JBB, JOO, JBO);
    outDir  = fullfile(baseDir, folder, 'ParamSweep_1_Output');
    posFile = fullfile(outDir, ['Pos_' posTag '.dat']);

    if exist(posFile, 'file')
        return;
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

    if ~exist(caseDir, 'dir'); mkdir(caseDir); end
    if exist(outDir, 'dir'); rmdir(outDir, 's'); end
    if ~exist(fresh, 'dir'); error('Solver did not create folder "%s".', fresh); end

    movefile(fresh, caseDir);
end

%=========================================================================

function histProb = getRoundnessHist(JBB, JOO, JBO, baseDir, alphaRad, nbins, posTag)
    folder  = sprintf('JBB%.2f_JOO%.5f_JBO%.2f', JBB, JOO, JBO);
    dataDir = fullfile(baseDir, folder, 'ParamSweep_1_Output');

    posFile = fullfile(dataDir, ['Pos_' posTag '.dat']);
    typFile = fullfile(dataDir, ['Types_' posTag '.dat']);

    if ~exist(posFile, 'file') || ~exist(typFile, 'file')
        posFile = fullfile(dataDir, 'Pos_0500000.dat');
        typFile = fullfile(dataDir, 'Types_0500000.dat');
        if ~exist(posFile, 'file') || ~exist(typFile, 'file')
            error('No position/type file found in %s', dataDir);
        end
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
