%========================================================
% File: bisectTransition_JOO.m
% Description:
%   Implements the "Initialization" paragraph from Zhao et al.
%   We fix  JBB = 0.05,  JBO = 0.10
%   and localize the transition along the JOO‑axis
%   by bisection on the Wasserstein distance of the
%   area‑distribution feature (orange cells, alpha‑shape).
%========================================================
function [p0, logTable] = bisectTransition_JOO()

% --- constants ----------------------------------------------------------
JBB_val   = 0.05;
JBO_val   = 0.10;

JOO_low   = 0.005;     % initial bracket  (column‑1 & column‑2)
JOO_high  = 0.010;

tol       = 0.0005;    % stop when interval width ≤ tol
maxIter   = 10;        % safety cap

alphaRad  = 1.5;       % alpha‑shape radius
nbins     = 30;        % histogram bins

baseDir   = 'LineSweep_JBB0p05';   % where simulations are / will be stored
% ------------------------------------------------------------------------

iter = 0;
logTable = [];

while (JOO_high - JOO_low) > tol && iter < maxIter
    iter = iter + 1;
    
    % Midpoint
    JOO_mid = 0.5*(JOO_low + JOO_high);
    
    % --- Ensure simulations exist (run if missing) ----------------------
    runIfNeeded(JBB_val, JOO_low,  JBO_val, baseDir);
    runIfNeeded(JBB_val, JOO_mid,  JBO_val, baseDir);
    runIfNeeded(JBB_val, JOO_high, JBO_val, baseDir);
    
    % --- Compute area‑distribution histograms ---------------------------
    h_low  = getAreaHist(JBB_val,JOO_low, JBO_val, baseDir, alphaRad, nbins);
    h_mid  = getAreaHist(JBB_val,JOO_mid, JBO_val, baseDir, alphaRad, nbins);
    h_high = getAreaHist(JBB_val,JOO_high,JBO_val, baseDir, alphaRad, nbins);
    
    % --- Wasserstein distances (1‑D, L1 of CDF diff) --------------------
    W1 = wasserstein1D(h_low , h_mid );
    W2 = wasserstein1D(h_mid , h_high);
    
    % --- Log ------------------------------------------------------------
    logTable = [logTable ; iter JOO_low JOO_mid JOO_high W1 W2];
    fprintf('Iter %d   interval = [%.5f , %.5f]   W1=%.4f  W2=%.4f\n',...
             iter, JOO_low, JOO_high, W1, W2);
    
    % --- Choose half‑interval with larger distance ----------------------
    if W1 > W2
        JOO_high = JOO_mid;   % keep lower half
    else
        JOO_low  = JOO_mid;   % keep upper half
    end
end

p0 = 0.5*(JOO_low + JOO_high);
fprintf('\nLocalized transition at JOO ≈ %.5f (width %.5f)\n',...
        p0, JOO_high-JOO_low);

% convert log to table for convenience
logTable = array2table(logTable,...
    'VariableNames',{'iter','JOO_low','JOO_mid','JOO_high','W_low_mid','W_mid_high'});
end
%=========================================================================


%------------------------------------------------------------------
% Helper: run simulation if folder / Pos_*.dat missing
%------------------------------------------------------------------
function runIfNeeded(JBB,JOO,JBO,baseDir)
    folder = sprintf('JBB%.2f_JOO%.5f_JBO%.2f',JBB,JOO,JBO);
    outDir = fullfile(baseDir,folder,'ParamSweep_1_Output');
    posTag = '0500000';                          % your chosen iteration tag
    if exist(fullfile(outDir,['Pos_' posTag '.dat']),'file')
        return          % already done
    end
    % otherwise build modelpar and run once
    fprintf('  -> running simulation for JOO = %.5f ...\n',JOO);
    modelpar.model          = 'Bhaskar';
    modelpar.RR             = JOO;
    modelpar.RG             = JBO;
    modelpar.GG             = JBB;
    modelpar.pol_val        = 0.005;
    modelpar.cell_pop_prop  = [0.4,0.6];
    
    model(modelpar);
    if ~exist(folder,'dir'), mkdir(fullfile(baseDir,folder)); end
    movefile('ParamSweep_1_Output', fullfile(baseDir,folder));
end

%------------------------------------------------------------------
% Helper: load / compute histogram of areas
%------------------------------------------------------------------
function histProb = getAreaHist(JBB,JOO,JBO,baseDir,alphaRad,nbins)
    folder = sprintf('JBB%.2f_JOO%.5f_JBO%.2f',JBB,JOO,JBO);
    dataDir = fullfile(baseDir,folder,'ParamSweep_1_Output');
    posFile = dir(fullfile(dataDir,'Pos_*.dat'));   posFile = fullfile(dataDir,posFile(1).name);
    typFile = dir(fullfile(dataDir,'Types_*.dat')); typFile = fullfile(dataDir,typFile(1).name);
    
    % load complex positions
    txt  = fileread(posFile);   toks = strsplit(txt,',');
    posC = str2double(toks).';  X = real(posC);  Y = imag(posC);
    types = load(typFile);
    
    keep = (types == 1);        % orange cells (type 1)
    X = X(keep);  Y = Y(keep);
    
    shp = alphaShape(X,Y,alphaRad);
    nReg = numRegions(shp);
    A = zeros(nReg,1);
    for j=1:nReg
        A(j) = area(shp,j);
    end
   
    edges = linspace(0,max(A)+eps,nbins+1);
    histProb = histcounts(A,edges,'Normalization','probability');
end

%------------------------------------------------------------------
% 1‑D 1‑Wasserstein distance between two probability histograms
%------------------------------------------------------------------
function W = wasserstein1D(p,q)
    W = sum( abs( cumsum(p) - cumsum(q) ) );
end
