%===========================================================
% File : bisectTransition_JOO_roundness.m
% Goal : Localise the JOO transition at fixed JBB = 0.05, JBO = 0.10
%        using the roundness‑score feature for ORANGE cells.
%
% 
%===========================================================
function [p0, logTable] = bisectTransition_JOO_roundness()

% ---------- fixed parameters ------------------------------------------
JBB_val   = 0.05;
JBO_val   = 0.10;

JOO_low   = 0.020;          % <-- left end of initial bracket
JOO_high  = 0.045;          % <-- right end (change as desired)

tol       = 0.0005;         % stop when width ≤ tol
maxIter   = 12;             % safety cap

alphaRad  = 1.5;            % α‑shape radius for roundness
nbins     = 30;             % histogram bins

baseDir   = 'LineSweep_JBB0p05';   % parent folder that stores runs
posTag    = '0500000';             % iteration tag inside Pos_*.dat
% ----------------------------------------------------------------------

iter = 0;
logTable = [];

while (JOO_high - JOO_low) > tol && iter < maxIter
    iter = iter + 1;
    JOO_mid = 0.5*(JOO_low + JOO_high);

    % --- make sure the three simulations exist ------------------------
    runIfNeeded(JBB_val,JOO_low ,JBO_val,baseDir,posTag);
    runIfNeeded(JBB_val,JOO_mid ,JBO_val,baseDir,posTag);
    runIfNeeded(JBB_val,JOO_high,JBO_val,baseDir,posTag);

    % --- compute roundness histograms ---------------------------------
    h_low  = getRoundnessHist(JBB_val,JOO_low ,JBO_val,baseDir,alphaRad,nbins,posTag);
    h_mid  = getRoundnessHist(JBB_val,JOO_mid ,JBO_val,baseDir,alphaRad,nbins,posTag);
    h_high = getRoundnessHist(JBB_val,JOO_high,JBO_val,baseDir,alphaRad,nbins,posTag);

    % --- 1‑Wasserstein distances --------------------------------------
    W_low_mid  = wasserstein1D(h_low ,h_mid );
    W_mid_high = wasserstein1D(h_mid ,h_high);

    logTable = [logTable ; iter JOO_low JOO_mid JOO_high W_low_mid W_mid_high];
    fprintf('Iter %2d  interval=[%.5f , %.5f]  W1=%.3f  W2=%.3f\n',...
            iter, JOO_low, JOO_high, W_low_mid, W_mid_high);

    % --- bisection decision -------------------------------------------
    if W_low_mid > W_mid_high
        JOO_high = JOO_mid;     % keep left half
    else
        JOO_low  = JOO_mid;     % keep right half
    end
end

p0 = 0.5*(JOO_low + JOO_high);
fprintf('\nEstimated transition at  JOO ≈ %.5f   (final width %.5f)\n',...
        p0, JOO_high-JOO_low);

% return a nice table
logTable = array2table(logTable,...
    'VariableNames',{'iter','JOO_low','JOO_mid','JOO_high','W_low_mid','W_mid_high'});
end
%======================  END MAIN =======================================



%-----------------------------------------------------------------------
% runIfNeeded  – runs Bhaskar model once if the Pos_ file is missing
%-----------------------------------------------------------------------
function runIfNeeded(JBB,JOO,JBO,baseDir,posTag)

folder  = sprintf('JBB%.2f_JOO%.5f_JBO%.2f',JBB,JOO,JBO);
outDir  = fullfile(baseDir,folder,'ParamSweep_1_Output');
posFile = fullfile(outDir,['Pos_' posTag '.dat']);

if exist(posFile,'file');  return;  end     % nothing to do

fprintf('  -> running simulation for  JOO = %.5f\n',JOO);

% ---------- build the structure expected by model.m -------------------
modelpar              = struct();
modelpar.model        = 'Bhaskar';   % << mandatory
modelpar.sets         = 'Bhaskar';   % make model.m return early

modelpar.GG           = JBB;           % blue–blue
modelpar.RR           = JOO;           % orange–orange
modelpar.RG           = JBO;           % blue–orange

modelpar.pol_val      = 0.005;         % your usual polarity
modelpar.cell_pop_prop= [0.60 0.40];   % <- two entries, sum = 1


% save frequency so the desired Pos_ file is written
modelpar.iterSave     = str2double(posTag);   % 500 000 for “0500000”

% ---------- call the solver -------------------------------------------
model(modelpar);

fresh = 'ParamSweep_1_Output';                       % name written by solver
caseDir = fullfile(baseDir,folder);                  % e.g. …\JBB0.05_JOO0.01250_JBO0.10
outDir  = fullfile(caseDir,'ParamSweep_1_Output');   % desired final path

if ~exist(caseDir,'dir');  mkdir(caseDir);  end      % make parent

if exist(outDir,'dir')                               % stale result?
    rmdir(outDir,'s');                               % delete it (s = recursive)
end

if ~exist(fresh,'dir')
    error('Solver did not create folder "%s".',fresh);
end
movefile(fresh,caseDir);                             % move once, no nesting
end


%-----------------------------------------------------------------------
% getRoundnessHist  – roundness‑score histogram for orange cells
%-----------------------------------------------------------------------
function histProb = getRoundnessHist(JBB,JOO,JBO,baseDir,alphaRad,nbins,posTag)
folder  = sprintf('JBB%.2f_JOO%.5f_JBO%.2f',JBB,JOO,JBO);
dataDir = fullfile(baseDir,folder,'ParamSweep_1_Output');

% --- load particle positions ------------------------------------------
posFile = fullfile(dataDir,['Pos_' posTag '.dat']);
typFile = fullfile(dataDir,['Types_' posTag '.dat']);
if ~exist(posFile,'file') || ~exist(typFile,'file')
    error('files not found in %s',dataDir);
end

txt  = fileread(posFile);  toks = strsplit(txt,',');
posC = str2double(toks).';  X = real(posC);  Y = imag(posC);
types = load(typFile);

orange = (types == 1);          % type‑1 = orange
X = X(orange);  Y = Y(orange);

if numel(X) < 3      % too few points
    histProb = zeros(1,nbins);  return;
end

% --- α‑shape and roundness scores --------------------------------------
shp = alphaShape(X,Y,alphaRad);
nr  = numRegions(shp);
roundness = zeros(nr,1);
for k = 1:nr
    A = area     (shp,k);
    P = perimeter(shp,k);
    roundness(k) = min(4*pi*A / P^2 , 1);   % 1 = perfect disk
end

edges = linspace(0,1,nbins+1);
histProb = histcounts(roundness,edges,'Normalization','probability');
end



%-----------------------------------------------------------------------
% wasserstein1D  – 1‑Wasserstein between two 1‑D histograms
%-----------------------------------------------------------------------
function W = wasserstein1D(p,q)
W = sum( abs( cumsum(p) - cumsum(q) ) );
end
%=======================================================================


%====================  OPTIONAL: use "bag of areas"  ===================
% If you prefer the area distribution instead of roundness,
%  * rename getRoundnessHist to getAreaBagHist
%  * inside the loop use that function
%
% getAreaBagHist is identical except you replace the roundness
% computation by:
%
%     nr  = numRegions(shp);
%     A   = zeros(nr,1);
%     for k = 1:nr
%         A(k) = area(shp,k);
%     end
%     edges = linspace(0,max(A)+eps,nbins+1);
%     histProb = histcounts(A,edges,'Normalization','probability');
%
%=======================================================================
